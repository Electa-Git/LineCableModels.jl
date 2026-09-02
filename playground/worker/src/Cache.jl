"""
    CachedResult

Store one normalized deterministic result and its scientific warnings.
"""
struct CachedResult
    result::Dict{String,Any}
    warnings::Vector{String}
end

"""
    ResultCache

Store immutable completed results as content-addressed JSON files.
"""
mutable struct ResultCache
    directory::String
    lock::ReentrantLock
end

function ResultCache(directory::AbstractString)
    mkpath(directory)
    return ResultCache(abspath(directory), ReentrantLock())
end

function result_cache_key(
        request::JobRequest,
        schema_version::AbstractString,
        engine_version::AbstractString,
        environment_fingerprint::AbstractString;
        parameters=request.parameters
    )
    scientific_input_hash = input_hash(
        request.operation,
        parameters;
        engine_constraint=request.engine_constraint,
        schema_version
    )
    return bytes2hex(SHA.sha256(join((
        request.operation,
        PROTOCOL_VERSION,
        string(schema_version),
        string(engine_version),
        string(environment_fingerprint),
        scientific_input_hash,
    ), '\0')))
end

result_cache_path(cache::ResultCache, key::AbstractString) =
    joinpath(cache.directory, "$key.json")

function cache_get(cache::ResultCache, key::AbstractString)
    path = result_cache_path(cache, key)
    return lock(cache.lock) do
        isfile(path) || return nothing
        document = JSON3.read(read(path, String))
        CachedResult(
            normalize_wire(document.result),
            String[string(item) for item in document.warnings]
        )
    end
end

function cache_put!(
        cache::ResultCache,
        key::AbstractString,
        value::CachedResult
    )
    path = result_cache_path(cache, key)
    lock(cache.lock) do
        isfile(path) && return path
        temporary, io = mktemp(cache.directory)
        try
            write(io, JSON3.write(Dict(
                "result" => value.result,
                "warnings" => value.warnings,
            )))
            flush(io)
            close(io)
            mv(temporary, path)
        catch
            isopen(io) && close(io)
            isfile(temporary) && rm(temporary; force=true)
            rethrow()
        end
    end
    return path
end

mutable struct PreparedEntry
    state::Symbol
    value::Any
    task::Union{Nothing,Task}
    expires_at::Float64
    last_error::Union{Nothing,Exception}
end

"""
    PreparedResourceCache(; ttl_seconds=900.0, max_entries=8)

Hold bounded worker-local scientific preparations. Entries are leased for
`ttl_seconds` \\[s\\] and concurrent requests for one key share a single task.
"""
mutable struct PreparedResourceCache
    entries::Dict{String,PreparedEntry}
    ttl_seconds::Float64
    max_entries::Int
    lock::ReentrantLock
end

function PreparedResourceCache(; ttl_seconds::Real=900.0, max_entries::Integer=8)
    ttl_seconds > 0 || throw(ArgumentError("prepared-resource TTL must be positive"))
    max_entries > 0 || throw(ArgumentError("prepared-resource limit must be positive"))
    return PreparedResourceCache(
        Dict{String,PreparedEntry}(),
        Float64(ttl_seconds),
        Int(max_entries),
        ReentrantLock()
    )
end

function prune_prepared!(cache::PreparedResourceCache)
    now = time()
    lock(cache.lock) do
        filter!(pair -> last(pair).state == :warming || last(pair).expires_at > now,
            cache.entries)
        while length(cache.entries) > cache.max_entries
            candidates = filter(pair -> last(pair).state != :warming,
                collect(cache.entries))
            isempty(candidates) && break
            victim = first(sort!(candidates; by=pair -> last(pair).expires_at))
            delete!(cache.entries, first(victim))
        end
    end
    return cache
end

"""
    prepare_resource!(builder, cache, key)

Return one worker-local prepared resource. Concurrent callers for the same key
share one construction task; completed entries expire after the configured
lease interval.
"""
function prepare_resource!(builder, cache::PreparedResourceCache, key::AbstractString)
    prune_prepared!(cache)
    task = lock(cache.lock) do
        entry = get(cache.entries, string(key), nothing)
        if !isnothing(entry) && entry.state == :hot && entry.expires_at > time()
            entry.expires_at = time() + cache.ttl_seconds
            return entry.value
        elseif !isnothing(entry) && entry.state == :warming
            return entry.task
        elseif !isnothing(entry) && entry.state == :failed &&
                entry.expires_at > time()
            throw(something(entry.last_error, ErrorException(
                "prepared-resource construction failed"
            )))
        end

        resource_task = @async builder()
        cache.entries[string(key)] = PreparedEntry(
            :warming,
            nothing,
            resource_task,
            Inf,
            nothing
        )
        return resource_task
    end
    task isa Task || return task
    try
        value = fetch(task)
        lock(cache.lock) do
            cache.entries[string(key)] = PreparedEntry(
                :hot,
                value,
                nothing,
                time() + cache.ttl_seconds,
                nothing
            )
        end
        prune_prepared!(cache)
        return value
    catch error
        lock(cache.lock) do
            cache.entries[string(key)] = PreparedEntry(
                :failed,
                nothing,
                nothing,
                time() + min(30.0, cache.ttl_seconds),
                error
            )
        end
        rethrow()
    end
end

function prepared_status(cache::PreparedResourceCache, key::AbstractString)
    prune_prepared!(cache)
    return lock(cache.lock) do
        entry = get(cache.entries, string(key), nothing)
        isnothing(entry) ? :cold : entry.state
    end
end
