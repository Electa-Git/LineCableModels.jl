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
    "Estimated retained payload size in bytes; excludes transient solver memory."
    retained_bytes::Int
end

"""
    PreparedResourceCache(; ttl_seconds=900, max_entries=8,
                            max_bytes=268435456, clock=() -> time_ns() / 1e9)

Keep single-flight preparations within one executor. Retain at most
`max_entries` keys and `max_bytes` estimated bytes of completed values.
`ttl_seconds` \\[s\\] expires idle entries using a monotonic clock.

The byte limit bounds retained cache values, not transient solver allocation,
package memory or peak process memory. The process supervisor owns those
resource and deadline limits. A full cache of in-progress builders rejects
new keys instead of creating an unbounded queue.
"""
mutable struct PreparedResourceCache
    "Entries owned by this executor."
    entries::Dict{String,PreparedEntry}
    "Idle expiry interval \\[s\\]."
    ttl_seconds::Float64
    "Maximum resident and in-progress keys."
    max_entries::Int
    "Maximum estimated retained value bytes."
    max_bytes::Int
    "Monotonic local time source."
    clock::Function
    "Serialize admission, completion and eviction."
    lock::ReentrantLock
end

function PreparedResourceCache(; ttl_seconds::Real=900.0, max_entries::Integer=8,
        max_bytes::Integer=256 * 1024^2, clock::Function=() -> time_ns() / 1e9)
    !(ttl_seconds isa Bool) && isfinite(ttl_seconds) && 0 < ttl_seconds <= 86400 ||
        throw(ArgumentError("prepared-resource TTL must be in (0, 86400] seconds"))
    !(max_entries isa Bool) && 1 <= max_entries <= 256 ||
        throw(ArgumentError("prepared-resource limit must be in 1:256"))
    !(max_bytes isa Bool) && 0 < max_bytes <= typemax(Int) ||
        throw(ArgumentError("prepared-resource byte limit must be positive"))
    return PreparedResourceCache(Dict{String,PreparedEntry}(), Float64(ttl_seconds),
        Int(max_entries), Int(max_bytes), clock, ReentrantLock())
end

function prune_prepared!(cache::PreparedResourceCache)
    lock(cache.lock) do
        now = cache.clock()
        filter!(pair -> last(pair).state == :warming || last(pair).expires_at > now,
            cache.entries)
    end
    return cache
end

# Caller holds the cache lock. Pending work and the current builder cannot be
# evicted to make an apparent free slot; only completed values are replaceable.
function evict_prepared!(cache::PreparedResourceCache; except=nothing)
    candidates = [pair for pair in cache.entries
        if last(pair).state != :warming && last(pair) !== except]
    isempty(candidates) && return false
    victim = first(sort!(candidates; by=pair -> (last(pair).expires_at, first(pair))))
    delete!(cache.entries, first(victim))
    return true
end

"""
    prepare_resource!(builder, cache, key)

Return a prepared value, sharing one admitted builder per key. Successful values
renew their idle TTL on reuse. Failures have a short retry delay; cache eviction
or explicit cleanup may remove that history. Admission and retained-byte checks
are serialized. Completion after cleanup cannot repopulate the cache.
"""
function prepare_resource!(builder, cache::PreparedResourceCache, key::AbstractString)
    id = String(key)
    0 < ncodeunits(id) <= 256 || throw(ArgumentError("prepared-resource key must contain 1–256 bytes"))
    choice = lock(cache.lock) do
        prune_prepared!(cache)
        entry = get(cache.entries, id, nothing)
        if entry !== nothing
            if entry.state == :hot
                entry.expires_at = cache.clock() + cache.ttl_seconds
                return (:value, entry.value)
            elseif entry.state == :warming
                return (:task, entry.task)
            elseif entry.state == :failed
                throw(something(entry.last_error, ErrorException("prepared-resource construction failed")))
            end
        end
        while length(cache.entries) >= cache.max_entries
            evict_prepared!(cache) || throw(ArgumentError("prepared-resource capacity is occupied"))
        end
        entry = PreparedEntry(:warming, nothing, nothing, Inf, nothing, 0)
        cache.entries[id] = entry
        entry.task = @async begin
            try
                value = builder()
                bytes = Base.summarysize(value)
                bytes <= cache.max_bytes || throw(ArgumentError("prepared resource exceeds retained-byte limit"))
                lock(cache.lock) do
                    get(cache.entries, id, nothing) === entry ||
                        throw(OperationCanceled())
                    prune_prepared!(cache)
                    while sum(e.retained_bytes for e in values(cache.entries); init=0) + bytes > cache.max_bytes
                        evict_prepared!(cache; except=entry) ||
                            throw(ArgumentError("prepared-resource byte capacity is occupied"))
                    end
                    entry.state = :hot
                    entry.value = value
                    entry.task = nothing
                    entry.retained_bytes = bytes
                    entry.expires_at = cache.clock() + cache.ttl_seconds
                end
                return value
            catch error
                lock(cache.lock) do
                    if get(cache.entries, id, nothing) === entry
                        entry.state = :failed
                        entry.value = nothing
                        entry.task = nothing
                        entry.retained_bytes = 0
                        entry.last_error = error isa Exception ? error : ErrorException("prepared-resource construction failed")
                        entry.expires_at = cache.clock() + min(30.0, cache.ttl_seconds)
                    end
                end
                rethrow()
            end
        end
        return (:task, entry.task)
    end
    return first(choice) == :value ? last(choice) : fetch(last(choice))
end

function prepared_status(cache::PreparedResourceCache, key::AbstractString)
    return lock(cache.lock) do
        prune_prepared!(cache)
        entry = get(cache.entries, string(key), nothing)
        entry === nothing ? :cold : entry.state
    end
end

"""
    clear_prepared!(cache)

Invalidate this executor's retained values and pending cache identities. A builder
that finishes later cannot restore a cleared entry. This does not interrupt
numerical code; cancellation and hard termination remain supervisor responsibilities.
"""
function clear_prepared!(cache::PreparedResourceCache)
    lock(cache.lock) do
        empty!(cache.entries)
    end
    return nothing
end
