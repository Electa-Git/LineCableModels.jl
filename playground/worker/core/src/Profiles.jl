"""
    AbstractScientificProfile

Extend four hooks for an approved numerical environment: `operation_registry`,
`validate_preparation`, `prepare!`, and `cleanup!`. The host agent retains lease,
resource admission, process lifetime, timeout and result-authorization decisions.
"""
abstract type AbstractScientificProfile end
const PreparationParameters = Dict{String,Any}

@required AbstractScientificProfile begin
    operation_registry(::AbstractScientificProfile)
    validate_preparation(::AbstractScientificProfile, ::PreparationParameters)
    prepare!(::AbstractScientificProfile, ::ExecutionContext, ::PreparationParameters)
    cleanup!(::AbstractScientificProfile, ::PreparedResourceCache)
end

"""
    operation_registry(profile) -> OperationRegistry

Return only the operations provided by this numerical environment. Constructing
the registry must not prepare a model or perform numerical work.
"""
function operation_registry end

"""
    load_profile!(profile)

Load the profile's fixed, approved numerical dependencies inside its executor.
The default does nothing. Loading must be idempotent and must not prepare a
model. The command reader calls this only for explicit preparation or execution,
under that operation's deadline, never while announcing bootstrap availability.
"""
load_profile!(::AbstractScientificProfile) = nothing

"""
    validate_preparation(profile, parameters) -> Dict{String,Any}

Validate and normalize preparation inputs. The normalized values determine the
single-flight preparation key within this disposable executor's private cache.
"""
function validate_preparation end

"""
    prepare!(profile, context, parameters) -> PreparedWorkload

Run the declared representative workload using validated inputs and return its
evidence. Package import alone is insufficient. Evidence describes the actual
workload, not a guarantee that every future Julia specialization is compiled.
"""
function prepare! end

"""
    cleanup!(profile, cache)

Release profile-owned preparations on graceful executor exit. The parent still
owns hard termination and scratch cleanup when numerical code cannot cooperate.
"""
function cleanup! end

"""
    PreparedWorkload(evidence; resources=())

Record a completed representative workload and the private cache keys it needs
to remain valid. The common executor checks these dependencies on reuse; expired
or evicted models cannot be represented by a retained ready label.
"""
struct PreparedWorkload
    "Normalized, browser-safe representative-workload evidence."
    evidence::Dict{String,Any}
    "Required model keys in this executor's cache."
    resources::Tuple{Vararg{String}}
    function PreparedWorkload(evidence::Dict{String,Any}; resources=())
        keys = Tuple(String.(resources))
        length(keys) <= 64 && allunique(keys) && all(k -> 0 < ncodeunits(k) <= 256, keys) ||
            throw(ArgumentError("invalid preparation resource keys"))
        return new(normalize_wire(evidence), keys)
    end
end

function preparation_dependencies!(cache::PreparedResourceCache, value::PreparedWorkload)
    now = cache.clock()
    all(key -> begin
        entry = get(cache.entries, key, nothing)
        entry !== nothing && entry.state == :hot && entry.expires_at > now
    end, value.resources) || return false
    for key in value.resources
        cache.entries[key].expires_at = now + cache.ttl_seconds
    end
    return true
end

function checked_profile_registry(profile::AbstractScientificProfile)
    RequiredInterfaces.check_interface_implemented(AbstractScientificProfile, typeof(profile)) === true ||
        throw(ArgumentError("scientific profile is missing required hooks"))
    registry = operation_registry(profile)
    registry isa OperationRegistry && !isempty(registry.operations) ||
        throw(ArgumentError("scientific profile requires a nonempty operation registry"))
    all(spec -> spec.execution_mode == :supervised, values(registry.operations)) ||
        throw(ArgumentError("profile operations must run in the supervised executor"))
    return registry
end

function prepare_profile!(profile::AbstractScientificProfile, context::ExecutionContext,
        parameters::Dict{String,Any})
    validated = validate_preparation(profile, parameters)
    validated isa Dict{String,Any} || throw(ArgumentError("preparation validation must return normalized inputs"))
    key = input_hash("runtime.preparation", validated)
    cache, cache_key = context.prepared_cache, "preparation:" * key
    prior = lock(cache.lock) do
        prune_prepared!(cache)
        entry = get(cache.entries, cache_key, nothing)
        if entry !== nothing && entry.state == :hot &&
                !(entry.value isa PreparedWorkload && preparation_dependencies!(cache, entry.value))
            delete!(cache.entries, cache_key)
            entry = nothing
        end
        entry === nothing ? "miss" : entry.state == :hot ? "hit" : entry.state == :warming ? "shared" : "failed"
    end
    workload = prepare_resource!(cache, cache_key) do
        check_canceled(context)
        value = prepare!(profile, context, validated)
        check_canceled(context)
        value isa PreparedWorkload || throw(ArgumentError("preparation must return PreparedWorkload"))
        value
    end
    lock(cache.lock) do
        if !preparation_dependencies!(cache, workload)
            delete!(cache.entries, cache_key)
            throw(ArgumentError("preparation model was not retained within its cache budget"))
        end
    end
    return Dict{String,Any}("preparation_input_hash" => key, "evidence" => workload.evidence,
        "cache_status" => prior, "idle_ttl_seconds" => cache.ttl_seconds)
end

"""
    preparation_snapshot(cache, key)

Inspect one normalized preparation identity without extending its idle lifetime.
Report ready only while its evidence and every dependent model remain retained.
The remaining lifetime is the earliest expiry among those required entries.
"""
function preparation_snapshot(cache::PreparedResourceCache, key::String)
    occursin(r"^[a-f0-9]{64}$", key) || throw(ArgumentError("invalid preparation identity"))
    return lock(cache.lock) do
        prune_prepared!(cache)
        entry = get(cache.entries, "preparation:" * key, nothing)
        ready = entry !== nothing && entry.state == :hot && entry.value isa PreparedWorkload
        remaining = 0.0
        if ready
            now = cache.clock()
            dependencies = [get(cache.entries, name, nothing) for name in entry.value.resources]
            ready = entry.expires_at > now && all(model -> model !== nothing &&
                model.state == :hot && model.expires_at > now, dependencies)
            if ready
                remaining = minimum([entry.expires_at; [model.expires_at for model in dependencies]]) - now
            else
                delete!(cache.entries, "preparation:" * key)
            end
        end
        return Dict{String,Any}("ready"=>ready, "preparation_input_hash"=>key,
            "remaining_seconds"=>remaining)
    end
end

"""
    profile_executor_main(profile; cache=PreparedResourceCache())

Serve the existing local JSON framing protocol for one approved profile. A
`prepare` command explicitly validates and runs its representative workload.
Bootstrap only announces the command reader; the host must correlate preparation
evidence with the live process, environment fingerprint and assignment generation.
Each invocation owns its cache; no prepared model is shared between processes.
"""
function profile_executor_main(profile::AbstractScientificProfile; cache=PreparedResourceCache())
    registry = checked_profile_registry(profile)
    return executor_main(registry, _ -> load_profile!(profile); prepared_cache=cache,
        prepare=(context, parameters) -> begin
            progress!(context, 0.0, "loading_environment"; message="Loading approved numerical environment")
            load_profile!(profile)
            check_canceled(context)
            # Importing an approved package may add methods after this reader
            # task began. Enter its preparation hook in the current world.
            Base.invokelatest(prepare_profile!, profile, context, parameters)
        end,
        preparation_status=key -> preparation_snapshot(cache, key),
        cleanup=() -> cleanup!(profile, cache))
end
