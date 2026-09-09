"""
    ControlConfig(endpoint, profiles, workers; limits=AssignmentLimits(), artifacts=nothing)

Hold passive, operator-provisioned worker identities and environments. This does
not enroll or approve workers, connect to NATS, create streams, or load engines.
Changing an enrolled identity's trust binding requires explicit reconciliation.
"""
struct ControlConfig
    "Coordinator-only broker credentials."
    endpoint::BrokerEndpoint
    "Approved passive profile definitions."
    profiles::ProfileRegistry
    "Provisioned identities eligible for administrator enrollment."
    workers::Dict{String,WorkerTrust}
    "Transactional assignment admission bounds."
    limits::AssignmentLimits
    "Optional private artifact reader, never the publisher's public digest store."
    artifacts::Union{Nothing,AbstractRuntimeArtifacts}
    function ControlConfig(endpoint::BrokerEndpoint, profiles::ProfileRegistry, workers;
            limits::AssignmentLimits=AssignmentLimits(),artifacts::Union{Nothing,AbstractRuntimeArtifacts}=nothing)
        1 <= length(profiles.definitions) <= 64 || throw(ArgumentError("control requires 1:64 profiles"))
        trusts = collect(workers)
        1 <= length(trusts) <= 128 && all(w -> w isa WorkerTrust, trusts) ||
            throw(ArgumentError("control requires 1:128 provisioned workers"))
        allunique(w.worker_id for w in trusts) && allunique(w.credential_ref for w in trusts) ||
            throw(ArgumentError("duplicate provisioned worker identity or credential reference"))
        for profile in values(profiles.definitions)
            Protocol.runtime_token(profile.id)
        end
        all(w -> all(id -> haskey(profiles.definitions, id), w.profiles), trusts) ||
            throw(ArgumentError("worker references an unregistered profile"))
        return new(endpoint, profiles, Dict(w.worker_id => w for w in trusts), limits, artifacts)
    end
end

Base.show(io::IO, config::ControlConfig) = print(io,
    "ControlConfig(", length(config.workers), " provisioned workers, server-owned credentials)")

function configuration_records(data, name; maximum)
    entries = get(data, name, nothing)
    entries isa AbstractVector && 1 <= length(entries) <= maximum &&
        all(item -> item isa AbstractDict, entries) ||
        throw(ArgumentError("$name must be a nonempty bounded array of TOML tables"))
    return entries
end

function configured_broker(data, base)
    table = strict_keys(config_table(data, "broker"),
        ("url", "password_file", "ca_file", "certificate_file", "key_file", "server_name",
         "allow_loopback_plaintext"), "broker")
    kwargs = Dict{Symbol,Any}()
    for name in ("ca_file", "certificate_file", "key_file")
        haskey(table, name) && (kwargs[Symbol(name)] = config_path(base, table[name]))
    end
    for name in ("server_name", "allow_loopback_plaintext")
        haskey(table, name) && (kwargs[Symbol(name)] = table[name])
    end
    get(table, "url", nothing) isa AbstractString || throw(ArgumentError("broker url is required"))
    endpoint = BrokerEndpoint(table["url"], config_path(base, get(table, "password_file", "")); kwargs...)
    broker_file(endpoint.password_file; private=true, max_bytes=4096)
    for (path, private) in ((endpoint.ca_file, false), (endpoint.certificate_file, false), (endpoint.key_file, true))
        path === nothing || broker_file(path; private)
    end
    return endpoint
end

function configured_profiles(data, base)
    profiles = ProfileRegistry()
    for entry in configuration_records(data, "profiles"; maximum=64)
        strict_keys(entry, ("id", "version", "kind", "isolation", "environment", "fingerprint",
            "operations", "preparation", "budget", "protocol_version"), "profiles")
        isolation = get(entry, "isolation", "trusted_process")
        isolation in ("trusted_process", "container") || throw(ArgumentError("invalid profile isolation"))
        kind = get(entry, "kind", "scientific")
        kind in ("scientific", "terminal") || throw(ArgumentError("invalid profile kind"))
        version = get(entry, "version", "1.0.0")
        version isa AbstractString || throw(ArgumentError("profile version must be a string"))
        environment = get(entry, "environment", "")
        isolation == "trusted_process" && (environment = config_path(base, environment))
        budget = strict_keys(config_table(entry, "budget"),
            ("cpus", "memory_bytes", "pids", "scratch_bytes", "prepare_seconds", "job_seconds"), "profile budget")
        operations = get(entry, "operations", String[])
        operations isa AbstractVector && all(v -> v isa AbstractString, operations) ||
            throw(ArgumentError("profile operations must be an array of identifiers"))
        register!(profiles, ProfileDefinition(get(entry, "id", ""), environment, get(entry, "fingerprint", "");
            version=VersionNumber(version), kind=Symbol(kind), isolation=Symbol(isolation), operations,
            preparation=get(entry, "preparation", "default"),
            budget=ResourceBudget(; (Symbol(k) => v for (k,v) in budget)...),
            protocol_version=get(entry, "protocol_version", 2)))
    end
    return profiles
end

"""
    read_control_config(path) -> ControlConfig

Read strict schema-1 TOML and validate private credential files without reading
their contents or connecting. Relative paths resolve beside this file. Native
environment existence and actual executor preparation are agent responsibilities;
the coordinator never imports or opens a scientific environment.
"""
function read_control_config(path::AbstractString)
    data = TOML.parsefile(path)
    strict_keys(data, ("schema_version", "broker", "profiles", "workers", "assignments", "artifacts"), "control root")
    get(data, "schema_version", nothing) === 1 || throw(ArgumentError("unsupported control configuration version"))
    base = dirname(abspath(path))
    endpoint = configured_broker(data, base)
    profiles = configured_profiles(data, base)
    workers = WorkerTrust[]
    for entry in configuration_records(data, "workers"; maximum=128)
        strict_keys(entry, ("id", "credential_ref", "profiles", "capacity"), "workers")
        ids = get(entry, "profiles", nothing)
        ids isa AbstractVector && all(v -> v isa AbstractString, ids) ||
            throw(ArgumentError("worker profiles must be an array of identifiers"))
        push!(workers, WorkerTrust(get(entry, "id", ""), get(entry, "credential_ref", ""), ids;
            capacity=get(entry, "capacity", 1)))
    end
    limits = strict_keys(config_table(data, "assignments"), ("total", "per_owner", "per_run"), "assignments")
    return ControlConfig(endpoint, profiles, workers;
        limits=AssignmentLimits(; (Symbol(k) => v for (k,v) in limits)...),artifacts=configured_runtime_artifacts(data,base))
end
