"""
    AgentConfig(worker_id, endpoint, profiles, scratch_root;
        capacity=1, container_runtime=:auto, artifacts=nothing)

Hold a provisioned agent identity and approved local environments. Runtime
selection is operator-owned (:auto, :podman, or :docker), not browser input.
This configuration does not claim that any profile passes its local preflight.
"""
struct AgentConfig
    "Broker-bound worker subject identity."
    worker_id::String
    "Worker-only broker credentials."
    endpoint::BrokerEndpoint
    "Approved profiles; the resource supervisor may advertise a verified subset."
    profiles::ProfileRegistry
    "Dedicated directory for owned resource receipts and scratch."
    scratch_root::String
    "Maximum occupied assignments, including cleanup."
    capacity::Int
    "Approved container command selection, never a socket passed to user code."
    container_runtime::Symbol
    "Optional job-scoped artifact writer; credentials never enter executors."
    artifacts::Union{Nothing,AbstractRuntimeArtifacts}
    function AgentConfig(worker_id::AbstractString, endpoint::BrokerEndpoint,
            profiles::ProfileRegistry, scratch_root::AbstractString;
            capacity=1, container_runtime::Symbol=:auto,artifacts::Union{Nothing,AbstractRuntimeArtifacts}=nothing)
        id = Protocol.runtime_token(worker_id)
        1 <= length(profiles.definitions) <= 64 || throw(ArgumentError("agent requires 1:64 approved profiles"))
        foreach(p -> Protocol.runtime_token(p.id), values(profiles.definitions))
        capacity isa Integer && !(capacity isa Bool) && 1 <= capacity <= 256 ||
            throw(ArgumentError("agent capacity must be in 1:256"))
        container_runtime in (:auto, :podman, :docker) || throw(ArgumentError("unsupported container runtime"))
        root = abspath(scratch_root)
        root in ("/", homedir(), pwd(), dirname(pwd()), tempdir()) &&
            throw(ArgumentError("agent scratch_root must be a dedicated directory"))
        return new(id, endpoint, profiles, root, capacity, container_runtime, artifacts)
    end
end

Base.show(io::IO, config::AgentConfig) = print(io, "AgentConfig(", config.worker_id, ", server-owned resources)")

"""
    read_agent_config(path) -> AgentConfig

Read strict schema-1 agent TOML without connecting, importing engines or starting
resources. Broker/profile tables have the same grammar as control configuration.
The agent table owns worker_id, scratch_root, capacity and container_runtime.
"""
function read_agent_config(path::AbstractString)
    data = TOML.parsefile(path)
    strict_keys(data, ("schema_version", "agent", "broker", "profiles", "artifacts"), "agent root")
    get(data, "schema_version", nothing) === 1 || throw(ArgumentError("unsupported agent configuration version"))
    base = dirname(abspath(path))
    agent = strict_keys(config_table(data, "agent"),
        ("worker_id", "scratch_root", "capacity", "container_runtime"), "agent")
    runtime = get(agent, "container_runtime", "auto")
    runtime in ("auto", "podman", "docker") || throw(ArgumentError("unsupported container runtime"))
    root = config_path(base, get(agent, "scratch_root", "state/agent"))
    root in (base, dirname(base)) && throw(ArgumentError("agent scratch_root must be a dedicated directory"))
    return AgentConfig(get(agent, "worker_id", ""), configured_broker(data, base),
        configured_profiles(data, base), root;
        capacity=get(agent, "capacity", 1), container_runtime=Symbol(runtime),artifacts=configured_runtime_artifacts(data,base))
end
