"""
    ResourceBudget(; cpus=1, memory_bytes=1073741824, pids=128,
                   scratch_bytes=268435456, prepare_seconds=300, job_seconds=120)

Declare finite executor limits. These are requirements, not evidence that a
host enforces them; the host agent must pass its isolation preflight.
"""
struct ResourceBudget
    "Required CPU quota, in logical CPU units."
    cpus::Float64
    "Maximum memory in bytes."
    memory_bytes::Int
    "Maximum process count."
    pids::Int
    "Maximum disposable writable scratch space in bytes."
    scratch_bytes::Int
    "Preparation deadline in seconds."
    prepare_seconds::Float64
    "Scientific job deadline in seconds."
    job_seconds::Float64
    function ResourceBudget(; cpus=1, memory_bytes=1024^3, pids=128,
            scratch_bytes=256 * 1024^2, prepare_seconds=300, job_seconds=120)
        cpus isa Real && !(cpus isa Bool) && isfinite(cpus) && 0 < cpus <= 256 ||
            throw(ArgumentError("CPU quota must be finite and in (0, 256]"))
        all(v -> v isa Integer && !(v isa Bool) && 0 < v <= typemax(Int),
            (memory_bytes, pids, scratch_bytes)) ||
            throw(ArgumentError("memory, PID and scratch limits must be positive integers"))
        pids <= 65536 || throw(ArgumentError("PID limit exceeds 65536"))
        all(v -> v isa Real && !(v isa Bool) && isfinite(v) && 0 < v <= 3600,
            (prepare_seconds, job_seconds)) ||
            throw(ArgumentError("executor deadlines must be in (0, 3600] seconds"))
        return new(cpus, memory_bytes, pids, scratch_bytes, prepare_seconds, job_seconds)
    end
end

"""
    ProfileDefinition(id, environment, fingerprint; version=v"1.0.0",
        kind=:scientific, isolation=:trusted_process, operations=(),
        preparation="default", budget=ResourceBudget(), protocol_version=2)

Register an approved environment without importing it. Environment is an
operator-owned reference (native project or digest-pinned image); fingerprint
is its lowercase SHA-256 identity. Native profiles are trusted execution, not
sandboxes. Terminal profiles require container isolation and no job operations.
"""
struct ProfileDefinition
    "Stable approved profile identifier."
    id::String
    "Profile contract version."
    version::VersionNumber
    "Scientific executor or private terminal."
    kind::Symbol
    "Trusted process or restricted container."
    isolation::Symbol
    "Server-owned environment reference, never a browser launch argument."
    environment::String
    "Immutable environment fingerprint."
    fingerprint::String
    "Allowlisted scientific wire operation identifiers."
    operations::Tuple{Vararg{String}}
    "Registered preparation recipe identity."
    preparation::String
    "Required resource and time limits."
    budget::ResourceBudget
    "Supported assigned-job/control wire version."
    protocol_version::Int
    function ProfileDefinition(id, environment, fingerprint;
            version::VersionNumber=v"1.0.0", kind::Symbol=:scientific,
            isolation::Symbol=:trusted_process, operations=(), preparation="default",
            budget::ResourceBudget=ResourceBudget(), protocol_version=2)
        kind in (:scientific, :terminal) || throw(ArgumentError("invalid profile kind"))
        isolation in (:trusted_process, :container) || throw(ArgumentError("unsupported isolation policy"))
        protocol_version === 2 || throw(ArgumentError("profile requires assigned protocol v2"))
        environment isa AbstractString && 0 < ncodeunits(environment) <= 1024 &&
            !occursin(r"[\x00\r\n]", environment) || throw(ArgumentError("invalid environment reference"))
        fingerprint isa AbstractString && occursin(r"^[a-f0-9]{64}$", fingerprint) ||
            throw(ArgumentError("environment fingerprint must be a lowercase SHA-256 digest"))
        if isolation == :container
            occursin(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*@sha256:[a-f0-9]{64}$", environment) ||
                throw(ArgumentError("container environment must be digest-pinned"))
            endswith(environment, "@sha256:" * fingerprint) ||
                throw(ArgumentError("image digest and environment fingerprint must agree"))
        end
        capabilities = Tuple(checked_id.(operations))
        allunique(capabilities) || throw(ArgumentError("duplicate profile operations"))
        if kind == :terminal
            isolation == :container && isempty(capabilities) ||
                throw(ArgumentError("terminals require container isolation and a separate byte channel"))
        else
            !isempty(capabilities) || throw(ArgumentError("scientific profile must declare operations"))
        end
        return new(checked_id(id), version, kind, isolation, environment, fingerprint,
            capabilities, checked_id(preparation), budget, protocol_version)
    end
end

"""
    ProfileRegistry()

Store passive, operator-approved profiles. Registration does not open their
environment, pull images, connect to a broker, prepare or run scientific work.
"""
struct ProfileRegistry
    "Approved immutable descriptions keyed by profile ID."
    definitions::Dict{String,ProfileDefinition}
end
ProfileRegistry() = ProfileRegistry(Dict{String,ProfileDefinition}())

"""
    register!(registry::ProfileRegistry, definition::ProfileDefinition)

Add an approved profile. Duplicate IDs are errors even if versions differ;
upgrades require explicit deployment/reconciliation rather than silent replacement.
"""
function register!(registry::ProfileRegistry, definition::ProfileDefinition)
    haskey(registry.definitions, definition.id) && throw(ArgumentError("profile ID is already registered"))
    registry.definitions[definition.id] = definition
    return definition
end

"""
    validate_requirements(applications, profiles)

Reject catalogue roles that reference missing profiles before assigning scientific work.
This check is passive and never loads an application or executor environment.
"""
function validate_requirements(applications::ApplicationRegistry, profiles::ProfileRegistry)
    for app in values(applications.definitions), role in app.requirements, profile in role.profiles
        haskey(profiles.definitions, profile) ||
            throw(ArgumentError("application $(app.id) requires unregistered profile $profile"))
    end
    return nothing
end
