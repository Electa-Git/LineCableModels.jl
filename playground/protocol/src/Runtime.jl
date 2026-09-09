"""
    RUNTIME_PROTOCOL_VERSION

Version of assigned runtime-control messages. The v1 scientific payload remains
embedded unchanged; v2 transport and fencing are separate from legacy subjects.
"""
const RUNTIME_PROTOCOL_VERSION = "2.0"

"""Passive records used by the assigned runtime protocol; never executable code."""
abstract type RuntimeRecord end

"""
    WorkerProbe

Request a fresh worker report from a coordinator-only broker subject. The agent
must reconcile any prior coordinator's leases before answering a new incarnation.
"""
struct WorkerProbe <: RuntimeRecord
    "Assigned control protocol version."
    protocol_version::String
    "Intended approved worker identity."
    worker_id::String
    "Current coordinator incarnation UUID."
    coordinator_id::String
    "Fresh, single-use report challenge UUID."
    challenge::String
end

"""
    ProfileAdvertisement

Describe an installed, approved environment without claiming that an executor
has prepared it. Registration and matching remain coordinator responsibilities.
"""
struct ProfileAdvertisement <: RuntimeRecord
    "Stable profile identifier."
    profile_id::String
    "Profile contract version."
    version::String
    "SHA-256 environment identity."
    fingerprint::String
end

"""
    WorkerAnnouncement

Report a worker incarnation on its broker-authorized identity subject. A payload
worker ID is not authentication. Coordinator challenge and increasing sequence
must also match before this report establishes current presence.
"""
struct WorkerAnnouncement <: RuntimeRecord
    "Assigned control protocol version."
    protocol_version::String
    "Broker-authorized worker identity."
    worker_id::String
    "Fresh UUID for this host-agent process."
    boot_id::String
    "Current coordinator incarnation UUID."
    coordinator_id::String
    "Coordinator-issued freshness challenge UUID."
    challenge::String
    "Strictly increasing report sequence within this boot."
    sequence::Int
    "Total simultaneously reserved executor slots."
    capacity::Int
    "Installed profiles, not prepared process states."
    profiles::Vector{ProfileAdvertisement}
end

"""
    AssignmentFence

Bind all assigned operations, events and results to the same owner, run, role,
worker incarnation and executor generation. Possessing this passive record does
not grant authority: the coordinator and agent must validate the live lease.
"""
struct AssignmentFence <: RuntimeRecord
    "Lease UUID."
    lease_id::String
    "Application-run UUID."
    run_id::String
    "Authenticated principal; never accepted from browser ownership claims."
    owner::String
    "Application-local runtime role."
    role::String
    "Approved worker identity."
    worker_id::String
    "Exact host-agent incarnation UUID."
    worker_boot::String
    "Exact coordinator incarnation UUID."
    coordinator_id::String
    "Approved profile identifier."
    profile_id::String
    "Approved profile contract version."
    profile_version::String
    "Approved environment SHA-256 identity."
    fingerprint::String
    "Monotonically increasing executor generation for this run and role."
    generation::Int
end

"""
    LeaseControl

Request grant, renewal or release through a coordinator-only subject. Repeated
request IDs are idempotent; revisions cannot move backwards. Duration starts
from the receiving agent's monotonic clock, never a browser wall clock.
"""
struct LeaseControl <: RuntimeRecord
    "Assigned control protocol version."
    protocol_version::String
    "Idempotent control request UUID."
    request_id::String
    "One of grant, renew or release."
    action::String
    "Complete assignment identity."
    fence::AssignmentFence
    "Increasing grant/renew/release revision."
    revision::Int
    "Bounded lease duration in milliseconds; zero for release."
    duration_ms::Int
end

"""
    LeaseAcknowledgement

Acknowledge the exact lease control revision from the intended worker. A grant
is unusable until its accepted acknowledgement is validated by the coordinator.
"""
struct LeaseAcknowledgement <: RuntimeRecord
    "Assigned control protocol version."
    protocol_version::String
    "Original idempotent control request UUID."
    request_id::String
    "Complete assignment identity."
    fence::AssignmentFence
    "Acknowledged control revision."
    revision::Int
    "Whether the agent applied the requested control."
    accepted::Bool
    "Bounded public reason code; not a raw exception or private path."
    reason::String
end

"""
    PreparedExecution

Identify the exact prepared child selected for a scientific request. This passive
record is not authority: the agent must match all three fields against its live
assignment-owned process and retained preparation before admitting work.
"""
struct PreparedExecution <: RuntimeRecord
    "Exact scientific resource UUID, replaced when its process is replaced."
    executor_id::String
    "Supervisor generation which produced the preparation."
    executor_generation::Int
    "Normalized preparation-input SHA-256 digest."
    preparation_key::String
end

"""
    AssignedJob

Carry the existing validated scientific request under an explicit v2 assignment.
Its session correlation must equal the application-run UUID. Delivery uses only
the intended worker's subject; a shared v1 queue cannot enforce this contract.
"""
struct AssignedJob <: RuntimeRecord
    "Assigned control protocol version."
    protocol_version::String
    "Complete assignment identity."
    fence::AssignmentFence
    "Existing passive scientific operation and normalized inputs."
    request::JobRequest
    "Exact prepared child and model; delayed jobs cannot use a replacement."
    execution::PreparedExecution
end

"""
    AssignedResult

Carry a durable scientific outcome with its exact assignment fence. The inner
result retains existing scientific provenance; its worker and environment must
match the assigned profile. A consumer must additionally match the submitted job
identity and input hash and reject results from revoked or replaced assignments.
"""
struct AssignedResult <: RuntimeRecord
    "Assigned runtime protocol version."
    protocol_version::String
    "Exact run, owner, worker incarnation and generation."
    fence::AssignmentFence
    "Existing validated scientific result, artifact or failure."
    result::JobResult
    "Prepared child/model targeted by the original request."
    execution::PreparedExecution
end

"""
    ScientificCommand

Request explicit preparation, retained-state inspection or cancellation on one
live assignment. The coordinator assigns increasing revisions; retries preserve
the complete record. These transient commands do not carry executable code or
replace durable scientific jobs. Status never starts or prepares a process.
"""
struct ScientificCommand <: RuntimeRecord
    "Assigned runtime protocol version."
    protocol_version::String
    "Fresh command correlation UUID, retained for retries."
    request_id::String
    "Complete assignment identity; the agent must still verify its live lease."
    fence::AssignmentFence
    "Increasing command revision within this assignment."
    revision::Int
    "One of prepare, status, cancel (current request), or cancel_job (including future delivery)."
    action::String
    "Passive preparation inputs; empty for inspection and cancellation."
    parameters::Dict{String,Any}
    "Exact request/job to cancel; nothing for other commands."
    target_id::Union{Nothing,String}
end

"""
    ScientificReport

Return bounded scientific activity for an exact command and assignment. Ready
requires fresh child inspection and a finite validity interval; receivers must
subtract total request latency and recheck local lease authority. This is neither
a durable result nor a reusable grant. Raw output, private paths and exceptions
are excluded.
"""
struct ScientificReport <: RuntimeRecord
    "Assigned runtime protocol version."
    protocol_version::String
    "Exact command correlation UUID."
    request_id::String
    "Complete assignment identity."
    fence::AssignmentFence
    "Exact command revision."
    revision::Int
    "Whether the requested action was accepted, not whether preparation is ready."
    accepted::Bool
    "Fixed public reason token."
    reason::String
    "Idle, starting, preparing, executing, failed or closing."
    phase::String
    "Unknown, cold, preparing, ready or failed."
    preparation::String
    "Current executor identity, if one has been acquired."
    executor_id::Union{Nothing,String}
    "Current process generation; zero before first start."
    executor_generation::Int
    "Current scientific request correlation, if any."
    current_request_id::Union{Nothing,String}
    "Verified preparation input identity, only when ready."
    preparation_key::Union{Nothing,String}
    "Progress in thousandths, from zero to 1000."
    progress_milli::Int
    "Elapsed request time in milliseconds."
    elapsed_ms::Int
    "Bounded count of output lines; never their contents."
    output_lines::Int
    "Fixed failure token, if any."
    failure::Union{Nothing,String}
    "Maximum remaining readiness validity in milliseconds, capped at five seconds."
    valid_for_ms::Int
end

StructTypes.StructType(::Type{T}) where {T<:RuntimeRecord} = StructTypes.Struct()
==(left::T, right::T) where {T<:RuntimeRecord} =
    all(getfield(left, field) == getfield(right, field) for field in fieldnames(T))

function runtime_token(value::AbstractString)
    occursin(r"^[a-z][a-z0-9_-]{0,63}$", value) ||
        throw(ArgumentError("runtime token must contain 1–64 lowercase safe characters"))
    return String(value)
end

function runtime_uuid(value::AbstractString)
    parsed = tryparse(UUID, value)
    parsed !== nothing && string(parsed) == value ||
        throw(ArgumentError("runtime identity must be a lowercase UUID"))
    return value
end

function runtime_version(value::AbstractString)
    ncodeunits(value) <= 48 || throw(ArgumentError("profile version is too long"))
    parsed = tryparse(VersionNumber, value)
    parsed !== nothing && string(parsed) == value ||
        throw(ArgumentError("profile version must use explicit semantic version syntax"))
    return value
end

function runtime_fingerprint(value::AbstractString)
    occursin(r"^[a-f0-9]{64}$", value) ||
        throw(ArgumentError("environment identity must be a lowercase SHA-256 digest"))
    return value
end

function runtime_sequence(value::Integer)
    0 < value <= 9_007_199_254_740_991 ||
        throw(ArgumentError("runtime sequence must be a positive JSON-safe integer"))
    return value
end

function validate(profile::ProfileAdvertisement)
    runtime_token(profile.profile_id)
    runtime_version(profile.version)
    runtime_fingerprint(profile.fingerprint)
    return profile
end

function validate(probe::WorkerProbe)
    probe.protocol_version == RUNTIME_PROTOCOL_VERSION ||
        throw(ArgumentError("unsupported assigned runtime protocol"))
    runtime_token(probe.worker_id)
    foreach(runtime_uuid, (probe.coordinator_id, probe.challenge))
    return probe
end

function validate(report::WorkerAnnouncement)
    report.protocol_version == RUNTIME_PROTOCOL_VERSION ||
        throw(ArgumentError("unsupported assigned runtime protocol"))
    runtime_token(report.worker_id)
    foreach(runtime_uuid, (report.boot_id, report.coordinator_id, report.challenge))
    runtime_sequence(report.sequence)
    1 <= report.capacity <= 256 || throw(ArgumentError("worker capacity must be in 1–256"))
    length(report.profiles) <= 64 ||
        throw(ArgumentError("worker may advertise at most 64 installed profiles"))
    allunique(profile.profile_id for profile in report.profiles) ||
        throw(ArgumentError("duplicate profile advertisement"))
    foreach(validate, report.profiles)
    return report
end

function validate(fence::AssignmentFence)
    foreach(runtime_uuid, (fence.lease_id, fence.run_id, fence.worker_boot, fence.coordinator_id))
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.@-]{0,95}$", fence.owner) ||
        throw(ArgumentError("invalid assignment principal"))
    foreach(runtime_token, (fence.role, fence.worker_id, fence.profile_id))
    runtime_version(fence.profile_version)
    runtime_fingerprint(fence.fingerprint)
    runtime_sequence(fence.generation)
    return fence
end

function validate(control::LeaseControl)
    control.protocol_version == RUNTIME_PROTOCOL_VERSION ||
        throw(ArgumentError("unsupported assigned runtime protocol"))
    runtime_uuid(control.request_id)
    validate(control.fence)
    runtime_sequence(control.revision)
    control.action in ("grant", "renew", "release") ||
        throw(ArgumentError("unsupported lease action"))
    if control.action == "release"
        control.duration_ms == 0 || throw(ArgumentError("release cannot extend a lease"))
    else
        100 <= control.duration_ms <= 60_000 ||
            throw(ArgumentError("lease duration must be in 100–60000 milliseconds"))
    end
    return control
end

function validate(ack::LeaseAcknowledgement)
    ack.protocol_version == RUNTIME_PROTOCOL_VERSION ||
        throw(ArgumentError("unsupported assigned runtime protocol"))
    runtime_uuid(ack.request_id)
    validate(ack.fence)
    runtime_sequence(ack.revision)
    runtime_token(ack.reason)
    return ack
end

function validate(execution::PreparedExecution)
    runtime_uuid(execution.executor_id)
    runtime_sequence(execution.executor_generation)
    runtime_fingerprint(execution.preparation_key)
    return execution
end

function validate(job::AssignedJob)
    job.protocol_version == RUNTIME_PROTOCOL_VERSION ||
        throw(ArgumentError("unsupported assigned runtime protocol"))
    validate(job.fence)
    validate(job.request)
    validate(job.execution)
    runtime_uuid(job.request.job_id)
    job.request.session_id == job.fence.run_id ||
        throw(ArgumentError("scientific request must belong to its assigned run"))
    return job
end

function validate(outcome::AssignedResult)
    outcome.protocol_version == RUNTIME_PROTOCOL_VERSION ||
        throw(ArgumentError("unsupported assigned runtime protocol"))
    validate(outcome.fence)
    validate(outcome.result)
    validate(outcome.execution)
    runtime_uuid(outcome.result.job_id)
    outcome.result.worker_id == outcome.fence.worker_id &&
        outcome.result.environment_fingerprint == outcome.fence.fingerprint ||
        throw(ArgumentError("scientific result does not match its assigned worker environment"))
    return outcome
end

function validate(command::ScientificCommand)
    command.protocol_version == RUNTIME_PROTOCOL_VERSION || throw(ArgumentError("unsupported assigned runtime protocol"))
    runtime_uuid(command.request_id); validate(command.fence); runtime_sequence(command.revision)
    command.action in ("prepare","status","cancel","cancel_job") || throw(ArgumentError("unsupported scientific control action"))
    normalize_wire(command.parameters)
    ncodeunits(JSON3.write(command.parameters)) <= 65536 || throw(ArgumentError("preparation inputs exceed 64 KiB"))
    if command.action in ("cancel","cancel_job")
        command.target_id !== nothing && isempty(command.parameters) || throw(ArgumentError("cancel requires only a target request"))
        runtime_uuid(command.target_id)
    else
        command.target_id === nothing || throw(ArgumentError("unexpected cancellation target"))
        command.action == "prepare" || isempty(command.parameters) || throw(ArgumentError("status cannot carry preparation inputs"))
    end
    return command
end

function validate(report::ScientificReport)
    report.protocol_version == RUNTIME_PROTOCOL_VERSION || throw(ArgumentError("unsupported assigned runtime protocol"))
    runtime_uuid(report.request_id); validate(report.fence); runtime_sequence(report.revision); runtime_token(report.reason)
    report.phase in ("idle","starting","preparing","executing","failed","closing") || throw(ArgumentError("invalid scientific phase"))
    report.preparation in ("unknown","cold","preparing","ready","failed") || throw(ArgumentError("invalid preparation state"))
    report.executor_id === nothing || runtime_uuid(report.executor_id)
    report.current_request_id === nothing || runtime_uuid(report.current_request_id)
    report.preparation_key === nothing || runtime_fingerprint(report.preparation_key)
    report.failure === nothing || runtime_token(report.failure)
    0 <= report.executor_generation <= 9_007_199_254_740_991 &&
        0 <= report.progress_milli <= 1000 && 0 <= report.elapsed_ms <= 9_007_199_254_740_991 &&
        0 <= report.output_lines <= 1_000_000 && 0 <= report.valid_for_ms <= 5000 ||
        throw(ArgumentError("scientific report counters exceed supported bounds"))
    if report.preparation == "ready"
        report.accepted && report.phase == "idle" && report.executor_id !== nothing && report.executor_generation > 0 &&
            report.preparation_key !== nothing && report.valid_for_ms > 0 || throw(ArgumentError("ready report lacks finite evidence"))
    else
        report.preparation_key === nothing && report.valid_for_ms == 0 || throw(ArgumentError("non-ready report cannot retain readiness"))
    end
    return report
end

# Validate the JSON shape before typed decoding, which may otherwise coerce a
# boolean to an integer or silently ignore unrecognized authority-related keys.
function runtime_shape(::Type{T}, value) where {T}
    if T isa Union
        any(Base.uniontypes(T)) do variant
            try
                runtime_shape(variant, value)
                true
            catch error
                error isa ArgumentError || rethrow()
                false
            end
        end || throw(ArgumentError("runtime field has an incompatible type"))
    elseif T <: RuntimeRecord || T in (JobRequest, JobResult, ArtifactReference, FailureInfo)
        value isa JSON3.Object || throw(ArgumentError("runtime record must be an object"))
        expected = Set(fieldnames(T))
        length(value) == length(expected) && Set(keys(value)) == expected ||
            throw(ArgumentError("runtime record has missing, duplicate or unknown fields"))
        for name in fieldnames(T)
            runtime_shape(fieldtype(T, name), value[name])
        end
    elseif T <: AbstractVector
        value isa JSON3.Array || throw(ArgumentError("runtime collection must be an array"))
        foreach(item -> runtime_shape(eltype(T), item), value)
    elseif T <: AbstractDict
        value isa JSON3.Object || throw(ArgumentError("runtime mapping must be an object"))
        normalize_wire(value)
    elseif T === Bool
        value isa Bool || throw(ArgumentError("runtime boolean must be JSON true or false"))
    elseif T <: Integer
        value isa Integer && !(value isa Bool) && typemin(T) <= value <= typemax(T) ||
            throw(ArgumentError("runtime integer must be an exact JSON integer"))
    elseif T <: AbstractString
        value isa AbstractString || throw(ArgumentError("runtime text must be a string"))
    elseif T === Nothing
        value === nothing || throw(ArgumentError("runtime null field must be null"))
    else
        throw(ArgumentError("unsupported runtime wire field type"))
    end
    return nothing
end

"""
    decode_runtime_message(T, payload)

Decode a known passive runtime record with strict nested fields and bounded
payload size. Unknown types, extra fields, coercible booleans, malformed UUIDs
and unsupported versions are rejected. This function grants no authorization.
"""
function decode_runtime_message(::Type{T}, payload) where {T<:RuntimeRecord}
    bytes = payload isa AbstractString ? ncodeunits(payload) : length(payload)
    bytes <= MAX_WIRE_PAYLOAD_BYTES || throw(ArgumentError("runtime message exceeds 256 KiB"))
    object = try
        JSON3.read(payload)
    catch
        throw(ArgumentError("invalid runtime JSON"))
    end
    runtime_shape(T, object)
    return validate(JSON3.read(payload, T))
end

"""
    assigned_job_subject(fence)

Return the exact targeted v2 job subject. The broker permissions and consumer
filter must enforce this worker identity; this helper alone grants no access.
"""
function assigned_job_subject(fence::AssignmentFence)
    validate(fence)
    return "lcm.jobs.v2.$(fence.worker_id).$(fence.worker_boot).$(fence.lease_id).$(fence.generation)"
end

"""
    assigned_result_subject(fence, job_id)

Return one worker-scoped durable result subject with exact lease/generation and
job identity. Subject construction is not caller authorization.
"""
function assigned_result_subject(fence::AssignmentFence, job_id::AbstractString)
    validate(fence)
    runtime_uuid(job_id)
    return "lcm.results.v2.$(fence.worker_id).$(fence.worker_boot).$(fence.lease_id).$(fence.generation).$job_id"
end
