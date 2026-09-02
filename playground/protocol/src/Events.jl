const EVENT_TYPES = (
    "JobQueued",
    "JobClaimed",
    "JobStarted",
    "JobProgress",
    "JobLog",
    "JobSucceeded",
    "JobFailed",
    "JobCanceled",
    "JobExpired",
)
const TERMINAL_EVENTS = ("JobSucceeded", "JobFailed", "JobCanceled", "JobExpired")
const MEDIA_TYPE_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$"
const MAX_ARTIFACT_BYTES = Int64(1) << 40

"""
    JobEvent

Represent one ordered lifecycle event emitted by a worker.
"""
struct JobEvent
    "Implemented wire-contract version."
    protocol_version::String

    "Lifecycle event type."
    event_type::String

    "Job identifier."
    job_id::String

    "Monotonic sequence within this job."
    event_sequence::Int

    "UTC event timestamp in ISO-8601 form."
    timestamp::String

    "Worker identity, when a worker has claimed the job."
    worker_id::Union{Nothing,String}

    "Completed fraction in the closed interval [0, 1]."
    progress::Union{Nothing,Float64}

    "Concise execution stage."
    stage::Union{Nothing,String}

    "Human-readable status or log message."
    message::Union{Nothing,String}

    "Additional JSON-compatible event data."
    details::Dict{String,Any}
end

StructTypes.StructType(::Type{JobEvent}) = StructTypes.Struct()

function validate(event::JobEvent)
    event.protocol_version == PROTOCOL_VERSION || throw(ArgumentError(
        "unsupported protocol version: $(event.protocol_version)"
    ))
    event.event_type in EVENT_TYPES || throw(ArgumentError(
        "unsupported event type: $(event.event_type)"
    ))
    tryparse(UUID, event.job_id) === nothing && throw(ArgumentError(
        "event job identifier must be a UUID"
    ))
    event.event_sequence > 0 || throw(ArgumentError(
        "event_sequence must be positive"
    ))
    parse_utc_timestamp(event.timestamp)
    if !isnothing(event.worker_id)
        occursin(IDENTIFIER_PATTERN, event.worker_id) || throw(ArgumentError(
            "invalid event worker identifier"
        ))
    end
    if !isnothing(event.progress)
        0.0 <= event.progress <= 1.0 || throw(ArgumentError(
            "event progress must be between zero and one"
        ))
    end
    if !isnothing(event.stage)
        occursin(IDENTIFIER_PATTERN, event.stage) || throw(ArgumentError(
            "invalid event stage"
        ))
    end
    if !isnothing(event.message)
        validate_wire_text(event.message, "event message"; maximum=8 * 1024)
    end
    validate_wire_data(event.details)
    return event
end

isterminal(event::JobEvent) = event.event_type in TERMINAL_EVENTS

"""
    ArtifactReference

Describe a content-addressed result stored outside the broker message.
"""
struct ArtifactReference
    "Content-addressed artifact identifier."
    artifact_id::String

    "IANA media type."
    media_type::String

    "Artifact size [byte]."
    size::Int

    "SHA-256 digest of the stored bytes."
    sha256::String

    "Storage implementation identifier."
    storage_backend::String

    "Server-side retrieval reference."
    retrieval_reference::String
end

StructTypes.StructType(::Type{ArtifactReference}) = StructTypes.Struct()

function validate(artifact::ArtifactReference)
    occursin(r"^sha256:[0-9a-f]{64}$", artifact.artifact_id) ||
        throw(ArgumentError("invalid content-addressed artifact identifier"))
    0 <= artifact.size <= MAX_ARTIFACT_BYTES || throw(ArgumentError(
        "artifact size must be between zero and one TiB"
    ))
    occursin(r"^[0-9a-f]{64}$", artifact.sha256) || throw(ArgumentError(
        "artifact sha256 must be a lowercase digest"
    ))
    artifact.artifact_id == "sha256:$(artifact.sha256)" || throw(ArgumentError(
        "artifact identifier and digest do not match"
    ))
    occursin(MEDIA_TYPE_PATTERN, artifact.media_type) || throw(ArgumentError(
        "invalid artifact media type"
    ))
    occursin(IDENTIFIER_PATTERN, artifact.storage_backend) || throw(ArgumentError(
        "invalid artifact storage backend"
    ))
    artifact.retrieval_reference == "/artifacts/sha256/$(artifact.sha256)" ||
        throw(ArgumentError(
            "artifact retrieval reference must be its same-origin digest route"
        ))
    return artifact
end

"""
    FailureInfo

Represent a browser-safe terminal failure. Detailed stack traces remain in
worker logs and are correlated through `diagnostic_id`.
"""
struct FailureInfo
    "Stable failure category."
    category::String

    "Concise browser-safe message."
    message::String

    "Identifier correlating the failure with worker logs."
    diagnostic_id::String

    "Whether redelivery may succeed without changing inputs."
    retryable::Bool
end

StructTypes.StructType(::Type{FailureInfo}) = StructTypes.Struct()

function validate(failure::FailureInfo)
    occursin(IDENTIFIER_PATTERN, failure.category) || throw(ArgumentError(
        "invalid failure category"
    ))
    validate_wire_text(
        failure.message,
        "failure message";
        minimum=1,
        maximum=4 * 1024
    )
    tryparse(UUID, failure.diagnostic_id) === nothing && throw(ArgumentError(
        "failure diagnostic identifier must be a UUID"
    ))
    return failure
end

"""
    JobResult

Represent a durable terminal result and its complete computational provenance.
Exactly one of `inline_result`, `artifact`, or `failure` is populated.
"""
struct JobResult
    "Implemented wire-contract version."
    protocol_version::String

    "Job identifier."
    job_id::String

    "Registered operation name."
    operation::String

    "Operation result-schema version."
    schema_version::String

    "SHA-256 digest of normalized inputs."
    input_hash::String

    "LineCableModels version or source commit."
    engine_version::String

    "Digest identifying the locked worker environment."
    environment_fingerprint::String

    "Worker identity."
    worker_id::String

    "Cache disposition: `hit`, `miss`, or `bypass`."
    cache_status::String

    "UTC computation start timestamp in ISO-8601 form."
    started_at::String

    "UTC completion timestamp in ISO-8601 form."
    completed_at::String

    "Small JSON-compatible result, when returned inline."
    inline_result::Union{Nothing,Dict{String,Any}}

    "Reference to a larger result artifact."
    artifact::Union{Nothing,ArtifactReference}

    "Browser-safe failure information."
    failure::Union{Nothing,FailureInfo}

    "Warnings emitted by the scientific engine."
    warnings::Vector{String}
end

StructTypes.StructType(::Type{JobResult}) = StructTypes.Struct()

function validate(result::JobResult)
    result.protocol_version == PROTOCOL_VERSION || throw(ArgumentError(
        "unsupported protocol version: $(result.protocol_version)"
    ))
    result.cache_status in ("hit", "miss", "bypass") || throw(ArgumentError(
        "invalid cache status: $(result.cache_status)"
    ))
    tryparse(UUID, result.job_id) === nothing && throw(ArgumentError(
        "result job identifier must be a UUID"
    ))
    occursin(OPERATION_PATTERN, result.operation) || throw(ArgumentError(
        "invalid result operation: $(result.operation)"
    ))
    occursin(IDENTIFIER_PATTERN, result.schema_version) || throw(ArgumentError(
        "invalid result schema version"
    ))
    validate_wire_text(
        result.engine_version,
        "result engine version";
        minimum=1,
        maximum=256
    )
    occursin(r"^[0-9a-f]{64}$", result.input_hash) || throw(ArgumentError(
        "result input_hash must be a lowercase SHA-256 digest"
    ))
    occursin(r"^[0-9a-f]{64}$", result.environment_fingerprint) ||
        throw(ArgumentError("invalid worker environment fingerprint"))
    occursin(IDENTIFIER_PATTERN, result.worker_id) || throw(ArgumentError(
        "invalid result worker identifier"
    ))
    started_at = parse_utc_timestamp(result.started_at)
    completed_at = parse_utc_timestamp(result.completed_at)
    completed_at >= started_at || throw(ArgumentError(
        "result completion time precedes its start time"
    ))
    populated = count(!isnothing, (result.inline_result, result.artifact, result.failure))
    populated == 1 || throw(ArgumentError(
        "exactly one result payload must be populated"
    ))
    !isnothing(result.inline_result) && validate_wire_data(result.inline_result)
    !isnothing(result.artifact) && validate(result.artifact)
    !isnothing(result.failure) && validate(result.failure)
    length(result.warnings) <= 128 || throw(ArgumentError(
        "result contains too many warnings"
    ))
    foreach(result.warnings) do warning
        validate_wire_text(warning, "result warning"; maximum=4 * 1024)
    end
    return result
end

"""
    WorkerHeartbeat

Advertise one worker's currently available registered operations.
"""
struct WorkerHeartbeat
    "Implemented wire-contract version."
    protocol_version::String

    "Worker identity."
    worker_id::String

    "UTC heartbeat timestamp in ISO-8601 form."
    timestamp::String

    "Registered operation names."
    capabilities::Vector{String}

    "Worker version."
    worker_version::String

    "LineCableModels version or source commit."
    engine_version::String

    "Current number of jobs the worker can accept."
    available_slots::Int
end

StructTypes.StructType(::Type{WorkerHeartbeat}) = StructTypes.Struct()

function validate(heartbeat::WorkerHeartbeat)
    heartbeat.protocol_version == PROTOCOL_VERSION || throw(ArgumentError(
        "unsupported protocol version: $(heartbeat.protocol_version)"
    ))
    0 <= heartbeat.available_slots <= 1_024 || throw(ArgumentError(
        "available_slots must be between zero and 1024"
    ))
    occursin(IDENTIFIER_PATTERN, heartbeat.worker_id) || throw(ArgumentError(
        "invalid heartbeat worker identifier"
    ))
    parse_utc_timestamp(heartbeat.timestamp)
    length(heartbeat.capabilities) <= 256 || throw(ArgumentError(
        "worker capability list exceeds 256 operations"
    ))
    length(unique(heartbeat.capabilities)) == length(heartbeat.capabilities) ||
        throw(ArgumentError("worker capability list contains duplicates"))
    foreach(heartbeat.capabilities) do operation
        occursin(OPERATION_PATTERN, operation) || throw(ArgumentError(
            "invalid capability name: $operation"
        ))
    end
    validate_wire_text(
        heartbeat.worker_version,
        "worker version";
        minimum=1,
        maximum=256
    )
    validate_wire_text(
        heartbeat.engine_version,
        "worker engine version";
        minimum=1,
        maximum=256
    )
    return heartbeat
end
