"""
    LineCableModelsPlaygroundProtocol

Define the versioned, data-only contract shared by the playground publisher
and scientific workers. The package contains no engine implementation.
"""
module LineCableModelsPlaygroundProtocol

using Dates
using JSON3
using SHA
using StructTypes
using UUIDs

import Base: ==

include("Jobs.jl")
include("Events.jl")
include("Encoding.jl")

export ArtifactReference,
    FailureInfo,
    JobEvent,
    JobRequest,
    JobResult,
    WorkerHeartbeat,
    PROTOCOL_VERSION,
    PRIORITIES,
    decode_job_event,
    decode_job_request,
    decode_job_result,
    decode_worker_heartbeat,
    encode_message,
    input_hash,
    isterminal,
    new_job_request,
    normalize_wire,
    operation_subject_token,
    parse_utc_timestamp,
    utc_timestamp,
    validate

end
