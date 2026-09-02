"Version implemented by the current wire contract."
const PROTOCOL_VERSION = "1.0"

const PRIORITIES = ("normal", "high")
const OPERATION_PATTERN = r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)+$"
const IDENTIFIER_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$"

"""
    JobRequest

Represent one validated request for a registered worker operation.

The `parameters` field contains JSON-compatible data. It cannot contain Julia
code, functions, expressions, or serialized objects.
"""
struct JobRequest
    "Implemented wire-contract version."
    protocol_version::String

    "Globally unique job identifier."
    job_id::String

    "Registered operation name."
    operation::String

    "Browser-session correlation identifier."
    session_id::String

    "UTC submission timestamp in ISO-8601 form."
    submitted_at::String

    "UTC deadline in ISO-8601 form."
    deadline::String

    "SHA-256 digest of the normalized scientific inputs."
    input_hash::String

    "Operation-specific JSON-compatible values."
    parameters::Dict{String,Any}

    "Optional required engine version or commit."
    engine_constraint::Union{Nothing,String}

    "Scheduling priority: `normal` or `high`."
    priority::String
end

StructTypes.StructType(::Type{JobRequest}) = StructTypes.Struct()

"""
    new_job_request(operation, parameters; session_id, timeout=Minute(5),
                    engine_constraint=nothing, priority="normal")

Construct and validate a request with a fresh UUID and content hash.

# Arguments

- `operation`: Registered dotted operation name.
- `parameters`: JSON-compatible operation inputs.

# Keywords

- `session_id`: Browser-session correlation identifier.
- `timeout`: Interval between submission and deadline.
- `engine_constraint`: Optional required engine version or commit.
- `priority`: `"normal"` or `"high"`.

# Returns

- A validated `JobRequest`.
"""
function new_job_request(
        operation::AbstractString,
        parameters;
        session_id::AbstractString,
        timeout::Dates.Period=Dates.Minute(5),
        engine_constraint::Union{Nothing,AbstractString}=nothing,
        priority::AbstractString="normal"
    )
    submitted = Dates.now(Dates.UTC)
    deadline = submitted + timeout
    deadline > submitted || throw(ArgumentError("job timeout must be positive"))
    normalized = normalize_wire(parameters)
    request = JobRequest(
        PROTOCOL_VERSION,
        string(uuid4()),
        string(operation),
        string(session_id),
        utc_timestamp(submitted),
        utc_timestamp(deadline),
        input_hash(operation, normalized; engine_constraint),
        normalized,
        isnothing(engine_constraint) ? nothing : string(engine_constraint),
        string(priority)
    )
    return validate(request)
end

function validate(request::JobRequest)
    request.protocol_version == PROTOCOL_VERSION || throw(ArgumentError(
        "unsupported protocol version: $(request.protocol_version)"
    ))
    occursin(IDENTIFIER_PATTERN, request.job_id) || throw(ArgumentError(
        "invalid job identifier"
    ))
    tryparse(UUID, request.job_id) === nothing && throw(ArgumentError(
        "job identifier must be a UUID"
    ))
    occursin(OPERATION_PATTERN, request.operation) || throw(ArgumentError(
        "invalid operation name: $(request.operation)"
    ))
    ncodeunits(request.operation) <= 128 || throw(ArgumentError(
        "operation name cannot exceed 128 bytes"
    ))
    occursin(IDENTIFIER_PATTERN, request.session_id) || throw(ArgumentError(
        "invalid session identifier"
    ))
    submitted_at = parse_utc_timestamp(request.submitted_at)
    deadline = parse_utc_timestamp(request.deadline)
    deadline > submitted_at || throw(ArgumentError(
        "job deadline must follow submission time"
    ))
    deadline - submitted_at <= Dates.Hour(24) || throw(ArgumentError(
        "job deadline cannot exceed 24 hours after submission"
    ))
    request.priority in PRIORITIES || throw(ArgumentError(
        "priority must be normal or high"
    ))
    occursin(r"^[0-9a-f]{64}$", request.input_hash) || throw(ArgumentError(
        "input_hash must be a lowercase SHA-256 digest"
    ))
    validate_wire_data(request.parameters)
    expected_hash = input_hash(
        request.operation,
        request.parameters;
        engine_constraint=request.engine_constraint
    )
    request.input_hash == expected_hash || throw(ArgumentError(
        "input_hash does not match normalized job inputs"
    ))
    if !isnothing(request.engine_constraint)
        1 <= ncodeunits(request.engine_constraint) <= 256 || throw(ArgumentError(
            "engine_constraint must contain between 1 and 256 bytes"
        ))
    end
    return request
end

"""
    operation_subject_token(operation)

Convert a registered dotted operation name into one NATS subject token.
"""
function operation_subject_token(operation::AbstractString)
    occursin(OPERATION_PATTERN, operation) || throw(ArgumentError(
        "invalid operation name: $operation"
    ))
    ncodeunits(operation) <= 128 || throw(ArgumentError(
        "operation name cannot exceed 128 bytes"
    ))
    return replace(string(operation), '.' => '_')
end
