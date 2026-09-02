"Current UTC time encoded for wire messages."
function utc_timestamp(value::Dates.DateTime=Dates.now(Dates.UTC))
    return Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z"
end

const MAX_WIRE_PAYLOAD_BYTES = 256 * 1024
const MAX_WIRE_STRING_BYTES = 128 * 1024
const MAX_WIRE_KEY_BYTES = 256
const MAX_WIRE_FIELDS = 512
const MAX_WIRE_ARRAY_ELEMENTS = 100_000
const MAX_WIRE_DEPTH = 16

function validate_wire_text(value::AbstractString, description::AbstractString;
        minimum::Integer=0,
        maximum::Integer=MAX_WIRE_STRING_BYTES
    )
    size = ncodeunits(value)
    minimum <= size <= maximum || throw(ArgumentError(
        "$description must contain between $minimum and $maximum UTF-8 bytes"
    ))
    return value
end

"""
    parse_utc_timestamp(value)

Parse a wire timestamp in canonical UTC millisecond form.
"""
function parse_utc_timestamp(value::AbstractString)
    endswith(value, "Z") || throw(ArgumentError(
        "wire timestamp must use UTC suffix Z"
    ))
    parsed = try
        DateTime(chop(value), dateformat"yyyy-mm-ddTHH:MM:SS.sss")
    catch
        throw(ArgumentError("invalid wire timestamp: $value"))
    end
    utc_timestamp(parsed) == value || throw(ArgumentError(
        "wire timestamp must use canonical millisecond precision"
    ))
    return parsed
end

function normalize_wire(value; depth::Int=0)
    depth <= MAX_WIRE_DEPTH || throw(ArgumentError(
        "wire data exceeds maximum nesting depth"
    ))
    if value isa AbstractDict || value isa JSON3.Object || value isa NamedTuple
        pairs_iter = value isa NamedTuple ? pairs(value) : pairs(value)
        length(pairs_iter) <= MAX_WIRE_FIELDS || throw(ArgumentError(
            "wire object exceeds maximum field count"
        ))
        normalized = Dict{String,Any}()
        for (key, item) in pairs_iter
            key isa Union{AbstractString,Symbol} || throw(ArgumentError(
                "wire object keys must be strings or symbols"
            ))
            normalized_key = string(key)
            validate_wire_text(
                normalized_key,
                "wire object key";
                minimum=1,
                maximum=MAX_WIRE_KEY_BYTES
            )
            haskey(normalized, normalized_key) && throw(ArgumentError(
                "wire object contains duplicate key after normalization: $normalized_key"
            ))
            normalized[normalized_key] = normalize_wire(item; depth=depth + 1)
        end
        return normalized
    elseif value isa AbstractVector || value isa Tuple || value isa JSON3.Array
        length(value) <= MAX_WIRE_ARRAY_ELEMENTS || throw(ArgumentError(
            "wire array exceeds maximum element count"
        ))
        return Any[normalize_wire(item; depth=depth + 1) for item in value]
    elseif value isa Nothing || value isa Bool || value isa AbstractString
        if value isa AbstractString
            validate_wire_text(value, "wire string")
            return string(value)
        end
        return value
    elseif value isa Integer
        return Int64(value)
    elseif value isa AbstractFloat
        isfinite(value) || throw(ArgumentError("wire numbers must be finite"))
        return Float64(value)
    elseif value isa Real
        converted = Float64(value)
        isfinite(converted) || throw(ArgumentError("wire numbers must be finite"))
        return converted
    end
    throw(ArgumentError(
        "unsupported wire value $(typeof(value)); only JSON-compatible data is allowed"
    ))
end

validate_wire_data(value) = (normalize_wire(value); value)

function canonical_json(value)
    io = IOBuffer()
    canonical_json!(io, normalize_wire(value))
    return String(take!(io))
end

function canonical_json!(io::IO, value::Dict{String,Any})
    write(io, '{')
    ordered = sort!(collect(keys(value)))
    for (index, key) in enumerate(ordered)
        index > 1 && write(io, ',')
        JSON3.write(io, key)
        write(io, ':')
        canonical_json!(io, value[key])
    end
    return write(io, '}')
end

function canonical_json!(io::IO, value::Vector{Any})
    write(io, '[')
    for (index, item) in enumerate(value)
        index > 1 && write(io, ',')
        canonical_json!(io, item)
    end
    return write(io, ']')
end

canonical_json!(io::IO, value::Bool) = JSON3.write(io, value)
canonical_json!(io::IO, value::Real) = JSON3.write(io, Float64(value))
canonical_json!(io::IO, value) = JSON3.write(io, value)

"""
    input_hash(operation, parameters; engine_constraint=nothing,
               schema_version=PROTOCOL_VERSION)

Calculate the SHA-256 digest used for idempotency and result caching.
Dictionary key ordering does not affect the digest.
"""
function input_hash(
        operation::AbstractString,
        parameters;
        engine_constraint=nothing,
        schema_version::AbstractString=PROTOCOL_VERSION
    )
    document = Dict{String,Any}(
        "operation" => string(operation),
        "schema_version" => string(schema_version),
        "engine_constraint" => engine_constraint,
        "parameters" => normalize_wire(parameters),
    )
    return bytes2hex(SHA.sha256(canonical_json(document)))
end

"""
    encode_message(message)

Validate and encode a protocol message as UTF-8 JSON. Payloads larger than
256 KiB are rejected so large results cross the artifact boundary.
"""
function encode_message(message)
    validate(message)
    payload = JSON3.write(message)
    ncodeunits(payload) <= MAX_WIRE_PAYLOAD_BYTES || throw(ArgumentError(
        "wire payload exceeds 256 KiB; store it as an artifact"
    ))
    return String(payload)
end

function decode_message(::Type{T}, payload) where {T}
    payload_size = payload isa AbstractString ? ncodeunits(payload) : length(payload)
    payload_size <= MAX_WIRE_PAYLOAD_BYTES || throw(ArgumentError(
        "wire payload exceeds 256 KiB"
    ))
    return validate(JSON3.read(payload, T))
end

decode_job_request(payload) = decode_message(JobRequest, payload)
decode_job_event(payload) = decode_message(JobEvent, payload)
decode_job_result(payload) = decode_message(JobResult, payload)
decode_worker_heartbeat(payload) = decode_message(WorkerHeartbeat, payload)

function ==(left::T, right::T) where {
        T<:Union{
            JobRequest,
            JobEvent,
            ArtifactReference,
            FailureInfo,
            JobResult,
            WorkerHeartbeat,
        }
    }
    return all(
        getfield(left, index) == getfield(right, index)
        for index in 1:fieldcount(T)
    )
end
