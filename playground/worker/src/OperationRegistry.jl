struct PermanentOperationError <: Exception
    category::String
    message::String
end

Base.showerror(io::IO, error::PermanentOperationError) = print(io, error.message)

struct RetryableOperationError <: Exception
    message::String
end

Base.showerror(io::IO, error::RetryableOperationError) = print(io, error.message)

"""
    OperationSpec

Define one closed worker operation, including validation, execution mode,
resource limit, result schema, capability name, and cache policy.
"""
struct OperationSpec
    "Registered dotted operation name."
    name::String

    "Operation result-schema version."
    schema_version::String

    "Input validation function."
    validator::Function

    "Projection of request parameters onto scientifically cache-relevant inputs."
    cache_key_inputs::Function

    "Execution function accepting `(context, parameters)`."
    executor::Function

    "Maximum execution time [s]."
    timeout_seconds::Float64

    "Advertised capability name."
    capability::String

    "Cache policy: `content`, `prepared`, or `none`."
    cache_policy::Symbol

    "Execution mode: `daemon` or `supervised`."
    execution_mode::Symbol
end

function OperationSpec(
        name::AbstractString,
        validator,
        executor;
        schema_version::AbstractString="1.0",
        cache_key_inputs::Function=identity,
        timeout_seconds::Real=60.0,
        capability::AbstractString=name,
        cache_policy::Symbol=:content,
        execution_mode::Symbol=:daemon
    )
    operation_subject_token(name)
    operation_subject_token(capability)
    occursin(r"^[0-9]+(?:\.[0-9]+){1,2}$", schema_version) || throw(ArgumentError(
        "operation schema version must contain two or three numeric components"
    ))
    ncodeunits(schema_version) <= 32 || throw(ArgumentError(
        "operation schema version cannot exceed 32 bytes"
    ))
    0 < timeout_seconds <= 24 * 60 * 60 || throw(ArgumentError(
        "operation timeout must be between zero and 24 hours"
    ))
    cache_policy in (:content, :prepared, :none) || throw(ArgumentError(
        "unsupported cache policy: $cache_policy"
    ))
    execution_mode in (:daemon, :supervised) || throw(ArgumentError(
        "unsupported execution mode: $execution_mode"
    ))
    return OperationSpec(
        string(name),
        string(schema_version),
        validator,
        cache_key_inputs,
        executor,
        Float64(timeout_seconds),
        string(capability),
        cache_policy,
        execution_mode
    )
end

mutable struct OperationRegistry
    operations::Dict{String,OperationSpec}
end

OperationRegistry() = OperationRegistry(Dict{String,OperationSpec}())

function register!(registry::OperationRegistry, spec::OperationSpec)
    haskey(registry.operations, spec.name) && throw(ArgumentError(
        "operation already registered: $(spec.name)"
    ))
    registry.operations[spec.name] = spec
    return registry
end

function registered_operation(registry::OperationRegistry, name::AbstractString)
    return get(registry.operations, string(name)) do
        throw(PermanentOperationError(
            "unknown_operation",
            "Unknown registered operation: $name"
        ))
    end
end

capabilities(registry::OperationRegistry) = sort!(collect(keys(registry.operations)))

function required(parameters::AbstractDict, key::String, ::Type{T}) where {T}
    haskey(parameters, key) || throw(PermanentOperationError(
        "invalid_input",
        "Missing required parameter: $key"
    ))
    value = parameters[key]
    value isa T || throw(PermanentOperationError(
        "invalid_input",
        "Parameter $key must be $(T)"
    ))
    return value
end

function bounded_real(parameters, key, lower, upper)
    value = required(parameters, key, Real)
    lower <= value <= upper || throw(PermanentOperationError(
        "invalid_input",
        "Parameter $key must be between $lower and $upper"
    ))
    return Float64(value)
end

function bounded_integer(parameters, key, lower, upper)
    value = required(parameters, key, Integer)
    lower <= value <= upper || throw(PermanentOperationError(
        "invalid_input",
        "Parameter $key must be between $lower and $upper"
    ))
    return Int(value)
end
