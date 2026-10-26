"""
$(TYPEDEF)

Store the declared options of a formula or formulation owner `O` as a typed record.
Read option values through properties of the record.

$(TYPEDFIELDS)
"""
struct FormulationOptions{O, NT <: NamedTuple}
    "Options, each in the type of its default."
    values::NT
    function FormulationOptions{O}(values::NT) where {O, NT <: NamedTuple}
        return validate(new{O, NT}(values))
    end
end

"""
$(TYPEDSIGNATURES)

Build the options of `O` from [`defaults`](@ref) and the supplied keywords.
[`select`](@ref) converts each value to its declared type, then [`validate`](@ref)
checks the record.

# Arguments

- `O`: formula or formulation type that declares the options.

# Keywords

- `kw...`: option values that replace the defaults.

# Returns

- A `FormulationOptions` record with the declared keys in their declared order.

# Errors

- Throws `ArgumentError` for a keyword that is not an option of `O`.
"""
function FormulationOptions(::Type{O}; kw...) where {O}
    return FormulationOptions{O}(select(defaults(O), NamedTuple(kw)))
end

"""
$(TYPEDEF)

Store the declared options of a backend owner `B` for one computation.
Read option values through properties of the record.

$(TYPEDFIELDS)
"""
struct ComputationOptions{B, NT <: NamedTuple}
    "Options, each in the type of its default."
    values::NT
    function ComputationOptions{B}(values::NT) where {B, NT <: NamedTuple}
        return validate(new{B, NT}(values))
    end
end

"""
$(TYPEDSIGNATURES)

Build the options of `B` from [`defaults`](@ref) and the supplied keywords.
[`select`](@ref) converts each value to its declared type, then [`validate`](@ref)
checks the record.

# Arguments

- `B`: backend type that declares the options.

# Keywords

- `kw...`: option values that replace the defaults.

# Returns

- A `ComputationOptions` record with the declared keys in their declared order.

# Errors

- Throws `ArgumentError` for a keyword that is not an option of `B`.
"""
function ComputationOptions(::Type{B}; kw...) where {B}
    return ComputationOptions{B}(select(defaults(B), NamedTuple(kw)))
end

function Base.getproperty(
        options::Union{FormulationOptions, ComputationOptions}, name::Symbol)
    return getproperty(getfield(options, :values), name)
end
function Base.propertynames(options::Union{FormulationOptions, ComputationOptions})
    return propertynames(getfield(options, :values))
end
function Base.:(==)(a::FormulationOptions{O}, b::FormulationOptions{O}) where {O}
    return getfield(a, :values) == getfield(b, :values)
end
function Base.:(==)(a::ComputationOptions{B}, b::ComputationOptions{B}) where {B}
    return getfield(a, :values) == getfield(b, :values)
end

# The owner of a record adds the checks of its values.
validate(options::Union{FormulationOptions, ComputationOptions}) = options

# A scalar option takes the type of its default.
select(default, value) = convert(typeof(default), value)

# A complex option admits one value per frequency.
select(default::Complex, values::AbstractVector) = convert(Vector{typeof(default)}, values)

"""
$(TYPEDSIGNATURES)

Merge supplied options into their defaults, then convert each value through
[`select`](@ref).

# Arguments

- `declared`: option keys and default values.
- `supplied`: option values that replace defaults.

# Returns

- A `NamedTuple` with the declared keys, order, and value types.

# Errors

- Throws `ArgumentError` for a supplied key that `declared` lacks.
"""
function select(declared::NamedTuple, supplied::NamedTuple)
    unknown = setdiff(keys(supplied), keys(declared))
    isempty(unknown) || throw(ArgumentError(
        "unknown options $(Tuple(unknown)); the options are $(keys(declared))"))
    merged = merge(declared, supplied)
    return NamedTuple{keys(declared)}(map(select, Tuple(declared), Tuple(merged)))
end
