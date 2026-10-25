"""
$(TYPEDEF)

Store the options of a formula or a formulation. The owner `O`, a formula or formulation
type, selects the defaults and the checks. `NT` is the type of the options, each in its
declared order and type. The records of one owner have one type, whichever keywords built
them. An option reads as a property, such as `options.bundle`.

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

Build the options of the owner `O` from [`defaults`](@ref) of `O` and the keywords given.
Each value takes the type of its default through [`select`](@ref), and `validate` checks the
record.

# Errors

- Throws `ArgumentError` for a keyword that is not an option of `O`.
"""
function FormulationOptions(::Type{O}; kw...) where {O}
    return FormulationOptions{O}(select(defaults(O), NamedTuple(kw)))
end

"""
$(TYPEDEF)

Store the options of a backend for one computation, such as whether it keeps the
intermediate matrices. The owner `B` is the backend type, and the record has the form of
[`FormulationOptions`](@ref).

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

Build the options of the backend `B` from [`defaults`](@ref) of `B` and the keywords given.
Each value takes the type of its default through [`select`](@ref), and `validate` checks the
record.

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

Fit the `supplied` options to the `declared` ones. The result has the declared keys in their
order. A supplied value replaces its default, and [`select`](@ref) fits each value to the
type of its default.

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
