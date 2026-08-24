"""
$(TYPEDEF)

Evaluate every admitted configuration with `inner` while applying the selected
invalid-configuration policy.

$(TYPEDFIELDS)
"""
struct Combinatorial{F <: AbstractFormulation} <: AbstractFormulation
    "Formulation used for each materialized problem."
    inner::F
    "Invalid-configuration policy: `:error` or `:skip`."
    invalid::Symbol

    function Combinatorial(inner::F; invalid::Symbol = :error) where {F <:
                                                                      AbstractFormulation}
        invalid in (:error, :skip) || throw(ArgumentError(
            "invalid must be :error or :skip; got $(repr(invalid))",
        ))
        return new{F}(inner, invalid)
    end
end

"""
$(TYPEDEF)

Pair a lazy parameter space with computation options for a higher-order
calculation.

$(TYPEDFIELDS)
"""
struct ParametricProblem{S, O <: NamedTuple} <: AbstractProblemDefinition
    "Lazy space of complete primitive problems."
    space::S
    "Options forwarded to each primitive calculation."
    options::O
end

ParametricProblem(space) = ParametricProblem(space, (;))

function _primitive_result_type(::Type{T}) where {T}
    T <: Union{AbstractParametricResult, AbstractUncertaintyResult} && throw(ArgumentError(
        "composite results cannot contain another composite result as their primitive type",
    ))
    return T
end

"""
$(TYPEDEF)

Store ordered primitive results from a [`Combinatorial`](@ref) calculation.

$(TYPEDFIELDS)
"""
struct ParametricResult{
    T, F, S <: AbstractVector{<:NamedTuple}, D <: Dict{Symbol, NamedTuple}} <:
       AbstractParametricResult{T}
    "Resolved higher-order formulation."
    formulation::F
    "Successful primitive results in Gridspace order."
    values::Vector{T}
    "Successful resolved coordinates, aligned with `values`."
    space::S
    "Contingent failure and replay metadata."
    details::D

    function ParametricResult(formulation::F, values::Vector{T}, space::S,
            details::D) where {
            T, F, S <: AbstractVector{<:NamedTuple}, D <: Dict{Symbol, NamedTuple}}
        _primitive_result_type(T)
        length(space) == length(values) || throw(DimensionMismatch(
            "resolved coordinates must contain one entry per primitive result",
        ))
        return new{T, F, S, D}(formulation, values, space, details)
    end
end

Base.size(value::ParametricResult) = size(value.values)
Base.length(value::ParametricResult) = length(value.values)
Base.getindex(value::ParametricResult, index::Integer) = value.values[index]
Base.iterate(value::ParametricResult, state...) = iterate(value.values, state...)
Base.IndexStyle(::Type{<:ParametricResult}) = IndexLinear()

"Return the ordered primitive results of a composite calculation."
result(value::ParametricResult) = value.values

function observables(value::ParametricResult)
    return (
        result = result(value),
        details = value.details
    )
end
