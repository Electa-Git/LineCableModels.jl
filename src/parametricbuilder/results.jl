"""
$(TYPEDEF)

Evaluate every point in a [`ParametricProblem`](@ref) with `inner`.

$(TYPEDFIELDS)
"""
struct Combinatorial{F <: AbstractFormulation} <: AbstractFormulation
    "Formulation used for each primitive problem."
    inner::F
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
    "Options supplied to each primitive calculation."
    options::O
end

ParametricProblem(space) = ParametricProblem(space, (;))

function _primitive_result_type(::Type{T}) where {T}
    T <: Union{AbstractParametricResult, AbstractUncertaintyResult} && throw(
        ArgumentError(
        "composite results cannot contain another composite result as their primitive type",
    ),
    )
    return T
end

"""
$(TYPEDEF)

Store primitive results in Gridspace traversal order.

$(TYPEDFIELDS)
"""
struct ParametricResult{T, F} <: AbstractParametricResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Primitive results in Gridspace traversal order."
    values::Vector{T}

    function ParametricResult(formulation::F, values::Vector{T}) where {T, F}
        _primitive_result_type(T)
        return new{T, F}(formulation, values)
    end
end

Base.size(value::ParametricResult) = size(value.values)
Base.length(value::ParametricResult) = length(value.values)
Base.getindex(value::ParametricResult, index::Integer) = value.values[index]
Base.iterate(value::ParametricResult, state...) = iterate(value.values, state...)
Base.IndexStyle(::Type{<:ParametricResult}) = IndexLinear()

"Return the primitive results of a parametric calculation."
result(value::ParametricResult) = value.values
