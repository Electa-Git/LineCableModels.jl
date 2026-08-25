"""
$(TYPEDEF)

Evaluate every point in a [`ParametricProblem`](@ref) with `inner`.

$(TYPEDFIELDS)
"""
struct Combinatorial{F <: AbstractFormulation, O <: ComputationOptions} <: AbstractFormulation
    "Formulation used for each primitive problem."
    inner::F
    "Supplemental-output retention options owned by this traversal."
    options::O
end

function computation_options(
        ::Val{Combinatorial},
        options::NamedTuple
)::ComputationOptions
    unknown = filter(key -> key !== :retain_details, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown Combinatorial computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((retain_details = false,), options)
    normalized.retain_details isa Bool || throw(ArgumentError(
        "Combinatorial retain_details must be Bool",
    ))
    return (retain_details = normalized.retain_details,)
end

function Combinatorial(
        inner::F;
        options::NamedTuple = (;)
) where {F <: AbstractFormulation}
    normalized = computation_options(Val(Combinatorial), options)
    return Combinatorial{F, typeof(normalized)}(inner, normalized)
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
struct ParametricResult{T, F, D <: ComputationDetails} <: AbstractParametricResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Primitive results in Gridspace traversal order."
    values::Vector{T}
    "Typed supplemental output retained by the traversal."
    details::D

    function ParametricResult(
            formulation::F,
            values::Vector{T},
            details::D
    ) where {T, F, D <: ComputationDetails}
        _primitive_result_type(T)
        isempty(details) || keys(details) == (:points,) || throw(ArgumentError(
            "ParametricResult details must be empty or contain only points",
        ))
        isempty(details) || length(details.points) == length(values) ||
            throw(DimensionMismatch(
                "retained details must contain one entry per primitive result",
            ))
        return new{T, F, D}(formulation, values, details)
    end
end

ParametricResult(formulation, values) = ParametricResult(formulation, values, (;))

Base.size(value::ParametricResult) = size(value.values)
Base.length(value::ParametricResult) = length(value.values)
Base.getindex(value::ParametricResult, index::Integer) = value.values[index]
Base.iterate(value::ParametricResult, state...) = iterate(value.values, state...)
Base.IndexStyle(::Type{<:ParametricResult}) = IndexLinear()

"Return the primitive results of a parametric calculation."
result(value::ParametricResult) = value.values
details(value::ParametricResult) = value.details
