"""
$(TYPEDEF)

Evaluate every point in a [`ParametricProblem`](@ref) with `inner`.

`inner` may be one completed formulation, a deterministic target-bearing
formulation [`Gridspace`](@ref) or a deterministic [`Grid`](@ref) containing
completed formulations. Problem and formulation points always form a
Cartesian product. Uncertainty belongs to the problem space and is rejected
from the formulation source.

$(TYPEDFIELDS)
"""
struct Combinatorial{F, O <: ComputationOptions} <: AbstractTraversal
    "Formulation used for each core problem."
    inner::F
    "Supplemental-output retention options owned by this traversal."
    options::O

    Combinatorial(inner::F, options::O) where {
        F <: Union{AbstractFormulation, AbstractGrid, Gridspace}, O <: ComputationOptions} =
        validate(new{F, O}(inner, options))
end

function (::Type{T})(inner;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions(),
        keywords...) where {T <: AbstractTraversal}
    record = ComputationOptions(options)
    duplicates = filter(key -> haskey(record.data, key), keys(keywords))
    isempty(duplicates) || throw(ArgumentError(
        "$(nameof(T)) computation options supplied both as keywords and in options: " *
        "$(collect(duplicates))",
    ))
    merged = ComputationOptions(merge(record.data, values(keywords)))
    return T(inner, computation_options(T, merged))
end

# A traversal retains the details of its points or not.
function computation_options(
        ::Type{T},
        record::ComputationOptions
)::ComputationOptions where {T <: AbstractTraversal}
    options = record.data
    unknown = filter(key -> key !== :retain_details, keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown $(nameof(T)) computation options: $(sort!(collect(unknown)))",
    ))
    normalized = merge((retain_details = false,), options)
    normalized.retain_details isa Bool || throw(ArgumentError(
        "$(nameof(T)) retain_details must be Bool",
    ))
    return ComputationOptions(; retain_details = normalized.retain_details)
end

# A traversal stores a valid inner formulation and normalized options.
function validate(formulation::T) where {T <: AbstractTraversal}
    validate(formulation.inner, T)
    computation_options(T, formulation.options) === formulation.options ||
        throw(ArgumentError(
            "$(nameof(T)) options must be normalized; construct it with keywords",
        ))
    return formulation
end

# The inner formulation of a traversal: one completed formulation, or, for a Combinatorial
# computation, a deterministic space or Grid of them.
validate(formulation::AbstractFormulation, ::Type{<:AbstractTraversal}) = formulation
function validate(inner::Gridspace{Target}, ::Type{<:Combinatorial}) where {Target}
    Target <: AbstractFormulation || throw(ArgumentError(
        "Combinatorial Gridspace target must subtype AbstractFormulation; got $Target",
    ))
    has_uncertainty(inner) && throw(ArgumentError(
        "Combinatorial formulation spaces must be deterministic",
    ))
    return inner
end
function validate(inner::AbstractGrid, ::Type{<:Combinatorial})
    has_uncertainty(inner) && throw(ArgumentError(
        "Combinatorial formulation grids must be deterministic",
    ))
    foreach(value -> validate(value, Combinatorial), inner)
    return inner
end
function validate(value, ::Type{T}) where {T <: AbstractTraversal}
    throw(ArgumentError(
        "$(nameof(T)) formulation sources must contain completed AbstractFormulation " *
        "values; got $(nameof(typeof(value)))"
    ))
end

# A Combinatorial formulation is its own record of the traversal.
Combinatorial(formulation::Combinatorial) = formulation

# A Combinatorial formulation's points are those of its source, and it varies when its
# source does.
points(formulation::Combinatorial) = points(formulation.inner)
has_grid(formulation::Combinatorial) = has_grid(formulation.inner)

"""
$(TYPEDEF)

Pair a problem or a lazy parameter space of problems with the computation options of a
computation. A scalar problem stays a scalar, a run of one point. The problem itself need
not implement iteration.

$(TYPEDFIELDS)
"""
struct ParametricProblem{
    S <: ScalarOrGrid{AbstractProblemDefinition},
    O <: ComputationOptions
} <: AbstractProblemDefinition
    "A complete core problem or a lazy space of them."
    space::S
    "Options forwarded to each core computation for owner-specific normalization."
    options::O

    function ParametricProblem(space::S, options::O) where {
            S <: ScalarOrGrid{AbstractProblemDefinition},
            O <: ComputationOptions
    }
        return validate(new{S, O}(space, options))
    end
end

function ParametricProblem(space::ScalarOrGrid{AbstractProblemDefinition})
    return ParametricProblem(space, ComputationOptions())
end

function validate(problem::ParametricProblem)
    isempty(points(problem.space)) && throw(ArgumentError(
        "ParametricProblem.space must contain at least one problem point"
    ))
    return problem
end

"""
$(TYPEDEF)

Store the results of a parametric computation: the computed values, the gridpoints of the
space that computed them and the retained details. A result does not contain a problem
space, a formulation space or a run control.

The values follow column-major `(problem, formulation)` order: the problem index varies
fastest. Linear indexing and iteration run over the values. `gridpoints` gives the problem
and formulation indices of each point, and two-dimensional indexing accepts those indices.
Each value of a space computed by [`traverse`](@ref) also records its gridpoint identity in
its own details.

$(TYPEDFIELDS)
"""
struct ParametricResult{T} <: AbstractParametricResult{T}
    "Core results in problem-index-fastest formulation order."
    values::Vector{T}
    "Problem and formulation indices of the points of the space, in the order of `values`."
    gridpoints::CartesianIndices{2, NTuple{2, Base.OneTo{Int}}}
    "Supplemental output retained by the traversal."
    details::ComputationDetails

    function ParametricResult(values::Vector{T},
            gridpoints::CartesianIndices{2, NTuple{2, Base.OneTo{Int}}},
            details::ComputationDetails) where {T}
        validate(T, ParametricResult)
        length(gridpoints) == length(values) || throw(DimensionMismatch(
            "ParametricResult gridpoints must span every stored core result",
        ))
        isempty(details.data) || keys(details.data) == (:points,) ||
            throw(ArgumentError(
                "ParametricResult details must be empty or contain only points",
            ))
        isempty(details.data) || length(details.data.points) == length(values) ||
            throw(DimensionMismatch(
                "retained details must contain one entry per core result",
            ))
        return new{T}(values, gridpoints, details)
    end
end

Base.length(value::ParametricResult) = length(value.values)
Base.size(value::ParametricResult) = (length(value),)
Base.getindex(value::ParametricResult, index::Integer) = value.values[index]
Base.iterate(value::ParametricResult, state...) = iterate(value.values, state...)
Base.firstindex(value::ParametricResult) = firstindex(value.values)
Base.lastindex(value::ParametricResult) = lastindex(value.values)

"""
Return one result selected by its problem and formulation indices.
"""
function Base.getindex(value::ParametricResult, problem_index::Integer,
        formulation_index::Integer)
    return value.values[LinearIndices(value.gridpoints)[problem_index, formulation_index]]
end

details(value::ParametricResult) = value.details

"""
$(TYPEDSIGNATURES)

Transport every completed parametric result into one target-bearing downstream problem
point while preserving source cardinality and order.
"""
Gridspace{Target}(source::ParametricResult) where {Target} =
    Gridspace{Target}(Target, (source,))
