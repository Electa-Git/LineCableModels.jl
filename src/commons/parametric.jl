"""
$(TYPEDEF)

Evaluate every point in a [`ParametricProblem`](@ref) with `inner`.

`inner` is one core formulation, a deterministic [`Gridspace`](@ref) of them or a
deterministic [`Grid`](@ref) of them, a `ScalarOrGrid{AbstractFormulation}`.
Problem and formulation points always form a Cartesian product. Uncertainty belongs to the
problem space and is rejected from the formulation source. The study has no options of its
own: the computation options of its run are those of the parametric call.

$(TYPEDFIELDS)
"""
struct Combinatorial{F} <: AbstractTraversal
    "Formulation used for each core problem."
    inner::F

    Combinatorial(inner::F) where {F <: ScalarOrGrid{AbstractFormulation}} =
        validate(new{F}(inner))
end

# The computation options of the run of a study, whatever the study, are `verbosity`,
# `on_result` and `retain_details`.
function computation_options(
        ::Type{<:AbstractTraversal},
        record::ComputationOptions
)::ComputationOptions
    options = record.data
    defaults = (verbosity = (default = 0,), on_result = nothing, retain_details = false)
    unknown = filter(key -> !haskey(defaults, key), keys(options))
    isempty(unknown) || throw(ArgumentError(
        "unknown computation options of a parametric run: $(sort!(collect(unknown))). " *
        "The computation options of the backend for every point are the options of " *
        "ParametricProblem(space, options)",
    ))
    normalized = merge(defaults, options)
    normalized.retain_details isa Bool || throw(ArgumentError("retain_details must be Bool"))
    return ComputationOptions(; verbosity = verbosity(normalized.verbosity),
        normalized.on_result, normalized.retain_details)
end

# A study stores a valid inner formulation.
function validate(study::T) where {T <: AbstractTraversal}
    validate(study.inner, T)
    return study
end

# The inner formulation of a study: one core formulation, or, for a Combinatorial study, a
# deterministic space or Grid of them.
validate(formulation::AbstractFormulation, ::Type{<:AbstractTraversal}) = formulation
function validate(inner::Gridspace, ::Type{<:Combinatorial})
    has_uncertainty(inner) && throw(ArgumentError(
        "Combinatorial formulation spaces must be deterministic",
    ))
    return inner
end
function validate(inner::AbstractGrid, ::Type{<:Combinatorial})
    has_uncertainty(inner) && throw(ArgumentError(
        "Combinatorial formulation grids must be deterministic",
    ))
    return inner
end

# A Combinatorial study's formulations are the points of its source.
points(study::Combinatorial) = points(study.inner)

"""
$(TYPEDEF)

Pair a core problem or a lazy parameter space of core problems with the computation
options of the backend for every point. The owner of the formulations of a study normalizes
these options once, before the first point. A scalar problem is a space of one point, and
the problem itself need not implement iteration. As a problem space, an
[`AbstractProblemSpace`](@ref), a parametric problem contains core problems, never a
parametric problem. Its options never vary, and the computation options of the run are those
of the parametric call, `compute(problem, study; options)`.

$(TYPEDFIELDS)
"""
struct ParametricProblem{
    S <: ScalarOrGrid{AbstractProblemDefinition},
    O <: ComputationOptions
} <: AbstractProblemSpace
    "A complete core problem or a lazy space of them."
    space::S
    "Computation options of the backend for every point, which its owner normalizes."
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

Store the results of a parametric computation and the failed points of its space. A result
does not contain a problem space, a formulation space or a run control.

The values follow column-major `(problem, formulation)` order: the problem index varies
fastest. Linear indexing and iteration run over the values of the computed points. A failed
point has no value: `failures` holds its gridpoint identity and the exception it threw, in
the same order. `gridpoints` gives the problem and formulation indices of every point,
computed or failed, and two-dimensional indexing accepts those indices. For a failed point
it throws an error with the indices of the point and its exception. Each value of a space
computed by [`traverse`](@ref) also records its gridpoint identity in its own details.

$(TYPEDFIELDS)
"""
struct ParametricResult{T} <: AbstractParametricResult{T}
    "Core results of the computed points in problem-index-fastest formulation order."
    values::Vector{T}
    "Gridpoint identity and exception of each failed point, in the same order."
    failures::Vector{Pair{typeof(gridpoint_id()), Exception}}
    "Problem and formulation indices of the points of the space, computed or failed."
    gridpoints::CartesianIndices{2, NTuple{2, Base.OneTo{Int}}}
    "Supplemental output retained by the traversal."
    details::ComputationDetails

    function ParametricResult(values::Vector{T},
            failures::Vector{Pair{typeof(gridpoint_id()), Exception}},
            gridpoints::CartesianIndices{2, NTuple{2, Base.OneTo{Int}}},
            details::ComputationDetails) where {T}
        validate(T, ParametricResult)
        length(values) + length(failures) == length(gridpoints) || throw(DimensionMismatch(
            "ParametricResult gridpoints must span every computed and failed point",
        ))
        isempty(details.data) || keys(details.data) == (:points,) ||
            throw(ArgumentError(
                "ParametricResult details must be empty or contain only points",
            ))
        isempty(details.data) || length(details.data.points) == length(values) ||
            throw(DimensionMismatch(
                "retained details must contain one entry per core result",
            ))
        return new{T}(values, failures, gridpoints, details)
    end
end

Base.length(value::ParametricResult) = length(value.values)
Base.size(value::ParametricResult) = (length(value),)
Base.getindex(value::ParametricResult, index::Integer) = value.values[index]
Base.iterate(value::ParametricResult, state...) = iterate(value.values, state...)
Base.firstindex(value::ParametricResult) = firstindex(value.values)
Base.lastindex(value::ParametricResult) = lastindex(value.values)

"""
Return the result of a point selected by its problem and formulation indices. For a failed
point, throw an error with the indices of the point and the exception that it threw.
"""
function Base.getindex(value::ParametricResult, problem_index::Integer,
        formulation_index::Integer)
    index = LinearIndices(value.gridpoints)[problem_index, formulation_index]
    # The values leave out the failed points before the point.
    preceding = 0
    for (gridpoint, exception) in value.failures
        failed = LinearIndices(value.gridpoints)[gridpoint.problem_index,
            gridpoint.formulation_index]
        failed == index && throw(ErrorException("the point of problem $problem_index " *
            "and formulation $formulation_index failed: " * sprint(showerror, exception)))
        preceding += failed < index
    end
    return value.values[index - preceding]
end

details(value::ParametricResult) = value.details

"""
$(TYPEDSIGNATURES)

Transport every completed parametric result into one target-bearing downstream problem
point while preserving source cardinality and order.
"""
Gridspace{Target}(source::ParametricResult) where {Target} =
    Gridspace{Target}(Target, (source,))
