"""
$(TYPEDEF)

Evaluate every point in a [`ParametricProblem`](@ref) with `inner`.

`inner` may be one completed formulation, a deterministic target-bearing
formulation [`Gridspace`](@ref), or a deterministic [`Grid`](@ref) containing
completed formulations. Problem and formulation points always form a
Cartesian product. Uncertainty belongs to the problem space and is rejected
from the formulation source.

$(TYPEDFIELDS)
"""
struct Combinatorial{F, O <: ComputationOptions} <: AbstractFormulation
    "Formulation used for each core problem."
    inner::F
    "Supplemental-output retention options owned by this traversal."
    options::O
end

function computation_options(
        ::Type{Combinatorial},
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

function _combinatorial_source(inner::AbstractFormulation)
    return inner
end

function _combinatorial_source(inner::Gridspace{Target}) where {Target}
    Target <: AbstractFormulation || throw(ArgumentError(
        "Combinatorial Gridspace target must subtype AbstractFormulation; got $Target",
    ))
    has_uncertainty(inner) && throw(ArgumentError(
        "Combinatorial formulation spaces must be deterministic",
    ))
    return inner
end

function _combinatorial_source(inner::AbstractGrid)
    has_uncertainty(inner) && throw(ArgumentError(
        "Combinatorial formulation grids must be deterministic",
    ))
    all(value -> value isa AbstractFormulation, inner) || throw(ArgumentError(
        "Combinatorial formulation grids must contain completed AbstractFormulation values",
    ))
    return inner
end

function _combinatorial_source(inner)
    throw(ArgumentError(
        "Combinatorial requires a formulation, formulation Grid, or target-bearing formulation Gridspace; got $(typeof(inner))",
    ))
end

function Combinatorial(inner::F; options::NamedTuple = (;)) where {F}
    source = _combinatorial_source(inner)
    normalized = computation_options(Combinatorial, options)
    return Combinatorial{typeof(source), typeof(normalized)}(source, normalized)
end

"""
$(TYPEDEF)

Pair a lazy parameter space with computation options for a higher-order
calculation.

$(TYPEDFIELDS)
"""
struct ParametricProblem{S, O <: NamedTuple} <: AbstractProblemDefinition
    "Lazy space of complete core problems."
    space::S
    "Options supplied to each core computation."
    options::O

    function ParametricProblem(space::S, options::O) where {S, O <: NamedTuple}
        return validate(new{S, O}(space, options))
    end
end

ParametricProblem(space) = ParametricProblem(space, (;))

function validate(problem::ParametricProblem)
    problem.space isa Union{
        AbstractProblemDefinition,
        Gridspace{<:AbstractProblemDefinition}
    } || throw(ArgumentError(
        "ParametricProblem.space must be a completed problem or a Gridspace " *
        "targeting AbstractProblemDefinition; received $(typeof(problem.space))"
    ))
    if problem.space isa AbstractProblemDefinition
        validate(problem.space)
    else
        length(problem.space) > 0 || throw(ArgumentError(
            "ParametricProblem.space must contain at least one problem point"
        ))
    end
    return problem
end

"""
$(TYPEDEF)

Store core results over resolved problem and formulation axes.

Linear indexing remains compatible with ordinary result-space iteration.
Values use column-major `(problem, formulation)` order: the problem index
varies fastest. Two-dimensional indexing accepts those indices directly.

$(TYPEDFIELDS)
"""
struct ParametricResult{T, F, A <: NamedTuple, D <: ComputationDetails} <:
       AbstractParametricResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Core results in problem-index-fastest formulation order."
    values::Vector{T}
    "Problem and resolved-formulation sources aligned with `values`."
    axes::A
    "Typed supplemental output retained by the traversal."
    details::D

    function ParametricResult(
            formulation::F,
            values::Vector{T},
            axes::A,
            details::D
    ) where {T, F, A <: NamedTuple, D <: ComputationDetails}
        check_core_result(T)
        if !isempty(axes)
            keys(axes) == (:problems, :formulations) || throw(ArgumentError(
                "ParametricResult axes must contain problems and formulations",
            ))
            length(axes.problems) * length(axes.formulations) == length(values) ||
                throw(DimensionMismatch(
                    "ParametricResult axes must span every stored core result",
                ))
        end
        isempty(details) || keys(details) == (:points,) ||
            throw(ArgumentError(
                "ParametricResult details must be empty or contain only points",
            ))
        isempty(details) || length(details.points) == length(values) ||
            throw(DimensionMismatch(
                "retained details must contain one entry per core result",
            ))
        return new{T, F, A, D}(formulation, values, axes, details)
    end
end

function ParametricResult(formulation, values, details::ComputationDetails)
    ParametricResult(formulation, values, (;), details)
end
ParametricResult(formulation, values) = ParametricResult(formulation, values, (;), (;))

Base.length(value::ParametricResult) = length(value.values)
Base.size(value::ParametricResult) = (length(value),)
Base.getindex(value::ParametricResult, index::Integer) = value.values[index]
Base.iterate(value::ParametricResult, state...) = iterate(value.values, state...)
Base.firstindex(value::ParametricResult) = firstindex(value.values)
Base.lastindex(value::ParametricResult) = lastindex(value.values)

"Return one result selected by its problem and formulation indices."
function Base.getindex(
        value::ParametricResult,
        problem_index::Integer,
        formulation_index::Integer
)
    isempty(value.axes) && throw(ArgumentError(
        "this ParametricResult does not retain problem/formulation axes",
    ))
    problem_count = length(value.axes.problems)
    formulation_count = length(value.axes.formulations)
    1 <= problem_index <= problem_count || throw(BoundsError(
        value,
        (problem_index, formulation_index)
    ))
    1 <= formulation_index <= formulation_count || throw(BoundsError(
        value,
        (problem_index, formulation_index)
    ))
    index = problem_index + (formulation_index - 1) * problem_count
    return value.values[index]
end

details(value::ParametricResult) = value.details

const _GRIDSPACE_TRANSPORT_DOCUMENTATION = "https://electa-git.github.io/LineCableModels.jl/dev/gridspace/" *
                                           "#Transporting-completed-result-spaces"

function _transport_error(::Type{Target}, source::AbstractResultSpace) where {Target}
    target_name = nameof(Target)
    source_name = nameof(typeof(source))
    throw(ArgumentError(
        "Gridspace transport from $source_name to problem $target_name " *
        "is not defined. Define Gridspace{$target_name}(::$source_name) " *
        "for that result family. See `?Gridspace` and " *
        _GRIDSPACE_TRANSPORT_DOCUMENTATION,
    ))
end

"""
$(TYPEDSIGNATURES)

Transport every completed combinatorial result into one target-bearing
downstream problem point while preserving source cardinality and order.
"""
function Gridspace{Target}(
        source::ParametricResult{T, F}
) where {Target, T, F <: Combinatorial}
    return Gridspace{Target}(Target, (source,))
end

"""
$(TYPEDSIGNATURES)

Reject result-to-problem transports whose result owner has not defined how to
construct the requested target-bearing `Gridspace`.
"""
function Gridspace{Target}(source::AbstractResultSpace) where {Target}
    return _transport_error(Target, source)
end
