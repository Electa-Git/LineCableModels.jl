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

"""
$(TYPEDEF)

Record one explicitly skipped Gridspace configuration.

$(TYPEDFIELDS)
"""
struct ConfigurationFailure{C}
    "One-based configuration index."
    index::Int
    "Stable materialization coordinates."
    configuration::C
    "Concrete exception type name."
    exception_type::String
    "Rendered exception message."
    message::String
end

"""
$(TYPEDEF)

Record the stable inputs and SHA-256 identity of one composite calculation.

$(TYPEDFIELDS)
"""
struct CalculationManifest{R, P, F, S, E, O}
    "Successfully resolved parameter coordinates."
    resolved_parameterization::R
    "Original parameter-space assumptions."
    problem_assumptions::P
    "Resolved scientific formulation record."
    formulation::F
    "Primitive solver identity."
    solver::S
    "Higher-order formulation record."
    execution_policy::E
    "Primitive computation options."
    calculation_options::O
    "Stable SHA-256 digest."
    hash::String
end

function _primitive_result_type(::Type{T}) where {T}
    T <: Union{AbstractParametricResult, AbstractUncertaintyResult} && throw(ArgumentError(
        "composite results cannot contain another composite result as their primitive type",
    ))
    return T
end

const _DETAIL_KEYS = Set((:failures, :samples, :histograms, :random, :manifest))

function _validate_details(details::Dict{Symbol, NamedTuple})
    Set(keys(details)) == _DETAIL_KEYS || throw(ArgumentError(
        "result details must define exactly :failures, :samples, :histograms, " *
        ":random, and :manifest",
    ))
    return details
end

"""
$(TYPEDEF)

Store ordered primitive results from a [`Combinatorial`](@ref) calculation.

$(TYPEDFIELDS)
"""
struct ParametricResult{T, F, S, D <: Dict{Symbol, NamedTuple}} <:
       AbstractParametricResult{T}
    "Resolved higher-order formulation."
    formulation::F
    "Successful primitive results in Gridspace order."
    values::Vector{T}
    "Source parameter space."
    space::S
    "Failures, retained data, replay data, and manifest entries."
    details::D

    function ParametricResult(formulation::F, values::Vector{T}, space::S,
            details::D) where {T, F, S, D <: Dict{Symbol, NamedTuple}}
        _primitive_result_type(T)
        _validate_details(details)
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

"Return the deterministic calculation manifest stored by a composite result."
manifest(value::ParametricResult) = value.details[:manifest].value

function observables(value::ParametricResult)
    return (
        result = result(value),
        details = value.details,
        manifest = manifest(value)
    )
end
