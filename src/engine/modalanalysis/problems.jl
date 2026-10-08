"""
$(TYPEDEF)

Define modal analysis of one completed phase-domain frequency scan.

$(TYPEDFIELDS)
"""
struct ModalAnalysisProblem{P <: LineParameters} <: AbstractProblemDefinition
    "Source line parameters."
    parameters::P

    function ModalAnalysisProblem{P}(parameters::P) where {P <: LineParameters}
        return validate(new{P}(parameters))
    end
end

function validate(problem::ModalAnalysisProblem)
    validate(problem.parameters)
    parameters = problem.parameters
    size(parameters.Z, 1) > 0 || throw(ArgumentError("modal analysis requires at least one mode"))
    isempty(parameters.f) && throw(ArgumentError("modal analysis requires frequency samples"))
    all(>(zero(eltype(parameters.f))), parameters.f) ||
        throw(DomainError(parameters.f, "modal analysis requires positive frequencies"))
    all(>(zero(eltype(parameters.f))), diff(parameters.f)) ||
        throw(ArgumentError("modal analysis requires strictly increasing frequencies"))
    all(isfinite, parameters.Z.values) ||
        throw(DomainError(parameters.Z.values, "series impedance must be finite"))
    all(isfinite, parameters.Y.values) ||
        throw(DomainError(parameters.Y.values, "shunt admittance must be finite"))
    return problem
end

function ModalAnalysisProblem(
        parameters::LineParameters{T, U, PhaseDomain}
) where {T, U}
    return ModalAnalysisProblem{typeof(parameters)}(parameters)
end

Base.eltype(problem::ModalAnalysisProblem) = eltype(problem.parameters)
Base.eltype(::Type{ModalAnalysisProblem{P}}) where {P} = eltype(P)
