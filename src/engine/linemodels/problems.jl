"""
$(TYPEDEF)

Define fitting of one completed propagation segment on the real frequency axis.

$(TYPEDFIELDS)
"""
struct LineModelProblem{S <: PropagationParameters} <: AbstractProblemDefinition
    "Source propagation segment."
    segment::S

    function LineModelProblem(segment::S) where {S <: PropagationParameters}
        return validate(new{S}(segment))
    end
end

function validate(problem::LineModelProblem)
    segment = problem.segment
    f = frequencies(segment)
    isempty(f) && throw(ArgumentError("line models require frequency samples"))
    all(isfinite, f) || throw(DomainError(f, "frequencies must be finite"))
    all(>(zero(eltype(f))), f) ||
        throw(DomainError(f, "line models require positive frequencies"))
    all(>(zero(eltype(f))), diff(f)) ||
        throw(ArgumentError("line models require strictly increasing frequencies"))
    parameters = segment.parameters
    all(isfinite, parameters.Z.values) ||
        throw(DomainError(parameters.Z.values, "series impedance must be finite"))
    all(isfinite, parameters.Y.values) ||
        throw(DomainError(parameters.Y.values, "shunt admittance must be finite"))
    all(isfinite, parameters.domain.gamma) ||
        throw(DomainError(parameters.domain.gamma, "propagation constants must be finite"))
    all(isfinite, parameters.domain.operators.Tv) &&
    all(isfinite, parameters.domain.operators.Ti) ||
        throw(DomainError(parameters.domain.operators, "modal operators must be finite"))
    sigma(segment) == 0 ||
        throw(ArgumentError("line models are fitted on the real frequency axis"))
    line_length(segment) > 0 ||
        throw(DomainError(line_length(segment), "line models require positive line length"))
    return problem
end

Base.eltype(problem::LineModelProblem) = eltype(problem.segment.parameters)
Base.eltype(::Type{LineModelProblem{PropagationParameters{P, L, D}}}) where {P, L, D} = eltype(fieldtype(P, :Z))
