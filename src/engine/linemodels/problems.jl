"""
$(TYPEDEF)

Define fitting of one completed propagation segment on the real frequency axis.

$(TYPEDFIELDS)
"""
struct LineModelProblem{S <: PropagationParameters, L <: Union{Nothing, LineParametersProblem}} <: AbstractProblemDefinition
    "Source propagation segment."
    segment::S
    "Physical line problem used to compute the segment, when supplied."
    line::L

    function LineModelProblem(segment::S;
            line::Union{Nothing, LineParametersProblem}=nothing) where {S <: PropagationParameters}
        return validate(new{S, typeof(line)}(segment, line))
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
    validate(problem.line, problem)
    return problem
end

validate(line::Nothing, ::LineModelProblem) = line

function validate(line::LineParametersProblem, problem::LineModelProblem)
    iszero(line.sigma) || throw(ArgumentError(
        "the line problem must be on the real frequency axis"))
    line.frequencies == frequencies(problem.segment) || throw(ArgumentError(
        "the line problem frequencies must equal the segment frequencies"))
    return line
end

Base.eltype(problem::LineModelProblem) = eltype(problem.segment.parameters)
Base.eltype(::Type{LineModelProblem{PropagationParameters{P, L, D}, I}}) where {P, L, D, I} = eltype(fieldtype(P, :Z))

function validate(subject::NamedTuple, problem::LineModelProblem, expected::NamedTuple,
        path::Tuple)
    keys(subject) == keys(expected) || throw(ArgumentError(
        "line consistency mismatch at $(join(path, '.'))"))
    for name in keys(subject)
        validate(subject[name], problem, expected[name], (path..., name))
    end
    return subject
end

function validate(subject, problem::LineModelProblem, expected, path::Tuple)
    isequal(subject, expected) || throw(ArgumentError(
        "line consistency mismatch at $(join(path, '.'))"))
    return subject
end
