import ...DataModel: terminal_labels

"""
$(TYPEDEF)

Define transient analysis of an uncomputed line or fitted model with terminal assignments.
An end is `Matched()` or a tuple of assignments. Unlisted terminals are open.

$(TYPEDFIELDS)
"""
struct TransientAnalysisProblem{L <: Union{LineParametersProblem, LineModel}, S, R, W <: NamedTuple} <: AbstractProblemDefinition
    "Uncomputed line problem or fitted line model."
    line::L
    "Sending-end assignments."
    sending::S
    "Receiving-end assignments."
    receiving::R
    "Time window `(tmax, step)` \\[s\\]."
    window::W

    function TransientAnalysisProblem(line::L, sending::S, receiving::R, window::W) where {
            L <: Union{LineParametersProblem, LineModel}, S, R, W <: NamedTuple}
        return validate(new{L, S, R, W}(line, sending, receiving, window))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a transient problem, lifting grids in the line, terminal assignments or window.
"""
function TransientAnalysisProblem(line; sending, receiving, window)
    return parameterize(TransientAnalysisProblem, TransientAnalysisProblem,
        (line, sending, receiving, window))
end

validate(assignment::Matched, ::Type{<:TransientAnalysisProblem}, coordinates) = assignment

function validate(assignments::Tuple, ::Type{<:TransientAnalysisProblem}, coordinates)
    for assignment in assignments
        validate(assignment, TransientAnalysisProblem, coordinates)
    end
    allunique(map(assignment -> assignment.terminal, assignments)) ||
        throw(ArgumentError("a terminal may be assigned at most once per end"))
    return assignments
end

function validate(assignment::Union{TerminalImpedance, Waveform},
        ::Type{<:TransientAnalysisProblem}, coordinates)
    assignment.terminal in coordinates ||
        throw(ArgumentError("unknown terminal label $(assignment.terminal)"))
    return assignment
end

function validate(problem::TransientAnalysisProblem, line::LineParametersProblem)
    coordinates = terminal_labels(line.system)
    validate(problem.sending, TransientAnalysisProblem, coordinates)
    validate(problem.receiving, TransientAnalysisProblem, coordinates)
    return problem
end

function validate(problem::TransientAnalysisProblem, line::LineModel)
    validate(problem.sending, TransientAnalysisProblem, line.coordinates)
    validate(problem.receiving, TransientAnalysisProblem, line.coordinates)
    return problem
end

function validate(problem::TransientAnalysisProblem)
    validate(problem, problem.line)
    keys(problem.window) == (:tmax, :step) ||
        throw(ArgumentError("window must contain tmax and step"))
    window = problem.window
    isfinite(window.tmax) && window.tmax > 0 ||
        throw(DomainError(window.tmax, "tmax must be finite and positive"))
    isfinite(window.step) && 0 < window.step < window.tmax ||
        throw(DomainError(window.step, "step must be finite, positive and less than tmax"))
    return problem
end
