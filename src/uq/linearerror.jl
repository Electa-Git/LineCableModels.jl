function compute(problem::ParametricProblem, formulation::LinearError;
        modal=nothing, modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    modal === nothing || return compute(problem, formulation,
        LineCableModels.ModalAnalysisFormulation(modal); modal_options)
    isempty(ComputationOptions(modal_options).data) ||
        throw(ArgumentError("modal_options require a modal formulation"))
    levels = verbosity(get(problem.options.data, :verbosity, (default = 0,)))
    progress = problem.space isa
               LineCableModels.Gridspace{<:Engine.LineParametersProblem} &&
               get(levels, :progress, levels.default) > 0
    logger = haskey(problem.options.data, :verbosity) ?
             VerbosityLogger(Logging.current_logger(), levels) : Logging.current_logger()
    return Logging.with_logger(logger) do
        started = progress ? time_ns() : UInt64(0)
        progress &&
            @info "Linear error computation started" _group=:progress problems=length(problem.space)
        result = _linear_error(formulation.derivative, formulation, problem)
        progress &&
            @info "Linear error computation completed successfully" _group=:progress completed=length(result) total=length(result) elapsed_seconds=(time_ns()-started)*1e-9
        result
    end
end

function compute(problem::ParametricProblem, formulation::LinearError,
        modal::LineCableModels.ModalAnalysisFormulation;
        modal_options::Union{NamedTuple,ComputationOptions}=ComputationOptions())
    return _linear_error(formulation.derivative, formulation,
        problem, modal, ComputationOptions(modal_options))
end

function _linear_error(::MeasurementsPropagation, formulation::LinearError, problem)
    traversed = traverse(problem, formulation)
    values = map(traversed.values) do value
        id = get(Commons.observation_gridpoint(value), :id, nothing)
        Engine.retain_gridpoint(value, id;
            fields = (uncertainty = (
                estimator = :first_order, representation = :dependency_preserving),
                uncertainty_descriptions = _uncertainty_descriptions(LinearError)))
    end
    return LinearErrorResult(formulation, values, traversed.details)
end

function _linear_error(::MeasurementsPropagation, formulation::LinearError, problem,
        modal, modal_options)
    throw(ArgumentError("LinearError does not support modal composition; first-order mode tracking is missing"))
end

function _linear_error(::MeasurementsPropagation,
        formulation::LinearError{<:Engine.TransientAnalysis.TransientFormulation}, problem)
    throw(ArgumentError("LinearError does not support transient responses; first-order transient propagation is missing"))
end

function _linear_error(::MeasurementsPropagation,
        formulation::LinearError{<:LineCableModels.ModalAnalysisFormulation}, problem)
    throw(ArgumentError("LinearError does not support modal composition; first-order mode tracking is missing"))
end
