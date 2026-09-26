function compute(problem::ParametricProblem, formulation::LinearError)
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
        traversed = traverse(problem, formulation)
        values=map(traversed.values) do value
            id=get(Grammar.observation_gridpoint(value), :id, nothing)
            Engine.retain_gridpoint(value, id;
                fields = (
                    uncertainty = (
                        estimator = :first_order, representation = :dependency_preserving),
                    uncertainty_descriptions = _uncertainty_descriptions(LinearError)))
        end
        result = LinearErrorResult(formulation, values, traversed.details)
        progress &&
            @info "Linear error computation completed successfully" _group=:progress completed=length(result) total=length(result) elapsed_seconds=(time_ns()-started)*1e-9
        result
    end
end
