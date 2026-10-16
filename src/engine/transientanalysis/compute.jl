function validate(problem::TransientAnalysisProblem, selected::Formula)
    throw(ArgumentError("transient formula :$(formula_id(selected)) does not accept $(typeof(problem.line))"))
end

function computation_options(::Type{<:TransientFormulation}, record::ComputationOptions)::ComputationOptions
    options = record.data
    isempty(setdiff(keys(options), (:verbosity, :timing))) ||
        throw(ArgumentError("unknown transient action options: $(Tuple(keys(options)))"))
    timing = get(options, :timing, false)::Bool
    levels = verbosity(get(options, :verbosity, (default=0,)))
    return ComputationOptions(; verbosity=levels, timing)
end

"""
$(TYPEDSIGNATURES)

Compute a transient response using the selected formula. Execution options admit
`verbosity=(default=0,)` and `timing=false`. Both built-in solvers are placeholders.
"""
function compute(problem::TransientAnalysisProblem, formulation::TransientFormulation;
        options::Union{NamedTuple, ComputationOptions}=ComputationOptions())
    options = ComputationOptions(options)
    execution = computation_options(TransientFormulation, options)
    selected = formulation.formula
    validate(problem, selected)
    input = (line=problem.line, sending=problem.sending, receiving=problem.receiving,
        window=problem.window)
    coordinates = _transient_coordinates(problem.line)
    plan = (window=problem.window, coordinates=coordinates)
    T = float(promote_type(eltype(problem.line), typeof(problem.window.tmax),
        typeof(problem.window.step)))
    buffers = initialize_buffers(selected, T, input, plan, (;))
    workspace = TransientWorkspace(input, plan, buffers)
    values, timing = with_logger(VerbosityLogger(current_logger(), execution.data.verbosity)) do
        if execution.data.timing
            measured = Base.@timed respond!(selected, workspace, selected.parameters, selected.options)
            return measured.value, (wall_seconds=measured.time, bytes=measured.bytes,
                gc_seconds=measured.gctime, compile_seconds=measured.compile_time,
                recompile_seconds=measured.recompile_time)
        end
        return respond!(selected, workspace, selected.parameters, selected.options), nothing
    end
    t, voltages, currents = values
    retained = (formula=formula_id(selected), options=selected.options.data,
        window=problem.window, gridpoint=Commons.gridpoint_id(), coordinates=coordinates)
    execution.data.timing && (retained = merge(retained, (timing=timing,)))
    return TransientResponse(t, voltages, currents, coordinates, ComputationDetails(retained))
end

function computation_details(::Type{<:TransientFormulation}, result::TransientResponse)::ComputationDetails
    return details(result)
end

respond!(::Formula{:nlt}, workspace, parameters, options) = throw(ArgumentError("not implemented yet"))
respond!(::Formula{:td}, workspace, parameters, options) = throw(ArgumentError("not implemented yet"))
