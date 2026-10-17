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
`verbosity=(default=0,)` and `timing=false`.
"""
function compute(problem::TransientAnalysisProblem, formulation::TransientFormulation;
        options::Union{NamedTuple, ComputationOptions}=ComputationOptions())
    options = ComputationOptions(options)
    execution = computation_options(TransientFormulation, options)
    selected = formulation.formula
    validate(problem, selected)
    values, timing = with_logger(VerbosityLogger(current_logger(), execution.data.verbosity)) do
        if execution.data.timing
            measured = Base.@timed compute(selected, problem, execution)
            return measured.value, (wall_seconds=measured.time, bytes=measured.bytes,
                gc_seconds=measured.gctime, compile_seconds=measured.compile_time,
                recompile_seconds=measured.recompile_time)
        end
        return compute(selected, problem, execution), nothing
    end
    t, voltages, currents, coordinates, retained = values
    retained = merge((formula=formula_id(selected), options=selected.options.data,
        window=problem.window, gridpoint=Commons.gridpoint_id(), coordinates=coordinates), retained)
    execution.data.timing && (retained = merge(retained, (timing=timing,)))
    return TransientResponse(t, voltages, currents, coordinates, ComputationDetails(retained))
end

function computation_details(::Type{<:TransientFormulation}, result::TransientResponse)::ComputationDetails
    return details(result)
end

compute(::Formula{:td}, problem::TransientAnalysisProblem{<:LineModel}, execution) =
    throw(ArgumentError("load UniversalLineModel to run :td on a wideband line model"))
function compute(selected::Formula{:td},
        problem::TransientAnalysisProblem{<:LineModel{<:ModalModel}}, execution)
    workspace = TransientWorkspace(problem, selected)
    t, voltages, currents = respond!(selected, workspace, selected.parameters, selected.options)
    matched = workspace.buffers.matched
    termination = (impedance=:fitted_characteristic,
        poles=sum(section -> length(section.poles), problem.line.model.impedance))
    retained = (td=(interp=:single, enforce_sign=false,
        initial_solve=selected.options.data.initial_solve, steps=length(t),
        terminations=(sending=matched[1] ? termination : nothing,
            receiving=matched[2] ? termination : nothing)), line=details(problem.line))
    return t, voltages, currents, problem.line.coordinates, retained
end
