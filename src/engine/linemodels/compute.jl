import ..Engine: completion_details

function computation_options(::Type{<:LineModelFormulation}, record::ComputationOptions)::ComputationOptions
    options = record.data
    isempty(setdiff(keys(options), (:verbosity, :timing))) ||
        throw(ArgumentError("unknown line-model action options: $(Tuple(keys(options)))"))
    timing = get(options, :timing, false)::Bool
    levels = verbosity(get(options, :verbosity, (default=0,)))
    return ComputationOptions(; verbosity=levels, timing)
end

"""
$(TYPEDSIGNATURES)

Fit a propagation segment with the selected formula. Execution options admit
`verbosity=(default=0,)` and `timing=false`.
"""
function compute(problem::LineModelProblem, formulation::LineModelFormulation;
        options::Union{NamedTuple, ComputationOptions}=ComputationOptions())
    execution = computation_options(LineModelFormulation, ComputationOptions(options))
    selected = formulation.formula
    segment = problem.segment
    input = (segment=segment, coordinates=copy(details(segment).data.phase_coordinates),
        length=convert(eltype(frequencies(segment)), line_length(segment)), frequencies=copy(frequencies(segment)))
    plan = (;)
    buffers = initialize_buffers(selected, eltype(problem), input, plan, (;))
    workspace = LineModelWorkspace(input, plan, buffers)
    values, timing = with_logger(VerbosityLogger(current_logger(), execution.data.verbosity)) do
        if execution.data.timing
            measured = Base.@timed fit(selected, workspace, selected.parameters, selected.options)
            return measured.value, (wall_seconds=measured.time, bytes=measured.bytes,
                gc_seconds=measured.gctime, compile_seconds=measured.compile_time,
                recompile_seconds=measured.recompile_time)
        end
        return fit(selected, workspace, selected.parameters, selected.options), nothing
    end
    model, errors, passive, dc = values
    retained = (formula=formula_id(selected), options=selected.options.data,
        gridpoint=get(details(segment).data, :gridpoint, nothing),
        inputs=get(details(segment).data, :inputs, nothing))
    execution.data.timing && (retained = merge(retained, (timing=timing,)))
    return LineModel(model, input.coordinates, input.length, input.frequencies,
        errors, passive, dc, completion_details(retained))
end

function computation_details(::Type{<:LineModelFormulation}, result::LineModel)::ComputationDetails
    return details(result)
end

fit(::Formula{:ulm}, args...) = throw(ArgumentError("load UniversalLineModel to fit :ulm"))
fit(::Formula{:fdcm}, args...) = throw(ArgumentError("load UniversalLineModel to fit :fdcm"))
fit(::Formula{:marti}, args...) = throw(ArgumentError("not implemented yet"))
