"""
    resolve(context, objects...; options)

Return what `context` needs for the `objects` that apply to it. A geometric context resolves
one primitive or contextual layer against itself. A problem resolves its formulations.
`resolve(problem, formulations; options)` returns one [`Resolved`](@ref) per formulation, a
tuple for a tuple, and shares between them the parts that the formulations need equally.
[`traverse`](@ref) resolves each problem once, before its points. A backend that shares work
between the formulations of one problem adds a method. The default method gives each
formulation `Resolved(problem, (;))`.
"""
function resolve end

"""
$(TYPEDEF)

A problem and what one formulation needs from it, as [`resolve`](@ref) returns it. The parts
that the formulations of one problem need equally are shared between their `Resolved`
values, not copied. A backend computes one point with `compute(resolved, formulation;
options)`.

$(TYPEDFIELDS)
"""
struct Resolved{P <: AbstractProblemDefinition, I}
    "The problem."
    problem::P
    "What one formulation needs from the problem."
    input::I
end

resolve(problem::AbstractProblemDefinition, formulations; options) =
    map(Returns(Resolved(problem, (;))), formulations)

"""
    retain_gridpoint(result, gridpoint; fields = (;))

Return `result` with `gridpoint` as the identity of its point and `fields` merged into its
retained details. The owner of a result type adds a method. This default returns any other
value as it is.
"""
retain_gridpoint(result, gridpoint; fields = (;)) = result

"""
    timed(timing, call)

Return the value of `call()` and the fields that a point retains about the call.
`timed(Val(false), call)` returns `(call(), (;))`. `timed(Val(true), call)` returns the value
and `(timing = (; wall_seconds, bytes, gc_seconds, compile_seconds, recompile_seconds),)`,
measured with `Base.@timed`.
"""
timed(::Val{false}, call) = (call(), (;))
function timed(::Val{true}, call)
    statistics = Base.@timed call()
    return statistics.value, (timing = (wall_seconds = statistics.time,
        bytes = statistics.bytes, gc_seconds = statistics.gctime,
        compile_seconds = statistics.compile_time,
        recompile_seconds = statistics.recompile_time),)
end

# Without an owner method for its type, a formulation receives the options as they are.
computation_options(::Type{<:AbstractFormulation}, record::ComputationOptions) = record

"""
    traverse(problem, formulation; options)

Compute the Cartesian product of the points of `problem` and `formulation`, each a scalar or
a space ([`ScalarOrGrid`](@ref)). The traversal materializes and validates the problems and
the formulations before the first point. It computes problem by problem: it resolves each
problem once with [`resolve`](@ref), and each point is the backend's `compute(resolved,
formulation; options)`.

Before the first point, the owner of the formulations normalizes the options once with
`computation_options(F, options)`, where `F` is the element type of the formulations. A
formulation type without an owner method receives them as they are. The traversal writes the
work around the points once. The library emits its records only at the progress levels that
ask for them. When the options set `verbosity` or `on_result`, each point computes in a
`VerbosityLogger` at the `verbosity` levels, level 0 by default, which filters the records
of `on_result` too.
The traversal logs the start, the progress and the end of the run. With `timing = true`,
each value records `wall_seconds`, `bytes`, `gc_seconds`, `compile_seconds` and
`recompile_seconds`. Each value of a space receives its gridpoint. After each point,
`on_result(problem, index, value)` runs, with the position of the point in the space as
`index`.

A scalar problem and a scalar formulation are a run of one point, and `traverse` returns its
value. Otherwise it returns the values in the order of the space, the problem index fastest,
the retained details and the axes, `(; values, details, axes)`, from which a higher-order
formulation builds its result. A [`Combinatorial`](@ref) formulation sets the retention of
details.
"""
Base.@constprop :aggressive function traverse(
        problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    record = validate(ComputationOptions(options), traverse)
    # A problem's constructor validates it, and materializing a problem constructs it. A
    # scalar call keeps its one problem and formulation in tuples. A space keeps vectors,
    # which the steps of its points share by reference. Its formulations take the element type
    # of their source and are the record of the formulation axis.
    spaced = has_grid((problem, formulation))
    problems = map(materialize, points(problem))
    formulations = spaced ? collect(Grid(Combinatorial(formulation).inner)) :
                   map(source -> validate(materialize(source), Combinatorial),
        points(validate(formulation, Combinatorial)))
    isempty(problems) && throw(ArgumentError("a traversal needs at least one problem"))
    isempty(formulations) && throw(ArgumentError("a traversal needs at least one formulation"))
    # Constant options do not reach the steps defined in the run below. The owner of the
    # formulations normalizes them here, once, into a record in which the options that change
    # the type of a value are types.
    # Timing changes the type of the details, and its switch becomes a type here too.
    execution = computation_options(eltype(formulations), record)
    timing = get(record.data, :timing, false) ? Val(true) : Val(false)
    return traverse(problem, formulation, problems, formulations, execution, timing;
        options = record)
end

# The run on the decided parts. Options read at run time give a small union of record types,
# and each concrete record type compiles its own run.
function traverse(problem, formulation, problems, formulations, execution::ComputationOptions,
        timing::Val; options::ComputationOptions)
    levels = verbosity(get(options.data, :verbosity, (default = 0,)))
    on_result = get(options.data, :on_result, nothing)
    spaced = has_grid((problem, formulation))
    progress = get(levels, :progress, levels.default) > 0
    count = length(problems) * length(formulations)
    started = time_ns()
    # A point of a space takes the identity of the run and its indices. A scalar value keeps
    # the identity that its backend gave it, unless it gains fields to retain.
    identified = if spaced
        let source_id = gridpoint_id().source_id
            (value, fields, problem_index, formulation_index) -> retain_gridpoint(value,
                gridpoint_id(; source_id, problem_index, formulation_index); fields)
        end
    else
        (value, fields, _, _) -> isempty(fields) ? value :
                                 retain_gridpoint(value, gridpoint_id(); fields)
    end
    # A point computes its value with the backend, timed or not, then records its gridpoint
    # and calls `on_result`. The step receives the resolved input and the formulation of the
    # point. A space of formulations of different types dispatches once per point, and the
    # step compiles for each formulation type.
    point = function (problem_index, formulation_index, input, selected)
        result, fields = timed(timing, () -> compute(input, selected; options = execution))
        result = identified(result, fields, problem_index, formulation_index)
        index = problem_index + (formulation_index - 1) * length(problems)
        on_result === nothing || on_result(problems[problem_index], index, result)
        return result
    end
    # The library emits its records only at the progress levels that ask for them. Only a
    # logger scope filters the records of user code. When the options set `verbosity` or
    # `on_result`, each point runs in a `VerbosityLogger`, `on_result` included. The step
    # builds the logger from the levels, whose type is concrete, and not from a captured
    # current logger, whose type inference cannot know. This is decided once, before the
    # first point.
    logged = if haskey(options.data, :verbosity) || on_result !== nothing
        function (arguments...)
            logger = VerbosityLogger(Logging.current_logger(), levels)
            return Logging.with_logger(logger) do
                point(arguments...)
            end
        end
    else
        point
    end
    # Progress is decided once. With progress, the step reports it around each point.
    # Otherwise the step computes the point alone.
    step = if progress
        action = description(first(problems))
        completed, last_log, previous, average = Ref(0), Ref(started), Ref(started), Ref(0.0)
        function (arguments...)
            result = logged(arguments...)
            now = time_ns()
            completed[] += 1
            interval = (now - previous[]) * 1e-9
            average[] = completed[] == 1 ? interval : 0.2 * interval + 0.8 * average[]
            previous[] = now
            if now - last_log[] >= 5_000_000_000
                elapsed = (now - started) * 1e-9
                eta = (count - completed[]) * average[] / 3600
                @info "$action progress" _group = :progress completed = completed[] total = count elapsed_seconds = elapsed eta_hours = eta
                last_log[] = now
            end
            return result
        end
    else
        logged
    end
    close_run = function (result)
        if progress
            name = description(first(problems))
            elapsed = (time_ns() - started) * 1e-9
            @info "$name computation completed successfully" _group = :progress completed = count total = count elapsed_seconds = elapsed
        end
        return result
    end
    progress && @info "$(description(first(problems))) computation started" _group = :progress total = count
    first_resolved = resolve(first(problems), formulations; options = execution)
    value = step(1, 1, first(first_resolved), first(formulations))
    spaced || return close_run(value)
    values = Vector{typeof(value)}(undef, count)
    values[1] = value
    slots = reshape(values, length(problems), length(formulations))
    map!(index -> step(1, index, first_resolved[index], formulations[index]),
        view(slots, 1, 2:length(formulations)), 2:length(formulations))
    for problem_index in 2:length(problems)
        resolved = resolve(problems[problem_index], formulations; options = execution)
        map!(index -> step(problem_index, index, resolved[index], formulations[index]),
            view(slots, problem_index, :), eachindex(formulations))
    end
    details = if Combinatorial(formulation).options.data.retain_details
        first_record = computation_details(typeof(first(formulations)), first(values))
        records = Vector{typeof(first_record)}(undef, count)
        map!(index -> computation_details(
                typeof(formulations[cld(index, length(problems))]), values[index]),
            records, eachindex(values))
        ComputationDetails((points = records,))
    else
        ComputationDetails()
    end
    return close_run((; values, details, axes = (problems = Grid(problem), formulations)))
end

# A traversal reads its timing switch from the computation options.
function validate(options::ComputationOptions, ::typeof(traverse))
    timing = get(options.data, :timing, false)
    timing === true || timing === false || throw(ArgumentError("timing must be Bool"))
    return options
end

"""
    compute(problem, formulation; options)

Compute `problem` with `formulation`, each a scalar or a space, as the
[`ParametricProblem`](@ref) of `problem` and `options`.
"""
function compute(problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    return compute(ParametricProblem(problem, ComputationOptions(options)), formulation)
end

"""
    compute(problem::ParametricProblem, formulation)

Compute every point of the problem's space with `formulation`, a scalar or a space of
formulations or a [`Combinatorial`](@ref) formulation, through [`traverse`](@ref) with the
problem's options. A backend computes one point with `compute(resolved, formulation;
options)` on the [`Resolved`](@ref) value of [`resolve`](@ref). A scalar problem and a scalar
formulation return their value, and a space returns a [`ParametricResult`](@ref).
"""
function compute(problem::ParametricProblem,
        formulation::ScalarOrGrid{AbstractFormulation})
    has_grid((problem.space, formulation)) ||
        return traverse(problem.space, formulation; options = problem.options)
    combinatorial = Combinatorial(formulation)
    traversed = traverse(problem.space, combinatorial; options = problem.options)
    return ParametricResult(combinatorial, traversed.values, traversed.axes, traversed.details)
end
