"""
    resolve(context, objects...; options)

Return what `context` needs for the `objects` that apply to it. A geometric context resolves
one primitive or contextual layer against itself.
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
    retain_gridpoint(result; fields)

Return `result` with `fields` merged into its retained details. The gridpoint identity of
`result` is unchanged. The owner of a result type adds a method. This default returns any
other value as it is.
"""
retain_gridpoint(result; fields = (;)) = result

# Each formulation space contains the formulations of one backend, whose owner normalizes
# their options. This method refuses a space that mixes backends and a formulation type
# without an owner method.
function computation_options(::Type{F}, ::ComputationOptions) where {F <: AbstractFormulation}
    throw(ArgumentError(
        "no computation options are defined for formulations of type $F. A formulation " *
        "space contains the formulations of one backend: compute each backend separately. " *
        "A formulation type of your own defines computation_options(::Type{<:T}, options). " *
        "See `?computation_options`.",
    ))
end

"""
    traverse(problem, formulation; options)

Compute the Cartesian product of the points of `problem` and `formulation`, each a scalar or
a space ([`ScalarOrGrid`](@ref)). The traversal materializes and validates the problems and
the formulations before the first point. It computes problem by problem: it resolves each
problem once with [`resolve`](@ref), and each point is the backend's `compute(resolved,
formulation; options)`.

Before the first point, the owner of the formulations normalizes the options once with
`computation_options(F, options)`, where `F` is the element type of the formulations. A
formulation space contains the formulations of one backend, and a space that mixes backends
raises an `ArgumentError`. The library emits its records only at the levels that enable
them. When the options set `verbosity` or `on_result`, the run computes in one
`VerbosityLogger` at the `verbosity` levels, level 0 by default, which filters the records
of `on_result` too.

The traversal logs the start, the progress and the end of the run in the group
`:progress`. At level 1 or more in the group `:timing`, each point writes one `Info` record
in that group with its `gridpoint` and the `wall_seconds`, `bytes`, `gc_seconds`,
`compile_seconds` and `recompile_seconds` of its computation, measured with `Base.@timed`.
A result does not record timing. Each value of a space receives its gridpoint. After each
point, `on_result(problem, index, value)` runs, with the position of the point in the space
as `index`.

A scalar problem and a scalar formulation are a run of one point, and `traverse` returns its
value. Otherwise it returns the values in the order of the space, the problem index fastest,
the retained details and the gridpoints, the problem and formulation indices of the points,
`(; values, details, gridpoints)`, from which a higher-order formulation builds its result.
A [`Combinatorial`](@ref) formulation sets the retention of details.
"""
function traverse(
        problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    record = ComputationOptions(options)
    # A problem's constructor validates it, and materializing a problem constructs it. A
    # scalar call keeps its one problem and formulation in tuples. A space keeps vectors.
    # Its formulations take the element type of their source.
    spaced = has_grid((problem, formulation))
    problems = map(materialize, points(problem))
    formulations = spaced ? collect(Grid(Combinatorial(formulation).inner)) :
                   map(source -> validate(materialize(source), Combinatorial),
        points(validate(formulation, Combinatorial)))
    isempty(problems) && throw(ArgumentError("a traversal needs at least one problem"))
    isempty(formulations) && throw(ArgumentError("a traversal needs at least one formulation"))
    # The owner of the formulations normalizes the options here, once, into a record in
    # which the options that change the type of a value are types.
    execution = computation_options(eltype(formulations), record)
    levels = verbosity(get(record.data, :verbosity, (default = 0,)))
    on_result = get(record.data, :on_result, nothing)
    # The run's controls are data: the levels decide progress and timing once, and
    # `on_result` is the user's callable. The points of a run share its source identity.
    run = (; problems, formulations, execution, on_result,
        count = length(problems) * length(formulations),
        spaced, source_id = gridpoint_id().source_id,
        progress = get(levels, :progress, levels.default) > 0,
        timing = get(levels, :timing, levels.default) > 0, started = time_ns())
    # The types of a space are known only from its values. Its run is called in the latest
    # world, which inference does not enter: the run compiles from the values it receives,
    # and the compiler does not walk the computation with partly known types. A scalar run
    # infers.
    scoped = haskey(record.data, :verbosity) || on_result !== nothing
    spaced || return traverse(problem, formulation, run, scoped, levels)
    return Base.invokelatest(traverse, problem, formulation, run, scoped, levels)
end

# Only a logger scope filters the records of user code. When the options set `verbosity` or
# `on_result`, the run computes in one `VerbosityLogger`.
function traverse(problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation}, run::NamedTuple, scoped::Bool,
        levels::NamedTuple)
    scoped || return traverse(problem, formulation, run)
    return Logging.with_logger(() -> traverse(problem, formulation, run),
        VerbosityLogger(Logging.current_logger(), levels))
end

# The run computes problem by problem and resolves each problem once. A scalar run returns
# the value of its one point, and a space allocates its values from the type of the first.
function traverse(problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation}, run::NamedTuple)
    (; problems, formulations, execution) = run
    run.progress && @info "$(description(first(problems))) computation started" _group = :progress total = run.count
    resolved = resolve(first(problems), formulations; options = execution)
    value = traverse(run, 1, 1, first(resolved), first(formulations))
    result = if has_grid((problem, formulation))
        values = Vector{typeof(value)}(undef, run.count)
        values[1] = value
        traverse(values, run, resolved)
        details = if Combinatorial(formulation).options.data.retain_details
            first_record = computation_details(typeof(first(formulations)), first(values))
            records = Vector{typeof(first_record)}(undef, run.count)
            for index in eachindex(values)
                records[index] = computation_details(
                    typeof(formulations[cld(index, length(problems))]), values[index])
            end
            ComputationDetails((points = records,))
        else
            ComputationDetails()
        end
        (; values, details,
            gridpoints = CartesianIndices((length(problems), length(formulations))))
    else
        value
    end
    run.progress && @info "$(description(first(problems))) computation completed successfully" _group = :progress completed = run.count total = run.count elapsed_seconds = (time_ns() - run.started) * 1e-9
    return result
end

# The points after the first, into `values` at their positions in the space, the problem
# index fastest. Each later problem is resolved once.
function traverse(values::Vector, run::NamedTuple, resolved::AbstractVector)
    (; problems, formulations, execution) = run
    for formulation_index in 2:length(formulations)
        values[1 + (formulation_index - 1) * length(problems)] = traverse(run, 1,
            formulation_index, resolved[formulation_index], formulations[formulation_index])
    end
    for problem_index in 2:length(problems)
        inputs = resolve(problems[problem_index], formulations; options = execution)
        for formulation_index in eachindex(formulations)
            values[problem_index + (formulation_index - 1) * length(problems)] =
                traverse(run, problem_index, formulation_index, inputs[formulation_index],
                    formulations[formulation_index])
        end
    end
    return values
end

# The point computes its value with the backend and the formulation of the point. A point
# of a space takes its gridpoint, and then `on_result` receives the value. This method
# dispatches on the type of the formulation, once per point of a space of formulations of
# different types. A point that ends in a later five-second interval of the run than the
# one in which it started logs the progress from its position and the elapsed time.
function traverse(run::NamedTuple, problem_index::Int, formulation_index::Int, input,
        selected::AbstractFormulation)
    started = time_ns()
    gridpoint = gridpoint_id(; run.source_id, problem_index, formulation_index)
    value = if run.timing
        statistics = Base.@timed compute(input, selected; options = run.execution)
        @info "Point computed" _group = :timing gridpoint wall_seconds = statistics.time bytes = statistics.bytes gc_seconds = statistics.gctime compile_seconds = statistics.compile_time recompile_seconds = statistics.recompile_time
        statistics.value
    else
        compute(input, selected; options = run.execution)
    end
    value = run.spaced ? retain_gridpoint(value, gridpoint) : value
    index = problem_index + (formulation_index - 1) * length(run.problems)
    run.on_result === nothing || run.on_result(run.problems[problem_index], index, value)
    run.progress || return value
    now = time_ns()
    if fld(now - run.started, 5_000_000_000) > fld(started - run.started, 5_000_000_000)
        completed = (problem_index - 1) * length(run.formulations) + formulation_index
        elapsed = (now - run.started) * 1e-9
        @info "$(description(input.problem)) progress" _group = :progress completed total = run.count elapsed_seconds = elapsed eta_hours = (run.count - completed) * elapsed / completed / 3600
    end
    return value
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
    return ParametricResult(traversed.values, traversed.gridpoints, traversed.details)
end
