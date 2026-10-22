"""
    resolve(context, objects...; options)

Return what `context` needs for the `objects` that apply to it. A geometric context resolves
one primitive or contextual layer against itself.
`resolve(problem, formulations; options)` returns one [`Resolved`](@ref) per formulation, a
tuple for a tuple, and shares between them the parts that the formulations need equally.
[`traverse`](@ref) resolves each problem once, before its points. A backend that shares work
between the formulations of one problem adds a method. The default method gives each
formulation `Resolved(problem, (;))`. A Monte Carlo study builds a `Resolved` value for
each point before the run, with the uncertain point beside the nominal problem, and
`resolve` returns that value for each formulation.
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
resolve(resolved::Resolved, formulations; options) = map(Returns(resolved), formulations)
# The description of a `Resolved` value is the description of its problem.
description(resolved::Resolved; compact::Bool = false) = description(resolved.problem; compact)

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
    on_error(study, exception)

Return `exception` when `study` records it as the failure of a point, and stop `study`
otherwise. A study calls it inside the `catch` block of the point that threw, and records
that exception when it returns. A [`Combinatorial`](@ref) study records the failed point
beside its values and continues, and `LinearError` does so through the `Combinatorial`
study that it traverses. A `MonteCarlo` study draws the trial again, and the
`Combinatorial` study of its points records a point whose nominal design or whose draws are
infeasible. The fallback stops the study with `rethrow()`. The exception keeps its
backtrace, and a bug stops a study with its stack trace.

A `Combinatorial` study records a `DomainError`, an input outside the domain of a formula,
and an `OverlapError`, cable exteriors that overlap. The owner of an exception type declares
with a method what a study records of it. A user records an exception type of their own with
one method:

```julia
LineCableModels.on_error(::Combinatorial, exception::MyError) = exception
```
"""
on_error(study, exception) = rethrow()
on_error(::Combinatorial, exception::DomainError) = exception

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
computed point, `on_result(problem, index, value)` runs, with the position of the point in
the space as `index`.

A scalar problem and a scalar formulation are a run of one point, and `traverse` returns its
value. A space records a failed point and continues. The build or the resolve of its
problem, or its computation, fails the point when it throws an exception that
[`on_error`](@ref) returns for the [`Combinatorial`](@ref) study, and a failed build or
resolve fails every point of its problem. Any other exception stops the run with its
backtrace. A space returns the values of its computed points in the order of the space, the
problem index fastest, the failures in the same order, each the pair of the gridpoint
identity of a failed point and its exception, the retained details of the values and the
gridpoints, the problem and formulation indices of every point, `(; values, failures,
details, gridpoints)`, from which a higher-order formulation builds its result. The values
take the result type of the earliest point that computes. A space with no computed point
throws an error that gives the count of its points and shows the first recorded exception.
A `Combinatorial` formulation sets the retention of details.
"""
function traverse(
        problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    record = ComputationOptions(options)
    # A scalar call keeps its one problem and formulation in tuples, and a space keeps
    # vectors, its formulations with the element type of their source.
    spaced = has_grid((problem, formulation))
    sources = points(problem)
    formulations = spaced ? collect(Grid(Combinatorial(formulation).inner)) :
                   map(source -> validate(materialize(source), Combinatorial),
        points(validate(formulation, Combinatorial)))
    isempty(sources) && throw(ArgumentError("a traversal needs at least one problem"))
    isempty(formulations) && throw(ArgumentError("a traversal needs at least one formulation"))
    # The owner of the formulations normalizes the options here, once, into a record in
    # which the options that change the type of a value are types.
    execution = computation_options(eltype(formulations), record)
    levels = verbosity(get(record.data, :verbosity, (default = 0,)))
    on_result = get(record.data, :on_result, nothing)
    # The run's controls are data: the levels decide progress and timing once, and
    # `on_result` is the user's callable. The points of a run share its source identity,
    # and its gridpoints are their problem and formulation indices.
    gridpoints = CartesianIndices((length(sources), length(formulations)))
    source_id = gridpoint_id().source_id
    progress = get(levels, :progress, levels.default) > 0
    timing = get(levels, :timing, levels.default) > 0
    started = time_ns()
    scoped = haskey(record.data, :verbosity) || on_result !== nothing
    # A problem's constructor validates it, and materializing a problem constructs it. A
    # scalar run infers.
    spaced || return traverse(problem, formulation,
        (; problems = map(materialize, sources), formulations, execution, on_result, gridpoints,
            spaced, source_id, progress, timing, started), scoped, levels)
    # A space builds every problem before its first point. When `on_error` returns the
    # exception of a build for the study, every point of its problem fails and the exception
    # takes the place of the problem. The problems that built keep the type that `map` gives
    # them, with their indices in the space.
    study = Combinatorial(formulation)
    failures = fieldtype(ParametricResult, :failures)()
    attempts = has_grid(problem) ? map(enumerate(sources)) do (problem_index, source)
        try
            materialize(source)
        catch exception
            on_error(study, exception)
            append!(failures, (gridpoint_id(; source_id, problem_index, formulation_index) =>
                exception for formulation_index in eachindex(formulations)))
            exception
        end
    end : map(materialize, sources)
    problem_indices = isempty(failures) ? eachindex(attempts) :
                      setdiff!(collect(eachindex(attempts)),
        (gridpoint.problem_index for (gridpoint, _) in failures))
    problems = isempty(failures) ? attempts :
               map(Base.Fix1(getindex, attempts), problem_indices)
    # The types of a space are known only from its values. Its run is called in the latest
    # world, which inference does not enter: the run compiles from the values it receives,
    # and the compiler does not walk the computation with partly known types.
    return Base.invokelatest(traverse, problem, formulation,
        (; problems, problem_indices, study, failures, formulations, execution, on_result,
            gridpoints, spaced, source_id, progress, timing, started), scoped, levels)
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

# A scalar run returns the value of its one point, and its exception propagates. A space
# searches for its first computed point.
function traverse(problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::ScalarOrGrid{AbstractFormulation}, run::NamedTuple)
    has_grid((problem, formulation)) && return traverse(run)
    (; problems, formulations, execution) = run
    run.progress && @info "$(description(first(problems))) computation started" _group = :progress total = 1
    resolved = resolve(first(problems), formulations; options = execution)
    value = traverse(run, 1, 1, first(resolved), first(formulations))
    run.progress && @info "$(description(first(problems))) computation completed successfully" _group = :progress completed = 1 total = 1 elapsed_seconds = (time_ns() - run.started) * 1e-9
    return value
end

# The points of the problems that built, problem by problem, each problem resolved once,
# until one computes: its result type is the element type of the values. A resolve or a
# point whose exception `on_error` returns for the study fails, and a failed resolve fails
# every point of its problem.
function traverse(run::NamedTuple)
    (; problems, problem_indices, formulations, failures) = run
    run.progress && !isempty(problems) && @info "$(description(first(problems))) computation started" _group = :progress total = length(run.gridpoints)
    for (position, problem) in enumerate(problems)
        problem_index = problem_indices[position]
        resolved = try
            resolve(problem, formulations; options = run.execution)
        catch exception
            on_error(run.study, exception)
            append!(failures, (gridpoint_id(; run.source_id, problem_index, formulation_index) =>
                exception for formulation_index in eachindex(formulations)))
            continue
        end
        for formulation_index in eachindex(formulations)
            value = try
                traverse(run, problem_index, formulation_index, resolved[formulation_index],
                    formulations[formulation_index])
            catch exception
                on_error(run.study, exception)
                push!(failures,
                    gridpoint_id(; run.source_id, problem_index, formulation_index) => exception)
                continue
            end
            values = Vector{typeof(value)}(undef, length(run.gridpoints))
            values[LinearIndices(run.gridpoints)[problem_index, formulation_index]] = value
            return traverse(values, run, position, formulation_index, resolved)
        end
    end
    gridpoint, recorded = first(failures)
    throw(ErrorException("all $(length(run.gridpoints)) points of the space failed. The " *
        "first recorded, at problem $(gridpoint.problem_index) and formulation " *
        "$(gridpoint.formulation_index): " * sprint(showerror, recorded)))
end

# The points after the point of the problem at `position` and the formulation at `from`,
# into `values` at their positions in the space, each later problem resolved once. The
# failures take the order of the space, the failed points leave the values, and the
# retained details are those of the values.
function traverse(values::Vector, run::NamedTuple, position::Int, from::Int, resolved)
    (; problems, problem_indices, formulations, failures, gridpoints) = run
    for current in position:length(problems)
        problem_index = problem_indices[current]
        if current > position
            resolved = try
                resolve(problems[current], formulations; options = run.execution)
            catch exception
                on_error(run.study, exception)
                append!(failures, (gridpoint_id(; run.source_id, problem_index,
                    formulation_index) => exception for formulation_index in eachindex(formulations)))
                continue
            end
        end
        for formulation_index in (current == position ? from + 1 : 1):length(formulations)
            try
                values[LinearIndices(gridpoints)[problem_index, formulation_index]] =
                    traverse(run, problem_index, formulation_index,
                        resolved[formulation_index], formulations[formulation_index])
            catch exception
                on_error(run.study, exception)
                push!(failures,
                    gridpoint_id(; run.source_id, problem_index, formulation_index) => exception)
            end
        end
    end
    computed = eachindex(values)
    if !isempty(failures)
        positions = Vector{Int}(undef, length(failures))
        for (index, (gridpoint, _)) in enumerate(failures)
            positions[index] =
                LinearIndices(gridpoints)[gridpoint.problem_index, gridpoint.formulation_index]
        end
        order = sortperm(positions)
        permute!(failures, order)
        failed = positions[order]
        computed = deleteat!(collect(computed), failed)
        deleteat!(values, failed)
    end
    details = if run.study.options.data.retain_details
        first_record = computation_details(
            typeof(formulations[cld(first(computed), size(gridpoints, 1))]), first(values))
        records = Vector{typeof(first_record)}(undef, length(values))
        for (index, point) in enumerate(computed)
            records[index] = computation_details(
                typeof(formulations[cld(point, size(gridpoints, 1))]), values[index])
        end
        ComputationDetails((points = records,))
    else
        ComputationDetails()
    end
    run.progress && @info "$(description(first(problems))) computation completed successfully" _group = :progress completed = length(values) total = length(gridpoints) elapsed_seconds = (time_ns() - run.started) * 1e-9
    return (; values, failures, details, gridpoints)
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
    index = LinearIndices(run.gridpoints)[problem_index, formulation_index]
    run.on_result === nothing || run.on_result(input.problem, index, value)
    run.progress || return value
    now = time_ns()
    if fld(now - run.started, 5_000_000_000) > fld(started - run.started, 5_000_000_000)
        completed = (problem_index - 1) * length(run.formulations) + formulation_index
        elapsed = (now - run.started) * 1e-9
        @info "$(description(input.problem)) progress" _group = :progress completed total = length(run.gridpoints) elapsed_seconds = elapsed eta_hours = (length(run.gridpoints) - completed) * elapsed / completed / 3600
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
formulation return their value, and a space returns a [`ParametricResult`](@ref) with the
failed points of the study.
"""
function compute(problem::ParametricProblem,
        formulation::ScalarOrGrid{AbstractFormulation})
    has_grid((problem.space, formulation)) ||
        return traverse(problem.space, formulation; options = problem.options)
    combinatorial = Combinatorial(formulation)
    traversed = traverse(problem.space, combinatorial; options = problem.options)
    return ParametricResult(traversed.values, traversed.failures, traversed.gridpoints,
        traversed.details)
end
