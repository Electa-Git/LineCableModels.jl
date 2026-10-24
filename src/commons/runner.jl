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
    resolve(study, problem, formulations; options)

Return the input of each point of `problem`, a problem of the space of `study`, one per
formulation: what the backend's `resolve(problem, formulations; options)` returns. A study
adds a method when its points compute something other than the problem. A `MonteCarlo`
point computes its trials, and the backend resolves the drawn problem of each trial.
"""
resolve(::AbstractTraversal, problem, formulations; options) =
    resolve(problem, formulations; options)
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
beside its values and continues, and a `LinearError` study records what `Combinatorial`
records. A `MonteCarlo` study draws a trial again when its draw is infeasible, and records a
point whose nominal design or whose draws are infeasible. The fallback stops the study with
`rethrow()`. The exception keeps its backtrace, and a bug stops a study with its stack trace.

A `Combinatorial` study records a `DomainError`, an input outside the domain of a formula,
and an `OverlapError`, cable exteriors that overlap. The owner of an exception type declares
with a method what a study records of it. A user records an exception type of their own with
one method for each study:

```julia
LineCableModels.on_error(::Combinatorial, exception::MyError) = exception
```
"""
on_error(study, exception) = rethrow()
on_error(::Combinatorial, exception::DomainError) = exception

"""
    traverse(problem::ParametricProblem, formulation; options)
    traverse(problem::ParametricProblem, study; options)

Compute the points of a run, whose `options` are the computation options of the run, and
the problem's options those of the backend. The [`ParametricProblem`](@ref) of a core
problem and a core formulation are a run of one point: `traverse` returns its value, and its
exception propagates. A `ParametricProblem` and a study, an `AbstractTraversal`, are a run
of the points of the problem's space with the formulations of the study.

Before the first point, each owner normalizes its options once: the backend, the owner of
the formulations, normalizes the options of the points with `computation_options(F,
options)`, where `F` is the element type of the formulations, and a study's run normalizes
its own with `computation_options(T, options)` for its type `T`. A formulation space contains
the formulations of one backend, and a space that mixes backends raises an `ArgumentError`.
Without `verbosity` or `on_result` among the options, the run does not filter log records,
and the caller's logger decides. With either, the run computes in one `VerbosityLogger`,
which filters the records, those of `on_result` too, by the `verbosity` levels, level 0 by
default. The runner's own progress and timing records appear only when `verbosity` asks for
them.

The traversal logs the start, the progress and the end of the run in the group
`:progress`. At level 1 or more in the group `:timing`, each point writes one `Info` record
in that group with its `gridpoint` and the `wall_seconds`, `bytes`, `gc_seconds`,
`compile_seconds` and `recompile_seconds` of its computation, measured with `Base.@timed`.
A result does not record timing. After each computed point, `on_result(problem, index,
value)` runs, with the position of the point in the space as `index`.

A study builds every problem of its space before the first point and computes problem by
problem. It builds each point with `materialize(study, source)` and resolves each problem
once with `resolve(study, problem, formulations; options)`, both of which the type of the
study selects. The value of a point is the backend's `compute(resolved, formulation;
options)`, or a `traverse` method on the type of a study whose points compute otherwise. A
`MonteCarlo` study takes the nominal problem of a point as `Resolved(problem, (;))`, and its
points compute their trials. Each value receives its gridpoint. A study records a failed
point and continues. The build or the resolve of its problem, or its computation, fails the
point when it throws an exception that [`on_error`](@ref) returns for the study, and a
failed build or resolve fails every point of its problem. Any other exception stops the run
with its backtrace. A study returns the values of its computed points in the order of the
space, the problem index fastest, the failures in the same order, each the pair of the
gridpoint identity of a failed point and its exception, the retained details of the values
and the gridpoints, the problem and formulation indices of every point, `(; values,
failures, details, gridpoints)`, from which the study builds its result. The values take the
result type of the earliest point that computes. A study with no computed point throws an
error that gives the count of its points and shows the first recorded exception. The option
`retain_details` of the run sets the retention of details.
"""
function traverse(problem::ParametricProblem{<:AbstractProblemDefinition},
        formulation::AbstractFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    given = ComputationOptions(options)
    # A run of one point keeps its problem and formulation in tuples, and it infers. A
    # problem's constructor validates it. The run normalizes its own options once, as the
    # run of a study does. The backend normalizes those of the point.
    formulations = (formulation,)
    record = computation_options(AbstractTraversal, given)
    execution = computation_options(typeof(formulation), problem.options)
    return traverse(problem.space, formulation,
        (; problems = (problem.space,), formulation, formulations, execution, record,
            gridpoints = CartesianIndices((1, 1)), spaced = false,
            source_id = gridpoint_id().source_id, started = time_ns()), given)
end

function traverse(problem::ParametricProblem, study::AbstractTraversal;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    given = ComputationOptions(options)
    (; space) = problem
    # The points of the problem's space and the formulations of the study's points,
    # materialized. A scalar problem or formulation is a source of one point.
    sources = points(space)
    formulations = map(materialize, points(study))
    isempty(formulations) && throw(ArgumentError("a traversal needs at least one formulation"))
    # Each owner normalizes its options here, once. The study normalizes the options of its
    # run. The owner of the formulations normalizes those of the points into a record in
    # which the options that change the type of a value are types.
    record = computation_options(typeof(study), given)
    execution = computation_options(eltype(formulations), problem.options)
    gridpoints = CartesianIndices((length(sources), length(formulations)))
    source_id = gridpoint_id().source_id
    # A study builds every problem before its first point. When `on_error` returns the
    # exception of a build for the study, every point of its problem fails and the exception
    # takes the place of the problem. The problems that built keep the type that `map` gives
    # them, with their indices in the space.
    failures = fieldtype(ParametricResult, :failures)()
    attempts = map(enumerate(sources)) do (problem_index, source)
        try
            materialize(study, source)
        catch exception
            on_error(study, exception)
            append!(failures, (gridpoint_id(; source_id, problem_index, formulation_index) =>
                exception for formulation_index in eachindex(formulations)))
            exception
        end
    end
    problem_indices = isempty(failures) ? eachindex(attempts) :
                      setdiff!(collect(eachindex(attempts)),
        (gridpoint.problem_index for (gridpoint, _) in failures))
    problems = isempty(failures) ? attempts :
               map(Base.Fix1(getindex, attempts), problem_indices)
    # The types of a space are known only from its values. Its run is called in the latest
    # world, which inference does not enter: the run compiles from the values it receives,
    # and the compiler does not walk the computation with partly known types.
    # The formulation of the run is the study, the formulation of the parametric problem.
    return Base.invokelatest(traverse, space, study,
        (; problems, problem_indices, formulation = study, failures, formulations, execution,
            record, gridpoints, spaced = true, source_id, started = time_ns()), given)
end

# Only a logger scope filters the records of user code. When the options given to the run
# set `verbosity` or `on_result`, the run computes in one `VerbosityLogger`.
function traverse(problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::Union{AbstractFormulation, AbstractTraversal}, run::NamedTuple,
        given::ComputationOptions)
    haskey(given.data, :verbosity) || get(given.data, :on_result, nothing) !== nothing ||
        return traverse(problem, formulation, run)
    return Logging.with_logger(() -> traverse(problem, formulation, run),
        VerbosityLogger(Logging.current_logger(),
            verbosity(get(given.data, :verbosity, (default = 0,)))))
end

# A study searches for its first computed point. A run of one point returns its value, and
# its exception propagates. The run's controls are its record: the levels decide progress
# and timing, and `on_result` is the user's callable. A study whose points read the space
# adds a method.
traverse(::ScalarOrGrid{AbstractProblemDefinition}, ::AbstractTraversal, run::NamedTuple) =
    traverse(run)
function traverse(::AbstractProblemDefinition, formulation::AbstractFormulation,
        run::NamedTuple)
    (; problems, formulations, execution, record) = run
    progress = verbosity(record, :progress) > 0
    progress && @info "$(description(first(problems))) computation started" _group = :progress total = 1
    resolved = resolve(first(problems), formulations; options = execution)
    value = traverse(run, 1, 1, first(resolved), first(formulations))
    progress && @info "$(description(first(problems))) computation completed successfully" _group = :progress completed = 1 total = 1 elapsed_seconds = (time_ns() - run.started) * 1e-9
    return value
end

# The points of the problems that built, problem by problem, each problem resolved once,
# until one computes: its result type is the element type of the values. A resolve or a
# point whose exception `on_error` returns for the study fails, and a failed resolve fails
# every point of its problem.
function traverse(run::NamedTuple)
    (; problems, problem_indices, formulation, formulations, failures) = run
    verbosity(run.record, :progress) > 0 && !isempty(problems) && @info "$(description(first(problems))) computation started" _group = :progress total = length(run.gridpoints)
    for (position, problem) in enumerate(problems)
        problem_index = problem_indices[position]
        resolved = try
            resolve(formulation, problem, formulations; options = run.execution)
        catch exception
            on_error(formulation, exception)
            append!(failures, (gridpoint_id(; run.source_id, problem_index, formulation_index) =>
                exception for formulation_index in eachindex(formulations)))
            continue
        end
        for formulation_index in eachindex(formulations)
            value = try
                traverse(run, problem_index, formulation_index, resolved[formulation_index],
                    formulations[formulation_index])
            catch exception
                on_error(formulation, exception)
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
    (; problems, problem_indices, formulation, formulations, failures, gridpoints) = run
    for current in position:length(problems)
        problem_index = problem_indices[current]
        if current > position
            resolved = try
                resolve(formulation, problems[current], formulations; options = run.execution)
            catch exception
                on_error(formulation, exception)
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
                on_error(formulation, exception)
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
    details = if run.record.data.retain_details
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
    verbosity(run.record, :progress) > 0 && @info "$(description(first(problems))) computation completed successfully" _group = :progress completed = length(values) total = length(gridpoints) elapsed_seconds = (time_ns() - run.started) * 1e-9
    return (; values, failures, details, gridpoints)
end

# The formulation of the run, the core formulation of a run of one point or the study,
# decides the value of a point. A point of a study then takes its gridpoint, and `on_result`
# receives the value. This method dispatches on the type of the
# formulation of the point, once per point of a space of formulations of different types. A
# point that ends in a later five-second interval of the run than the one in which it
# started logs the progress from its position and the elapsed time.
function traverse(run::NamedTuple, problem_index::Int, formulation_index::Int, input,
        selected::AbstractFormulation)
    started = time_ns()
    gridpoint = gridpoint_id(; run.source_id, problem_index, formulation_index)
    value = if verbosity(run.record, :timing) > 0
        statistics = Base.@timed traverse(run.formulation, run, gridpoint, input, selected)
        @info "Point computed" _group = :timing gridpoint wall_seconds = statistics.time bytes = statistics.bytes gc_seconds = statistics.gctime compile_seconds = statistics.compile_time recompile_seconds = statistics.recompile_time
        statistics.value
    else
        traverse(run.formulation, run, gridpoint, input, selected)
    end
    value = run.spaced ? retain_gridpoint(value, gridpoint) : value
    index = LinearIndices(run.gridpoints)[problem_index, formulation_index]
    on_result = get(run.record.data, :on_result, nothing)
    on_result === nothing || on_result(input.problem, index, value)
    verbosity(run.record, :progress) > 0 || return value
    now = time_ns()
    if fld(now - run.started, 5_000_000_000) > fld(started - run.started, 5_000_000_000)
        completed = (problem_index - 1) * length(run.formulations) + formulation_index
        elapsed = (now - run.started) * 1e-9
        @info "$(description(input.problem)) progress" _group = :progress completed total = length(run.gridpoints) elapsed_seconds = elapsed eta_hours = (length(run.gridpoints) - completed) * elapsed / completed / 3600
    end
    return value
end

# The backend computes the value of a point from the resolved input with the formulation of
# the point. A study whose points compute otherwise adds a method.
function traverse(::Union{AbstractFormulation, AbstractTraversal}, run::NamedTuple,
        ::NamedTuple, input, selected::AbstractFormulation)
    return compute(input, selected; options = run.execution)
end

"""
    compute(problem, formulation; options)

Compute `problem`, a core problem or a space of them, with `formulation`, a core formulation,
a space of them or a study, as the [`ParametricProblem`](@ref) of `problem`. The `options` go
to the owner of each key: the computation options of the run, `verbosity`, `on_result` and
`retain_details`, the keys that the run's `computation_options` accepts, to the call, and the
others to the problem, the computation options of the backend. A core problem and a core
formulation are a run of one point, which returns its value. A space computes the study, and
a core formulation its default study, `Combinatorial`.
"""
function compute(problem::ScalarOrGrid{AbstractProblemDefinition},
        formulation::Union{ScalarOrGrid{AbstractFormulation}, AbstractTraversal};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    record = ComputationOptions(options).data
    # The keys of the run are those of its normalized record.
    accepted = keys(computation_options(AbstractTraversal, ComputationOptions()).data)
    run = filter(in(accepted), keys(record))
    backend = filter(!in(accepted), keys(record))
    return compute(ParametricProblem(problem, ComputationOptions(record[backend])), formulation;
        options = ComputationOptions(record[run]))
end

"""
    compute(problem::ParametricProblem, formulation; options)

Compute every point of the problem's space with `formulation`, a core formulation or a space
of them: the problem's options are those of the backend, and `options` those of the run. A
core problem and a core formulation are a run of one point through [`traverse`](@ref), which
returns its value. A space computes the default study, the [`Combinatorial`](@ref) study of
`formulation`.
"""
function compute(problem::ParametricProblem, formulation::ScalarOrGrid{AbstractFormulation};
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    return compute(problem, Combinatorial(formulation); options)
end
function compute(problem::ParametricProblem{<:AbstractProblemDefinition},
        formulation::AbstractFormulation;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    return traverse(problem, formulation; options)
end

"""
    compute(problem::ParametricProblem, study::Combinatorial; options)

Compute every point of the problem's space with each formulation of the study, through
[`traverse`](@ref): the problem's options are those of the backend for every point, and
`options` those of the run. Return a [`ParametricResult`](@ref) with the failed points of the
study. A scalar problem and a scalar formulation are a study of one point.
"""
function compute(problem::ParametricProblem, study::Combinatorial;
        options::Union{NamedTuple, ComputationOptions} = ComputationOptions())
    traversed = traverse(problem, study; options)
    return ParametricResult(traversed.values, traversed.failures, traversed.gridpoints,
        traversed.details)
end
