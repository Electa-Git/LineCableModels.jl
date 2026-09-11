# Observation is disposable runtime state. Campaign/checkpoint records remain the
# authority for execution, recovery and scientific results.
const _CAMPAIGN_TRACKER = Base.ScopedValues.ScopedValue{Any}(nothing)
const _PROGRESS_TERMINAL_STATES = ("complete", "failed", "skipped", "interrupted")

"""Forward diagnostics while coordinating with the campaign's terminal display."""
struct CampaignLogger{L} <: Logging.AbstractLogger
    "The task's original logger."
    parent::L
end
Logging.min_enabled_level(logger::CampaignLogger)=Logging.min_enabled_level(logger.parent)
Logging.catch_exceptions(logger::CampaignLogger)=Logging.catch_exceptions(logger.parent)
Logging.shouldlog(logger::CampaignLogger, args...)=Logging.shouldlog(logger.parent, args...)
function Logging.handle_message(logger::CampaignLogger, args...; kwargs...)
    LineCableModels.with_progress_output() do
        Logging.handle_message(logger.parent, args...; kwargs...)
    end
end

"""
    CampaignProgress(root, ids, session; mode=:auto, io=stderr,
        clock=_progress_clock, terminal=io isa Base.TTY)

Keep disposable observation state for the selected benchmark IDs. Durations and
monotonic clock readings are in seconds; `rates` retains at most 32 observations
of the active work scope. Completed calculation durations teach advisory estimates
for comparable work in this invocation. Numerical results are never retained.
Operand completion summaries retain execution and compute-call wall durations
\\[s\\]. The renderer publishes each summary once, outside measured samples.
"""
mutable struct CampaignProgress
    "Campaign directory containing session snapshots."
    root::String
    "Execution session identity."
    session::String
    "Terminal selection: `:auto`, `:plain` or `:off`."
    mode::Symbol
    "Progress output stream."
    io::IO
    "Whether cursor-addressed rendering is enabled."
    redraw::Bool
    "Monotonic clock returning seconds."
    clock::Function
    "Invocation start in monotonic seconds."
    started::Float64
    "One bounded lifecycle record per selected benchmark."
    rows::Vector{Dict{String, Any}}
    "Current identifying scope, stage and absolute counters."
    active::Dict{String, Any}
    "Recent monotonic times and completed-work counts."
    rates::Vector{Tuple{Float64, Int}}
    "Identity of the scope described by `rates`."
    rate_scope::Any
    "Comparable historical benchmark wall durations, in seconds."
    estimates::Dict{String, Float64}
    "Protects observation updates and snapshot copying."
    state_lock::ReentrantLock
    "Serializes rendering, diagnostics and performance-span transitions."
    render_lock::ReentrantLock
    "Whether a controlled sample has suspended observation."
    paused::Bool
    "Whether the invocation has ended."
    closed::Bool
    "Whether optional terminal output failed."
    output_failed::Bool
    "Whether optional snapshot persistence failed."
    snapshot_failed::Bool
    "Last terminal refresh in monotonic seconds."
    last_render::Float64
    "Last snapshot publication in monotonic seconds."
    last_snapshot::Float64
    "Last cooperative UI yield in monotonic seconds."
    last_yield::Float64
    "Last rendered stage identity for plain-output throttling."
    last_stage::String
    "Number of terminal lines owned by the last display."
    lines::Int
    "Optional four-hertz rendering timer."
    timer::Union{Nothing, Timer}
end

_progress_clock() = time_ns() * 1e-9
function _progress_mode(value)
    value in (:auto, :plain, :off) ? value :
    throw(ArgumentError("progress must be :auto, :plain or :off"))
end

function _progress_job_count(calculation)
    problem, formulation = calculation.problem, calculation.formulation
    points = problem isa LineCableModels.ParametricProblem ? length(problem.space) : 1
    forms = formulation isa Gridspace ? length(formulation) :
            formulation isa
            Union{LineCableModels.LinearError, LineCableModels.Combinatorial} &&
            formulation.inner isa Union{Gridspace, AbstractVector} ?
            length(formulation.inner) : 1
    return points * forms
end

_progress_backend(::Gridspace{Target}) where {Target} = _progress_backend(Target)
function _progress_backend(formulation)
    formulation isa MonteCarlo && return _progress_backend(formulation.inner) * " / Monte Carlo"
    formulation isa LinearError && return _progress_backend(formulation.inner) * " / LEP"
    formulation isa LineCableModels.Combinatorial && return _progress_backend(formulation.inner)
    if formulation isa AbstractVector
        return join(unique(_progress_backend.(formulation)), ", ")
    end
    backend_type = formulation isa Type ? formulation : typeof(formulation)
    backend_type <: Engine.LineCableModelsFEM && return "FEM"
    backend_type <: PSCAD.PSCADFormulation && return "PSCAD"
    backend_type <: Engine.LineParametersFormulation && return "Owned"
    return string(nameof(backend_type))
end

function _progress_history_key(definition)
    return repr((case = definition.case_id,
        reference = _numerical_record(calculation_record(definition.reference)),
        candidate = _numerical_record(calculation_record(definition.candidate)),
        environment = _performance_identity()))
end

function CampaignProgress(root, ids, session; mode = :auto, io = stderr,
        clock = _progress_clock, terminal = io isa Base.TTY)
    _progress_mode(mode)
    rows = [Dict{String, Any}("id"=>string(id), "state"=>"pending",
                "reference"=>Dict{String, Any}("state"=>"pending", "completed"=>0,
                    "total"=>0, "known"=>false, "reused"=>0),
                "candidate"=>Dict{String, Any}("state"=>"pending", "completed"=>0,
                    "total"=>0, "known"=>false, "reused"=>0)) for id in ids]
    now = clock()
    return CampaignProgress(abspath(root), string(session), mode, io,
        terminal && mode === :auto && get(ENV, "TERM", "") != "dumb", clock, now,
        rows, Dict{String, Any}(), Tuple{Float64, Int}[], nothing, Dict{String, Float64}(),
        ReentrantLock(), ReentrantLock(), false, false, false, false,
        -Inf, -Inf, now, "", 0, nothing)
end

function _progress_workload(calculation, problem)
    formulation = calculation.formulation
    propagation = :core
    trials = 1
    if formulation isa Union{MonteCarlo, LinearError, LineCableModels.Combinatorial}
        propagation = nameof(typeof(formulation))
        if formulation isa MonteCarlo
            trials = formulation.options.trials
            trials === nothing && return nothing
        end
        formulation = formulation.inner
    end
    # Inspect declared formulations only; never materialize a sampled problem.
    if formulation isa Gridspace
        formulation = first(formulation)
    elseif formulation isa AbstractVector
        isempty(formulation) && return nothing
        formulation = first(formulation)
    end
    core_options = calculation.problem isa ParametricProblem ?
        merge(calculation.problem.options, calculation.options) : calculation.options
    fem_controls = formulation isa Engine.LineCableModelsFEM ?
        LineCableModels.computation_options(Engine.LineCableModelsFEM, core_options) : nothing
    execution = fem_controls !== nothing ?
        (workers=fem_controls.frequency_workers,
         threads=fem_controls.solver_threads,
         maps=fem_controls.plot_field_maps,
         mesh_policy=fem_controls.mesh_policy) : nothing
    options = Base.structdiff(core_options,
        (; on_result=get(core_options, :on_result, nothing),
           log_file=get(core_options, :log_file, nothing)))
    frequencies = problem.frequencies
    isempty(frequencies) && return nothing
    group = repr((backend=nameof(typeof(formulation)), propagation, execution,
        physical=hasproperty(formulation, :options) ? formulation.options : nothing,
        options=_selection_value(options),
        frequency_bounds=extrema(frequencies)))
    # Matrix size and frequency/trial counts are workload proxies, not a solver
    # complexity claim. Mesh-aware live observations supersede these estimates.
    work = Float64(_progress_job_count(calculation)) * length(frequencies) *
        max(1, length(problem.system.terminal_order))^2 * trials
    # PSCAD includes native project compilation and remote transfer per call;
    # treating that fixed setup as a cost per matrix entry badly understates
    # smaller cases. Learn whole-call duration for this opaque backend.
    formulation isa PSCAD.PSCADFormulation && (work = Float64(_progress_job_count(calculation)))
    return (; group, work)
end

function _progress_learn!(tracker)
    observations = Dict{String, Vector{Float64}}()
    for row in tracker.rows, role in ("reference", "candidate", "finalization")
        item = get(row, role, nothing)
        item === nothing && continue
        seconds = get(item, "observed_seconds", nothing)
        work = get(item, "work", 0.0)
        group = get(item, "group", nothing)
        if group !== nothing && seconds isa Real && isfinite(seconds) &&
           seconds >= 0 && work > 0
            push!(get!(Vector{Float64}, observations, group), seconds/work)
        end
    end
    for row in tracker.rows, role in ("reference", "candidate", "finalization")
        item = get(row, role, nothing)
        item === nothing && continue
        get(item, "estimate_source", "") == "history" && continue
        samples = get(observations, get(item, "group", ""), Float64[])
        isempty(samples) && continue
        item["estimate_seconds"] = median(samples) * item["work"]
        item["estimate_source"] = "comparable calculations"
        item["estimate_samples"] = length(samples)
    end
end

function _progress_declarations!(tracker, definitions)
    tracker === nothing && return
    for definition in definitions
        key = _progress_history_key(definition)
        path = joinpath(tracker.root, string(definition.id), "state.toml")
        previous = isfile(path) ? TOML.parsefile(path) : Dict{String, Any}()
        totals=map(role->_progress_job_count(getproperty(definition, role)), (
            :reference, :candidate))
        lock(tracker.state_lock) do
            row = only(filter(row->row["id"]==string(definition.id), tracker.rows))
            row["case"] = string(definition.case_id)
            row["timing_key"] = key
            for (role, total) in zip((:reference, :candidate), totals)
                row[string(role)]["total"] = total
                row[string(role)]["known"] = true
                row[string(role)]["backend"] = _progress_backend(getproperty(definition, role).formulation)
                workload = _progress_workload(getproperty(definition, role), definition.model.nominal_problem)
                if workload !== nothing
                    row[string(role)]["group"] = workload.group
                    row[string(role)]["work"] = workload.work
                end
            end
            row["finalization"] = Dict{String, Any}(
                "group"=>repr((collection=definition.collection,
                    report=definition.comparison_settings,
                    performance=get(definition.tolerances, :performance, nothing))),
                "work"=>Float64(sum(totals) * length(definition.model.nominal_problem.frequencies) *
                    max(1, length(definition.model.port_order))^2))
            elapsed = get(previous, "fresh_wall_seconds", get(previous, "wall_seconds", nothing))
            previous_key=get(previous, "fresh_timing_key", get(previous, "timing_key", nothing))
            if previous_key == key && elapsed isa Real && isfinite(elapsed) &&
               elapsed > 0 &&
               (haskey(previous, "fresh_wall_seconds") || !get(previous, "reused", false))
                tracker.estimates[string(definition.id)] = Float64(elapsed)
                for role in ("reference", "candidate")
                    estimate=get(previous, "fresh_$(role)_seconds", nothing)
                    if estimate isa Real && isfinite(estimate) && estimate>=0
                        row[role]["estimate_seconds"] = Float64(estimate)
                        row[role]["estimate_source"] = "history"
                        row[role]["observed_seconds"] = Float64(estimate)
                    end
                end
                if all(haskey(row[role], "observed_seconds") for role in ("reference", "candidate"))
                    row["finalization"]["observed_seconds"] = max(0.0, elapsed -
                        row["reference"]["observed_seconds"] - row["candidate"]["observed_seconds"])
                end
            end
        end
    end
    lock(tracker.state_lock) do
        _progress_learn!(tracker)
    end
end

function _progress_event!(tracker::CampaignProgress, event)
    if get(event, :kind, nothing) === :native_console
        # Native libraries write directly to stdout. Keep explicit verbose output
        # intact by using append-only progress for the remainder of this run.
        lock(tracker.render_lock) do
            tracker.redraw=false
            tracker.lines=0
        end
        return nothing
    elseif get(event, :kind, nothing) === :output_begin
        lock(tracker.render_lock)
        # Leave the previous display above the diagnostic. The next refresh starts
        # a new display below it instead of erasing a warning/error.
        tracker.lines=0
        return nothing
    elseif get(event, :kind, nothing) === :output_end
        unlock(tracker.render_lock)
        return nothing
    end
    lock(tracker.state_lock) do
        (tracker.closed || tracker.paused) && return
        now = tracker.clock()
        id = string(get(event, :benchmark, get(tracker.active, "benchmark", "")))
        row_index = findfirst(row->row["id"]==id, tracker.rows)
        row = row_index === nothing ? nothing : tracker.rows[row_index]
        kind = get(event, :kind, :work)
        if kind === :benchmark
            empty!(tracker.active)
            empty!(tracker.rates)
            tracker.rate_scope = nothing
            tracker.active["benchmark"] = id
            tracker.active["started"] = now
            tracker.active["phase_started"] = now
            row === nothing || (row["state"] = string(event.state))
        end
        if row !== nothing && haskey(event, :attempt)
            row["attempt"] = string(event.attempt)
        end
        if kind === :operand && row !== nothing
            tracker.active["phase_started"] = now
            role = row[string(event.role)]
            haskey(event, :backend) && (role["backend"] = string(event.backend))
            role["state"] = string(event.state)
            if event.state === :complete
                role["announced"] = get(role, "announced", false)
                haskey(event, :seconds) && (role["execution_seconds"] = event.seconds)
                haskey(event, :compute_seconds) && event.compute_seconds !== nothing &&
                    (role["compute_seconds"] = event.compute_seconds)
                role["saved_result"] = get(event, :reused, false)
                role["completed"] = role["total"]
                role["timing_reused"] = get(event, :timing_reused, false)
                role["reused"] = get(event, :reused, false) ? role["total"] :
                                 clamp(get(event, :jobs_reused, 0), 0, role["total"])
                if role["reused"] == 0 && !get(event, :timing_reused, false) && haskey(event, :seconds)
                    role["observed_seconds"] = event.seconds
                    _progress_learn!(tracker)
                end
                haskey(tracker.active,"total") &&
                    (tracker.active["completed"]=tracker.active["total"])
                empty!(tracker.rates)
                tracker.rate_scope=nothing
            elseif event.state === :running
                role["started"] = now
                empty!(tracker.rates)
                tracker.rate_scope = nothing
                for name in ("unit", "completed", "total", "attempts", "rejected",
                    "workers", "queued", "recovered", "point", "formulation",
                    "backend", "estimate_throughput", "remaining_seconds", "eta_updated", "eta_source",
                    "seconds", "compute_seconds", "timing_reused", "reused", "jobs_reused")
                    pop!(tracker.active, name, nothing)
                end
            end
        end
        if kind === :benchmark && event.state === :complete && row !== nothing &&
           haskey(row, "finalization") && haskey(event, :finalization_seconds) &&
           !get(event, :reused, false) &&
           !any(get(row[r], "timing_reused", false) for r in ("reference", "candidate"))
            row["finalization"]["observed_seconds"] = event.finalization_seconds
            _progress_learn!(tracker)
        end
        if kind === :benchmark && event.state in (:failed, :interrupted)
            for role in ("reference", "candidate")
                row === nothing && continue
                item = row[role]
                if item["state"] == "running"
                    item["execution_seconds"] = max(0.0, now-get(item, "started", now))
                    item["announced"] = false
                end
                item["state"] = item["state"] == "running" ? string(event.state) :
                                item["state"] == "pending" ? "skipped" : item["state"]
            end
        end
        for (name, value) in pairs(event)
            if value === nothing
                pop!(tracker.active, string(name), nothing)
                continue
            end
            tracker.active[string(name)] = value isa Symbol ? string(value) : value
        end
        haskey(event, :remaining_seconds) && (tracker.active["eta_updated"] = now)
        active_role = get(tracker.active, "role", "")
        if row !== nothing && active_role in ("reference", "candidate")
            operand=row[active_role]
            haskey(operand, "backend") && (tracker.active["backend"] = operand["backend"])
            count = get(event, :jobs_completed, nothing)
            if count === nothing && get(event, :unit, nothing) === :formulations
                count=get(event, :completed, nothing)
            elseif count === nothing && get(event, :unit, nothing) === :points
                count=get(event, :completed, 0)*get(event, :formulations, 1)
            end
            count === nothing || (operand["completed"]=clamp(count, 0, operand["total"]))
        end
        scope = (id, get(event, :role, get(tracker.active, "role", "")),
            get(event, :point, get(tracker.active, "point", 0)),
            get(event, :formulation, get(tracker.active, "formulation", 0)),
            get(event, :unit, nothing), get(event, :stage, nothing))
        if haskey(event, :completed) && haskey(event, :unit)
            if scope != tracker.rate_scope ||
               (!isempty(tracker.rates) && event.completed < last(tracker.rates)[2])
                empty!(tracker.rates)
                tracker.rate_scope = scope
            end
            # Rates describe actual completed work, not heartbeat updates.
            if isempty(tracker.rates) || event.completed > last(tracker.rates)[2]
                push!(tracker.rates, (now, Int(event.completed)))
                length(tracker.rates) > 32 && popfirst!(tracker.rates)
            end
        end
    end
    # Allow the UI task to run during long single-threaded MC workloads. There is
    # no rendering/IO here, and quiet samples have no receiver at all.
    now = tracker.clock()
    if now - tracker.last_yield >= 0.25
        tracker.last_yield = now
        yield()
    end
    return nothing
end

function _progress_snapshot(tracker)
    lock(tracker.state_lock) do
        now = tracker.clock()
        rows = deepcopy(tracker.rows)
        for row in rows
            pop!(row, "timing_key", nothing)
            for role in ("reference", "candidate", "finalization")
                haskey(row, role) && pop!(row[role], "group", nothing)
            end
        end
        counts = Dict(state=>count(row->row["state"]==state, rows)
        for state in ("pending", "running", _PROGRESS_TERMINAL_STATES...))
        remaining = 0.0
        unknown = 0
        for row in rows
            row["state"] in _PROGRESS_TERMINAL_STATES && continue
            estimate = get(tracker.estimates, row["id"], nothing)
            if estimate === nothing
                uncovered = false
                for role in ("reference", "candidate", "finalization")
                    item = get(row, role, Dict{String, Any}())
                    get(item, "state", "pending") == "complete" && continue
                    seconds = get(item, "estimate_seconds", nothing)
                    active = row["state"] == "running" &&
                        (role == "finalization" ?
                         all(row[r]["state"] == "complete" for r in ("reference", "candidate")) :
                         item["state"] == "running")
                    elapsed = active ? now-get(tracker.active, "phase_started", now) : 0.0
                    if seconds isa Real && seconds >= elapsed
                        remaining += seconds-elapsed
                    else
                        uncovered = true
                    end
                end
                unknown += uncovered
            else
                decomposed=all(haskey(row[role], "estimate_seconds")
                for role in ("reference", "candidate"))
                if row["state"] == "running" && decomposed
                    estimate-=sum(
                        row[role]["estimate_seconds"]
                        for role in ("reference", "candidate")
                        if row[role]["state"] == "complete";
                        init = 0.0)
                end
                elapsed = row["state"] == "running" ?
                          now-get(tracker.active,
                    decomposed ? "phase_started" : "started", now) : 0.0
                if elapsed >= estimate
                    unknown += 1
                else
                    remaining += estimate-elapsed
                end
            end
        end
        stage_eta = -1.0
        eta_scope = "stage"
        if get(tracker.active, "stage", "") == "solving"
            estimate = get(tracker.active, "remaining_seconds", -1.0)
            elapsed = now-get(tracker.active, "eta_updated", now)
            estimate > elapsed && (stage_eta = estimate-elapsed)
        end
        total = get(tracker.active, "total", nothing)
        if stage_eta < 0 && length(tracker.rates) >= 3 && total isa Integer &&
           get(tracker.active, "estimate_throughput", true)
            start, stop = first(tracker.rates), last(tracker.rates)
            done = stop[2]-start[2]
            elapsed = now-start[1]
            if done > 0 && elapsed > 0 && total > stop[2]
                stage_eta = (total-stop[2])*elapsed/done
            end
        end
        if stage_eta < 0
            index=findfirst(row->row["id"]==get(tracker.active,"benchmark",""),rows)
            role=get(tracker.active,"role","")
            if index !== nothing && role in ("reference","candidate")
                operand=rows[index][role]
                estimate=get(operand,"estimate_seconds",nothing)
                elapsed=now-get(tracker.active,"phase_started",now)
                if operand["state"]=="running" && estimate isa Real && estimate>elapsed
                    stage_eta=estimate-elapsed
                    eta_scope="operand"
                end
            end
        end
        return Dict{String, Any}("schema"=>1, "session"=>tracker.session,
            "updated_unix_seconds"=>time(), "elapsed_seconds"=>max(0.0, now-tracker.started),
            "closed"=>tracker.closed, "measurement_active"=>tracker.paused,
            "benchmarks"=>rows, "counts"=>counts, "active"=>copy(tracker.active),
            "stage_eta_seconds"=>stage_eta,"active_eta_scope"=>eta_scope,
            "campaign_eta_seconds"=>unknown == 0 ? remaining : -1.0,
            "estimated_remaining_seconds"=>remaining,
            "unestimated_benchmarks"=>unknown)
    end
end

function _progress_duration(seconds)
    seconds isa Real && isfinite(seconds) && seconds >= 0 || return "estimating"
    minutes = round(Int, seconds/60)
    return seconds < 60 ? "$(round(Int,seconds))s" :
           minutes < 60 ? "$(minutes)m" : "$(div(minutes,60))h $(rem(minutes,60))m"
end

function _progress_lines(snapshot, width)
    counts, active = snapshot["counts"], snapshot["active"]
    total = length(snapshot["benchmarks"])
    done = sum(counts[state] for state in _PROGRESS_TERMINAL_STATES)
    n = clamp(width-43, 8, 28)
    fill = total == 0 ? n : clamp(floor(Int, n*done/total), 0, n)
    lines = [
        "Selected ["*repeat("=", fill)*repeat(" ", n-fill)*"] $done/$total benchmarks finished",
        "OK $(counts["complete"]) | Failed $(counts["failed"]) | Running $(counts["running"]) | Pending $(counts["pending"])"]
    counts["skipped"]+counts["interrupted"] > 0 && push!(lines,
        "Skipped $(counts["skipped"]) | Interrupted $(counts["interrupted"])")
    if !isempty(active)
        index=findfirst(row->row["id"]==get(active, "benchmark", ""), snapshot["benchmarks"])
        push!(lines,
            "Active $(something(index,"?"))/$total | "*join(
                filter(!isempty, [string(get(active, key, ""))
                                  for key in ("role", "backend", "stage")]), " | "))
        push!(lines, string(get(active, "benchmark", "Preparing declarations")))
        scope=["$key $(active[key])"
               for key in ("point", "formulation") if haskey(active, key)]
        isempty(scope) || push!(lines, join(scope, " | "))
        if snapshot["measurement_active"]
            push!(lines,
                "Performance sample $(get(active,"sample","warmup"))/$(get(active,"samples","?")) — display paused")
        elseif haskey(active, "unit")
            line = "$(get(active,"unit","work")): $(get(active,"completed",0))/$(get(active,"total","?"))"
            haskey(active, "attempts") &&
                (line *= " | Attempts $(active["attempts"]) | Rejected $(get(active,"rejected",0))")
            haskey(active, "workers") &&
                (line *= " | Workers $(active["workers"]) | Queued $(get(active,"queued",0)) | Recovered $(get(active,"recovered",0))")
            push!(lines, line)
        end
    end
    jobs = [row[role] for row in snapshot["benchmarks"]
            for role in ("reference", "candidate")]
    remaining=sum(job["total"]-job["completed"] for job in jobs
        if job["state"] in ("pending","running");init=0)
    push!(lines,
        "Calculation jobs $(sum(job["completed"] for job in jobs;init=0))/$(sum(job["total"] for job in jobs;init=0)) complete | $remaining remaining | Reused $(sum(job["reused"] for job in jobs;init=0))")
    unresolved=count(job->!get(job, "known", true), jobs)
    unresolved>0 && push!(lines, "Job totals unresolved for $unresolved operands")
    approximate(seconds) = seconds isa Real && isfinite(seconds) && seconds >= 0 ?
        "~" * _progress_duration(seconds) : "estimating"
    push!(lines,
        "Elapsed $(_progress_duration(snapshot["elapsed_seconds"])) | $(uppercasefirst(get(snapshot,"active_eta_scope","stage"))) ETA $(approximate(snapshot["stage_eta_seconds"])) | Campaign ETA $(approximate(snapshot["campaign_eta_seconds"]))")
    if snapshot["unestimated_benchmarks"] > 0
        covered = get(snapshot, "estimated_remaining_seconds", 0.0)
        prefix = covered > 0 ? "Estimated portion ~" * _progress_duration(covered) * "; " : ""
        push!(lines, prefix * "$(snapshot["unestimated_benchmarks"]) benchmarks include unestimated work")
    end
    return _progress_fit_lines(lines, width)
end

function _progress_fit_lines(lines, width)
    return map(lines) do line
        columns=0
        join(Iterators.takewhile(line) do character
            columns+=textwidth(character)
            columns<=max(8, width)
        end)
    end
end

function _progress_completion_lines(row, role, width)
    item = row[role]
    duration(seconds) = seconds isa Real && isfinite(seconds) && seconds >= 0 ?
        "$(round(seconds; sigdigits=5))s" : "unavailable"
    compute = get(item, "saved_result", false) ? "not run (saved result)" :
        duration(get(item, "compute_seconds", nothing))
    recovery = haskey(item, "timing_reused") ? (item["timing_reused"] ? "yes" : "no") : "unknown"
    state = item["state"] == "complete" ? "Completed" : uppercasefirst(item["state"])
    return _progress_fit_lines([
        "$state | $role | $(get(item, "backend", "unknown backend"))",
        row["id"],
        "Execution wall $(duration(get(item, "execution_seconds", nothing))) | Compute-call wall $compute",
        "Calculation jobs $(item["completed"])/$(item["total"]) | Reused $(item["reused"]) | Recovery $recovery",
    ], width)
end

function _progress_paint!(tracker; force = false)
    lock(tracker.render_lock) do
        tracker.paused && !force && return
        now = tracker.clock()
        snapshot = _progress_snapshot(tracker)
        completions = [(index, role) for (index, row) in enumerate(snapshot["benchmarks"])
            for role in ("reference", "candidate")
            if row[role]["state"] in ("complete", "failed", "interrupted") &&
                !get(row[role], "announced", false)]
        active=get(snapshot, "active", Dict())
        stage=join(string(get(active, key, "")) for key in ("benchmark", "role", "stage"))
        if !tracker.output_failed && (force || !isempty(completions) ||
            now-tracker.last_render >=
                (tracker.redraw ? 0.25 : stage != tracker.last_stage ? 1.0 : 10.0))
            try
                width = max(20, displaysize(tracker.io)[2])
                lines = _progress_lines(snapshot, width)
                if tracker.redraw && tracker.lines > 0
                    print(tracker.io, "\e[$(tracker.lines)A")
                    if !isempty(completions)
                        # Commit summaries above the live display. Future redraws
                        # own only the dashboard and cannot erase completions.
                        for _ in 1:tracker.lines
                            print(tracker.io, "\r\e[2K\n")
                        end
                        print(tracker.io, "\e[$(tracker.lines)A")
                        tracker.lines = 0
                    end
                end
                for (index, role) in completions
                    for line in _progress_completion_lines(snapshot["benchmarks"][index], role, width)
                        println(tracker.io, line)
                    end
                end
                for line in lines
                    tracker.redraw && print(tracker.io, "\r\e[2K")
                    println(tracker.io, line)
                end
                if tracker.redraw
                    for _ in 1:max(0, tracker.lines - length(lines))
                        print(tracker.io, "\r\e[2K\n")
                    end
                    extra = max(0, tracker.lines-length(lines))
                    extra > 0 && print(tracker.io, "\e[$(extra)A")
                end
                flush(tracker.io)
                lock(tracker.state_lock) do
                    for (index, role) in completions
                        tracker.rows[index][role]["announced"] = true
                        snapshot["benchmarks"][index][role]["announced"] = true
                    end
                end
                tracker.lines = length(lines)
                tracker.last_render = now
                tracker.last_stage = stage
            catch exception
                exception isa InterruptException && rethrow()
                tracker.output_failed = true
                try
                    @warn "Campaign progress output disabled" exception
                catch diagnostic_error
                    diagnostic_error isa InterruptException && rethrow()
                end
            end
        end
        if !tracker.snapshot_failed && (force || now-tracker.last_snapshot >= 1.0)
            try
                path=joinpath(tracker.root, "sessions", tracker.session*".progress.toml")
                mkpath(dirname(path))
                _write_toml(path, snapshot)
                tracker.last_snapshot=now
            catch exception
                exception isa InterruptException && rethrow()
                tracker.snapshot_failed=true
                try
                    @warn "Campaign progress snapshot disabled" exception
                catch diagnostic_error
                    diagnostic_error isa InterruptException && rethrow()
                end
            end
        end
    end
    return nothing
end

function _with_campaign_progress(f, directory, ids, session; progress = :auto)
    _progress_mode(progress)
    existing = _CAMPAIGN_TRACKER[]
    existing === nothing || return f(existing)
    progress === :off && return f(nothing)
    tracker = CampaignProgress(directory, ids, session; mode = progress)
    receiver = event->_progress_event!(tracker, event)
    tracker.timer = Timer(0.25; interval = 0.25) do _
        tracker.closed || _progress_paint!(tracker)
    end
    return Base.ScopedValues.with(_CAMPAIGN_TRACKER=>tracker) do
        try
            LineCableModels.with_progress(receiver) do
                Logging.with_logger(CampaignLogger(Logging.current_logger())) do
                    f(tracker)
                end
            end
        catch exception
            lock(tracker.state_lock) do
                for row in tracker.rows
                    row["state"] == "running" &&
                        (row["state"]=exception isa InterruptException ? "interrupted" :
                                      "failed")
                    row["state"] == "pending" && (row["state"]="skipped")
                    for role in ("reference", "candidate")
                        row[role]["state"] == "pending" && (row[role]["state"]="skipped")
                        row[role]["state"] == "running" && (row[role]["state"]=row["state"])
                    end
                end
            end
            rethrow()
        finally
            close(tracker.timer)
            lock(tracker.state_lock) do
                tracker.closed=true
                tracker.paused=false
            end
            _progress_paint!(tracker; force = true)
        end
    end
end

function _performance_span(f; sample, samples, role)
    tracker = _CAMPAIGN_TRACKER[]
    receiver = LineCableModels.progress_receiver()
    LineCableModels.report_progress(receiver,
        (kind = :work, stage = :performance, role, sample, samples))
    if tracker !== nothing
        lock(tracker.render_lock) do
            lock(tracker.state_lock) do
                tracker.paused=true
            end
            _progress_paint!(tracker; force = true)
        end
    end
    try
        return LineCableModels.with_performance_sample(f)
    finally
        tracker === nothing || lock(tracker.render_lock) do
            lock(tracker.state_lock) do
                tracker.paused=false
            end
        end
    end
end

function _campaign_watch_lines(root, width; now = time())
    rows=campaign_status(root; verify = false)
    output=["Campaign $(basename(root)): "*join(
        ["$state $(count(row->row.state===state,rows))"
         for state in
             (:complete, :running, :pending, :failed, :interrupted)], " | ")]
    # Read only sessions that currently own a running attempt. Old snapshots are
    # neither numerical evidence nor a reason to scan an ever-growing history.
    for row in rows
        row.state === :running || continue
        state=TOML.parsefile(joinpath(root, row.id, "state.toml"))
        session=get(state, "session", nothing)
        session isa String && basename(session)==session || continue
        path=joinpath(root, "sessions", session*".progress.toml")
        isfile(path) || continue
        snapshot=try
            TOML.parsefile(path)
        catch
            ;
            continue
        end
        get(snapshot, "schema", 0)==1 || continue
        get(snapshot, "closed", true) && continue
        active=get(snapshot, "active", Dict())
        get(active, "benchmark", nothing)==row.id || continue
        get(snapshot, "session", nothing)==session || continue
        get(state, "attempt", nothing)==get(active, "attempt", nothing) || continue
        age=max(0.0, now-get(snapshot, "updated_unix_seconds", 0.0))
        append!(output, _progress_lines(snapshot, width))
        age>5 && push!(output,
            get(snapshot, "measurement_active", false) ?
            "Observation paused during performance sample" :
            "Last observation $(round(Int,age))s ago; liveness not inferred")
    end
    return output
end

"""
    watch_campaign(directory; io=stderr, interval=1.0)

Display lightweight campaign state and live progress snapshots until interrupted.
`interval` is the polling interval in seconds. This read-only observer never loads
executable declarations or numerical results. Ctrl-C stops watching, not solving.
"""
function watch_campaign(directory; io = stderr, interval = 1.0)
    isfinite(interval) && interval > 0 ||
        throw(ArgumentError("watch interval must be positive"))
    root=abspath(directory)
    terminal=io isa Base.TTY && get(ENV, "TERM", "") != "dumb"
    lines=0
    try
        while true
            output=_campaign_watch_lines(root, max(20, displaysize(io)[2]))
            terminal && lines>0 && print(io, "\e[$(lines)A")
            for line in output
                terminal && print(io, "\r\e[2K")
                println(io, line)
            end
            for _ in 1:(terminal ? max(0, lines - length(output)) : 0)
                print(io, "\r\e[2K\n")
            end
            terminal && lines>length(output) && print(io, "\e[$(lines-length(output))A")
            lines=length(output)
            flush(io)
            sleep(interval)
        end
    catch exception
        exception isa InterruptException || rethrow()
    end
    return nothing
end

export watch_campaign
