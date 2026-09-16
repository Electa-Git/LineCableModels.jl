# Disposable observations never determine execution or recovery policy.
const _CAMPAIGN_TRACKER = Base.ScopedValues.ScopedValue{Any}(nothing)
const _PROGRESS_TERMINAL_STATES = ("complete", "failed", "skipped", "interrupted")
_progress_clock() = time_ns() * 1e-9
_progress_mode(mode) = mode in (:auto, :plain, :off) ? mode :
    throw(ArgumentError("progress must be :auto, :plain or :off"))
_valid_seconds(x) = x isa Real && isfinite(x) && x >= 0

"""
    CampaignProgress(root, ids, session; mode=:auto, io=stderr,
        clock=_progress_clock, wall_clock=time, selected=false)

Collect selected benchmark outcomes, accepted scans, and advisory duration
budgets. Clocks return seconds. Memory is bounded by selected benchmarks and
active scopes; execution owns all scientific and timing records. Publication
occurs at outer boundaries; this object has no timer or terminal renderer.
"""
Base.@kwdef mutable struct CampaignProgress
    "Campaign directory."
    root::String
    "Invocation identity."
    session::String
    "Execution output policy."
    mode::Symbol = :auto
    "Optional execution output."
    io::IO = stderr
    "Monotonic clock in seconds."
    clock::Function = _progress_clock
    "Wall clock for cross-process snapshot age in seconds."
    wall_clock::Function = time
    "Invocation start in local monotonic seconds."
    started::Float64
    "Final local monotonic time, or NaN while open."
    stopped::Float64 = NaN
    "Whether the invocation selects part of a larger inventory."
    selected::Bool = false
    "One compact record per selected benchmark."
    rows::Vector{Dict{String,Any}}
    "Selected benchmark lookup."
    indices::Dict{String,Int}
    "Current operand and primary scan scope."
    active::Dict{String,Any} = Dict{String,Any}()
    "Open scan scopes only."
    scopes::Set{Int} = Set{Int}()
    "Highest scope generation opened."
    generation::Int = 0
    "Latest source sequence accepted."
    sequence::Int = 0
    "EWMA duration evidence keyed only by backend, mode and timing scope."
    evidence::Dict{Tuple{String,String,String},Float64} = Dict{Tuple{String,String,String},Float64}()
    "Synchronizes observations and small payload copies."
    state_lock::ReentrantLock = ReentrantLock()
    "Serializes replacements and controlled-call transitions."
    publication_lock::ReentrantLock = ReentrantLock()
    "Whether optional observation is suspended."
    paused::Bool = false
    "Whether the selected invocation is closed."
    closed::Bool = false
    "Whether execution exhausted its selection or stopped early."
    termination::String = "running"
    "Whether optional file publication failed."
    snapshot_failed::Bool = false
    "Whether optional console output failed."
    output_failed::Bool = false
    "Last publication eligibility check that produced a snapshot."
    last_snapshot::Float64 = -Inf
    "Last plain execution summary."
    last_plain::Float64 = -Inf
    "Last actual work change in local monotonic seconds."
    work_time::Float64
    "Published snapshot revision."
    revision::Int = 0
end

function CampaignProgress(root, ids, session; mode=:auto, io=stderr,
        clock=_progress_clock, wall_clock=time, selected=false)
    _progress_mode(mode)
    rows = [Dict{String,Any}("id"=>string(id), "state"=>"pending",
        "reference"=>Dict{String,Any}("state"=>"pending", "total"=>-1, "completed"=>0, "reused"=>0),
        "candidate"=>Dict{String,Any}("state"=>"pending", "total"=>-1, "completed"=>0, "reused"=>0)) for id in ids]
    now = clock()
    return CampaignProgress(; root=abspath(root), session=string(session), mode, io,
        clock, wall_clock, selected, started=now, work_time=now, rows,
        indices=Dict(row["id"]=>i for (i,row) in enumerate(rows)))
end

function _duration_evidence!(tracker, key, seconds)
    _valid_seconds(seconds) && seconds > 0 || return
    previous = get(tracker.evidence, key, seconds)
    tracker.evidence[key] = 0.25seconds + 0.75previous
end

# No geometry, option serialization, source hashes, or frequency detail in keys.
function _duration_estimate(tracker, operand, scope)
    backend, mode = get(operand,"backend",""), get(operand,"mode","")
    exact = get(tracker.evidence, (backend,mode,scope), -1.0)
    exact >= 0 && return (exact, false)
    pooled = [v for (k,v) in tracker.evidence if k[3] == scope]
    return isempty(pooled) ? (-1.0, true) : (sum(pooled)/length(pooled), true)
end

function _read_progress_toml(path)
    try
        return TOML.parsefile(path)
    catch error
        error isa InterruptException && rethrow()
        return Dict{String,Any}()
    end
end

function _progress_declarations!(tracker, definitions, plan=nothing)
    tracker === nothing && return
    for (index,definition) in enumerate(definitions)
        row = tracker.rows[tracker.indices[string(definition.id)]]
        get(row,"resolved",false) && continue
        # The execution owner has already established compatibility in its plan.
        item = plan === nothing ? nothing : plan[index]
        previous = item === nothing || item.attempt === nothing ? Dict{String,Any}() :
            _read_progress_toml(joinpath(tracker.root,string(definition.id),"state.toml"))
        descriptors = map(role->_calculation_progress(getproperty(definition,role)), (:reference,:candidate))
        settings = _benchmark_performance_settings(definition.tolerances)
        lock(tracker.state_lock) do
            row["resolved"] = true
            row["case"] = string(definition.case_id)
            row["performance_limit"] = settings === nothing ? 0.0 : settings.seconds
            for (role,descriptor) in zip(("reference","candidate"),descriptors)
                operand = row[role]
                merge!(operand, Dict(string(k)=>v for (k,v) in pairs(descriptor)))
                operand["performance_left"] = settings === nothing ? 0 : settings.samples
                operand["warmup_left"] = settings !== nothing && descriptor.warmup
                operand["performance_elapsed"] = 0.0
                if !get(previous,"reused",true)
                    seconds = get(previous,"fresh_$(role)_seconds",-1.0)
                    if _valid_seconds(seconds) && seconds > 0
                        operand["estimate"] = seconds
                        _duration_evidence!(tracker,(operand["backend"],operand["mode"],"operation"),seconds)
                    end
                end
            end
            whole = get(previous,"wall_seconds",-1.0)
            if !get(previous,"reused",true) && _valid_seconds(whole) && whole > 0
                row["whole_estimate"] = whole
                _duration_evidence!(tracker,("","","benchmark"),whole)
            end
        end
    end
end

function _scan_rate!(tracker, now; final=false)
    a=tracker.active
    get(a,"rate_eligible",false) || return
    completed=get(a,"completed",0)-get(a,"reused",0)
    anchor=get(a,"rate_count",0)
    finished=get(a,"accepted_time",now)
    elapsed=finished-get(a,"rate_time",finished)
    completed > anchor && elapsed > 0 && (final || elapsed >= 1.0) || return
    seconds=elapsed/(completed-anchor)
    a["scan_seconds"] = 0.25seconds + 0.75get(a,"scan_seconds",seconds)
    a["rate_count"],a["rate_time"] = completed,finished
    row=tracker.rows[tracker.indices[a["benchmark"]]]
    operand=row[a["role"]]
    # Effective throughput is only used by this same operation/batch arrangement.
    operand["scan_seconds"] = a["scan_seconds"]
    operand["opaque_interval"] = false
    if get(a,"batch",1)==1
        _duration_evidence!(tracker,(operand["backend"],operand["mode"],"scan"),seconds)
    end
end

function _progress_event!(tracker::CampaignProgress, event)
    force = get(event,:kind,:work) in (:benchmark,:operand)
    accepted = lock(tracker.state_lock) do
        (tracker.closed || tracker.paused) && return
        get(event,:session,tracker.session)==tracker.session || return
        seq=get(event,:sequence,tracker.sequence+1)
        seq > tracker.sequence || return
        id=string(get(event,:benchmark,get(tracker.active,"benchmark","")))
        index=get(tracker.indices,id,0)
        index==0 && return
        row=tracker.rows[index]
        kind=get(event,:kind,:work)
        attempt=string(get(event,:attempt,get(row,"attempt","")))
        if kind !== :benchmark
            get(row,"state","pending")=="running" || return
            attempt==get(row,"attempt",attempt) || return
            kind === :operand || string(get(event,:role,get(tracker.active,"role","")))==get(tracker.active,"role","") || return
        elseif row["state"] in _PROGRESS_TERMINAL_STATES ||
                row["state"]=="running" && attempt!=get(row,"attempt",attempt)
            return
        end
        if kind === :operand
            role=string(event.role)
            role in ("reference","candidate") || return
            if event.state === :running
                row[role]["state"]=="pending" || return
            else
                role==get(tracker.active,"role","") && row[role]["state"]=="running" || return
            end
        end
        scope=get(event,:scan_scope,0)
        if kind === :scan_start
            scope > tracker.generation || return
            parent=get(event,:scan_parent,0)
            parent==0 || parent in tracker.scopes || return
            tracker.generation=scope
            push!(tracker.scopes,scope)
        elseif scope != 0 && !(scope in tracker.scopes)
            return
        end
        child=get(event,:child,0)
        child > 0 && child <= get(tracker.active,"children_completed",0) && return
        tracker.sequence=seq
        now=tracker.clock()
        a=tracker.active
        changed=false
        if kind === :benchmark
            (event.state === :running || get(a,"benchmark","")!=id) && empty!(a)
            empty!(tracker.scopes)
            a["benchmark"]=id; a["attempt"]=attempt; a["stage"]=string(get(event,:stage,event.state))
            row["attempt"]=attempt
            row["state"]=string(event.state)
            row["execution_state"]=string(event.state)
            if event.state === :running
                row["started"]=now
            else
                for role in ("reference","candidate")
                    operand=row[role]
                    if operand["state"]=="running"
                        operand["state"]=string(event.state)
                    elseif operand["state"]=="pending"
                        operand["state"]="skipped"
                    end
                end
                seconds=get(event,:seconds,-1.0)
                if event.state === :complete && !get(event,:reused,false) && !any(get(row[r],"recovered_execution",false) for r in ("reference","candidate")) && _valid_seconds(seconds)
                    _duration_evidence!(tracker,("","","benchmark"),seconds)
                    row["whole_estimate"]=seconds
                end
            end
            changed=true
        elseif kind === :operand
            role=string(event.role); operand=row[role]
            if event.state === :running
                empty!(a); empty!(tracker.scopes)
                merge!(a,Dict("benchmark"=>id,"attempt"=>attempt,"role"=>role,
                    "stage"=>"preparing","completed"=>0,"total"=>get(operand,"total",1),
                    "reused"=>0,"scope"=>0,"phase_started"=>now))
                operand["started"]=now
                for key in (:backend,:mode)
                    haskey(event,key) && (operand[string(key)]=string(event[key]))
                end
                a["backend"]=get(operand,"backend",""); a["mode"]=get(operand,"mode","")
            else
                _scan_rate!(tracker,now;final=true)
                empty!(tracker.scopes)
                seconds=get(event,:seconds,-1.0)
                if event.state === :complete && !get(event,:timing_reused,false) && !get(operand,"recovered_execution",false) && _valid_seconds(seconds)
                    operand["estimate"]=seconds
                    _duration_evidence!(tracker,(get(operand,"backend",""),get(operand,"mode",""),"operation"),seconds)
                end
                a["stage"]="validating"
                row["finalization_started"]=now
            end
            operand["state"]=string(event.state)
            changed=true
        elseif kind === :scan_start && get(event,:scan_parent,0)==0
            for key in ("scan_seconds","accepted_time","frequencies_completed","frequencies_total","workers","rejected","point","formulation")
                pop!(a,key,nothing)
            end
            merge!(a,Dict("scope"=>scope,"completed"=>0,"total"=>something(get(event,:total,nothing),-1),
                "reused"=>0,"rate_count"=>0,"rate_time"=>now,"rate_eligible"=>true,
                "batch"=>get(event,:batch,1),"children_completed"=>0,"stage"=>"solving"))
            changed=true
        elseif kind === :scan_end
            scope==get(a,"scope",0) && _scan_rate!(tracker,now;final=true)
            delete!(tracker.scopes,scope)
        elseif kind === :scan || kind === :scan_result && get(a,"scope",0)==0
            if kind === :scan && scope != get(a,"scope",0)
                return
            end
            count=get(event,:completed,get(a,"completed",0))
            count >= get(a,"completed",0) || return
            role=get(a,"role","")
            operand=row[role]
            delta=count-get(a,"completed",0)
            if delta>0 && get(a,"child_recovered",false)
                a["reused"]=get(a,"reused",0)+delta
                a["child_recovered"]=false
            end
            operand["completed"]=get(operand,"completed",0)+delta
            if delta>0
                operand["count_updated"]=now
                a["accepted_time"]=now
            end
            for key in (:completed,:total,:reused,:children_completed)
                haskey(event,key) || continue
                value=something(event[key],-1)
                changed |= get(a,string(key),nothing)!=value
                a[string(key)]=value
            end
            if get(event,:reused,0)>0 || get(event,:partial_recovery,false)
                a["rate_count"]=count-get(event,:reused,0); a["rate_time"]=now
                a["rate_eligible"]=false
            end
            if get(operand,"total",-1)<0 && get(event,:total,nothing)!==nothing && get(event,:points,1)==1
                operand["total"]=event.total
            end
            if haskey(event,:batch)
                a["batch"]=event.batch
            end
            if get(a,"interval_recovered",false) && delta>0
                a["rate_count"]=count-get(a,"reused",0); a["rate_time"]=now
                a["rate_eligible"]=true; a["interval_recovered"]=false
            else
                _scan_rate!(tracker,now;final=count==get(a,"total",-1))
            end
        end
        # Heartbeats/contact are not work observations or duration evidence.
        get(event,:recovered,false) === true && (a["child_recovered"]=true)
        if get(event,:partial_recovery,false)
            a["rate_eligible"]=false
            a["interval_recovered"]=true
            pop!(a,"scan_seconds",nothing)
            role=get(a,"role","")
            if role in ("reference","candidate")
                pop!(row[role],"scan_seconds",nothing)
                row[role]["opaque_interval"]=true
                row[role]["recovered_execution"]=true
            end
        end
        if haskey(event,:capacity) && get(a,"capacity",event.capacity)!=event.capacity
            a["rate_count"]=get(a,"completed",0)-get(a,"reused",0)
            a["rate_time"]=now
            pop!(a,"scan_seconds",nothing)
            role=get(a,"role","")
            if role in ("reference","candidate")
                pop!(row[role],"scan_seconds",nothing)
                row[role]["opaque_interval"]=true
            end
        end
        haskey(event,:capacity) && (a["capacity"]=event.capacity)
        if haskey(event,:stage) && string(event.stage)!=get(a,"stage","")
            a["stage"]=string(event.stage)
            for key in ("rejected","workers","frequencies_completed","frequencies_total")
                pop!(a,key,nothing)
            end
            changed=true
        end
        for key in (:rejected,:frequencies_completed,:frequencies_total,:workers,:point,:formulation)
            haskey(event,key) || continue
            value=event[key]
            changed |= get(a,string(key),nothing)!=value
            a[string(key)]=value
        end
        changed && (tracker.work_time=now)
        return true
    end
    accepted === true && _publish_progress!(tracker;force)
    return nothing
end

# Called only after publication eligibility, never from spinner frames.
function _progress_budget(tracker,now)
    active_budget=-1.0; active_age=0.0; queued=0.0; provisional=false; unknown=false
    whole_seed=get(tracker.evidence,("","","benchmark"),-1.0)
    forecasts=Tuple{Dict{String,Any},Float64,Float64,Float64,Bool}[]
    for row in tracker.rows
        row["state"] in _PROGRESS_TERMINAL_STATES && continue
        if !get(row,"resolved",false)
            push!(forecasts,(row,-1.0,-1.0,0.0,true)); continue
        end
        remaining=0.0; whole=0.0; residual=-1.0; age=0.0; rough=false; known=true
        for role in ("reference","candidate")
            op=row[role]
            estimate=get(op,"estimate",-1.0)
            fallback=false
            if estimate < 0
                estimate,fallback=_duration_estimate(tracker,op,"operation")
            end
            scan=get(op,"scan_seconds",-1.0)
            if scan < 0 && get(op,"batch",1)==1
                scan,scan_rough=_duration_estimate(tracker,op,"scan")
                fallback |= scan_rough
            end
            total=get(op,"total",-1)
            get(op,"opaque_interval",false) && (scan=-1.0)
            # A batch's effective scan rate stays within that very operation.
            if scan >= 0 && total >= 0 && !get(op,"opaque_interval",false)
                estimate=scan*total
                fallback=true # Unmeasured operation preparation/persistence remains an allowance.
            end
            if estimate < 0
                known=false; continue
            end
            whole+=estimate
            rough |= fallback
            if op["state"] in ("pending","running")
                budget=scan>=0 && total>=0 ? scan*max(0,total-get(op,"completed",0)) : estimate
                if op["state"]=="running"
                    anchor=scan>=0 && total>=0 ? get(op,"count_updated",get(op,"started",now)) : get(op,"started",now)
                    residual=budget; age=max(0.0,now-anchor)
                else
                    remaining+=budget
                end
            end
            n=get(op,"performance_left",0)
            warm=get(op,"warmup_left",false)
            if n > 0 || warm
                cost,perf_rough=_duration_estimate(tracker,op,"controlled")
                if cost < 0
                    cost=estimate; perf_rough=true # Explicit whole-call policy assumption, never a sample.
                end
                limit=max(0.0,get(row,"performance_limit",0.0)-get(op,"performance_elapsed",0.0))
                possible=limit/max(cost,eps())
                calls=possible>=n ? n : max(1,ceil(Int,possible))
                performance_budget=(calls+Int(warm))*cost
                if tracker.paused && get(tracker.active,"benchmark","")==row["id"] && get(tracker.active,"role","")==role
                    residual=cost; age=max(0.0,now-get(tracker.active,"sample_started",now))
                    performance_budget=max(0.0,performance_budget-cost)
                end
                remaining+=performance_budget
                whole+=(calls+Int(warm))*cost
                rough |= perf_rough
            end
        end
        if known
            # Establish once from the whole forecast; never shrink with its countdown.
            overhead=get!(row,"overhead",max(1.0,0.1whole))
            whole+=overhead
            if row["state"]=="running" && all(row[r]["state"] in _PROGRESS_TERMINAL_STATES for r in ("reference","candidate"))
                if !tracker.paused
                    residual=overhead; age=max(0.0,now-get(row,"finalization_started",now))
                else
                    remaining+=overhead
                end
            else
                remaining+=overhead
            end
            rough=true # The uncovered-work allowance is explicitly provisional.
            whole_seed < 0 && (whole_seed=whole)
        else
            whole=get(row,"whole_estimate",-1.0)
            remaining=whole; residual=-1.0
            if row["state"]=="running" && whole>=0
                residual=whole; remaining=0.0; age=max(0.0,now-get(row,"started",now))
            end
        end
        push!(forecasts,(row,known ? remaining : whole,residual,age,rough || !known))
    end
    for (row,remaining,residual,age,rough) in forecasts
        if remaining < 0
            remaining=whole_seed
            rough=true
            if row["state"]=="running" && remaining>=0
                residual=remaining; remaining=0.0; age=max(0.0,now-get(row,"started",now))
            end
        end
        remaining < 0 && (unknown=true; continue)
        queued+=remaining
        residual>=0 && (active_budget=residual; active_age=age)
        provisional |= rough
    end
    return Dict("known"=>!unknown,"active_seconds"=>active_budget,
        "anchor_age_seconds"=>active_age,"queued_seconds"=>queued,"provisional"=>provisional)
end

"""Evaluate cached remaining wall-time budgets in seconds without learning or IO."""
function _remaining_budget(snapshot,advance=0.0)
    get(snapshot,"closed",false) && return (seconds=0.0,provisional=false,revising=false)
    budget=get(snapshot,"eta",Dict())
    get(budget,"known",false) || return (seconds=-1.0,provisional=false,revising=false)
    b=get(budget,"active_seconds",-1.0)
    a=max(0.0,get(budget,"anchor_age_seconds",0.0)+advance)
    residual=b<0 ? 0.0 : max(b-a,1.0+0.25max(a-b,0.0))
    seconds=residual+get(budget,"queued_seconds",0.0)
    return (seconds=_valid_seconds(seconds) ? seconds : -1.0,
        provisional=get(budget,"provisional",false) || b>=0 && a>b,
        revising=b>=0 && a>b)
end

function _progress_snapshot(tracker)
    lock(tracker.state_lock) do
        now=tracker.closed ? tracker.stopped : tracker.clock()
        rows=[Dict(k=>v for (k,v) in row if k in ("id","case","state","attempt","execution_state")) for row in tracker.rows]
        counts=Dict(state=>count(row->row["state"]==state,tracker.rows)
            for state in ("pending","running",_PROGRESS_TERMINAL_STATES...))
        keys=("benchmark","attempt","role","backend","mode","stage","completed","total","reused",
            "rejected","workers","frequencies_completed","frequencies_total","point","formulation","sample","samples","sample_outcome")
        active=Dict(k=>v for (k,v) in tracker.active if k in keys)
        return Dict{String,Any}("schema"=>2,"session"=>tracker.session,"revision"=>tracker.revision,
            "selected"=>tracker.selected,"updated_unix_seconds"=>tracker.wall_clock(),
            "elapsed_seconds"=>max(0.0,now-tracker.started),"work_age_seconds"=>max(0.0,now-tracker.work_time),
            "closed"=>tracker.closed,"termination"=>tracker.termination,"measurement_active"=>tracker.paused,
            "benchmarks"=>rows,"counts"=>counts,"active"=>active,"eta"=>_progress_budget(tracker,now))
    end
end

function _write_progress(path,snapshot)
    mkpath(dirname(path))
    temporary=tempname(dirname(path))
    try
        open(temporary,"w") do io
            TOML.print(io,snapshot;sorted=true)
        end
        Base.rename(temporary,path)
    finally
        isfile(temporary) && rm(temporary)
    end
end

function _progress_diagnostic(message,error)
    error isa InterruptException && throw(error)
    try
        @warn message exception=error
    catch diagnostic
        diagnostic isa InterruptException && rethrow()
    end
end

function _publish_progress!(tracker; force=false, boundary=false)
    tracker.mode===:off && return
    lock(tracker.publication_lock) do
        now=tracker.clock()
        tracker.paused && !boundary && return
        due=force || now-tracker.last_snapshot>=1.0
        due || return
        tracker.snapshot_failed && (tracker.output_failed || tracker.mode!==:plain && !tracker.closed) && return
        tracker.last_snapshot=now
        tracker.revision+=1
        snapshot=_progress_snapshot(tracker)
        if !tracker.snapshot_failed
            try
                _write_progress(joinpath(tracker.root,"sessions",tracker.session*".progress.toml"),snapshot)
            catch error
                tracker.snapshot_failed=true
                _progress_diagnostic("Campaign progress snapshot disabled",error)
            end
        end
        if !tracker.output_failed && (tracker.closed || tracker.mode===:plain && (now-tracker.last_plain>=10.0))
            try
                counts=snapshot["counts"]
                finished=sum(counts[k] for k in _PROGRESS_TERMINAL_STATES)
                println(tracker.io,"Benchmarks $finished/$(length(tracker.rows)) finished | Complete $(counts["complete"]) | Failed $(counts["failed"]) | ",
                    tracker.closed ? "ETA "*(tracker.termination=="exhausted" ? "done" : "stopped") : get(snapshot["active"],"stage","preparing"))
                flush(tracker.io)
                tracker.last_plain=now
            catch error
                tracker.output_failed=true
                _progress_diagnostic("Campaign progress output disabled",error)
            end
        end
    end
    return nothing
end

_shell_quote(value) = "'"*replace(string(value),"'"=>"'\"'\"'")*"'"

function _with_campaign_progress(f,directory,ids,session;progress=:auto,selected=false,io=stderr)
    _progress_mode(progress)
    existing=_CAMPAIGN_TRACKER[]
    existing===nothing || return f(existing)
    progress===:off && return f(nothing)
    tracker=CampaignProgress(directory,ids,session;mode=progress,selected,io)
    _publish_progress!(tracker;force=true)
    try
        println(io,"Watch in another terminal: ./gauntlet/lcm gauntlet status --directory ",
            _shell_quote(tracker.root)," --watch --session ",_shell_quote(tracker.session))
    catch error
        tracker.output_failed=true
        _progress_diagnostic("Campaign progress output disabled",error)
    end
    return Base.ScopedValues.with(_CAMPAIGN_TRACKER=>tracker) do
        try
            value=LineCableModels.with_progress(event->_progress_event!(tracker,event);scope=(session=tracker.session,)) do
                f(tracker)
            end
            tracker.termination="exhausted"
            return value
        catch error
            tracker.termination="aborted"
            # The invocation owner has aborted its selected traversal.
            lock(tracker.state_lock) do
                for row in tracker.rows
                    row["state"]=="running" && (row["state"]=error isa InterruptException ? "interrupted" : "failed")
                    row["state"]=="pending" && (row["state"]="skipped")
                    for role in ("reference","candidate")
                        row[role]["state"]=="running" && (row[role]["state"]=row["state"])
                        row[role]["state"]=="pending" && (row[role]["state"]="skipped")
                    end
                end
            end
            rethrow()
        finally
            lock(tracker.state_lock) do
                tracker.closed=true; tracker.paused=false
                tracker.stopped=tracker.clock(); tracker.work_time=tracker.stopped
            end
            _publish_progress!(tracker;force=true)
        end
    end
end

function _performance_span(f;sample,samples,role,record=(_->nothing))
    tracker=_CAMPAIGN_TRACKER[]
    if tracker===nothing
        value=LineCableModels.with_performance_sample(f)
        record(value)
        return value
    end
    saved=Dict{String,Any}()
    lock(tracker.publication_lock) do
        lock(tracker.state_lock) do
            saved=copy(tracker.active)
            merge!(tracker.active,Dict("stage"=>"performance","role"=>string(role),
                "sample"=>sample,"samples"=>samples,"sample_started"=>tracker.clock()))
            tracker.paused=true; tracker.work_time=tracker.clock()
        end
        _publish_progress!(tracker;force=true,boundary=true)
    end
    outcome="complete"
    try
        value=LineCableModels.with_performance_sample(f)
        record(value)
        return value
    catch error
        outcome=error isa InterruptException ? "interrupted" : "failed"
        rethrow()
    finally
        lock(tracker.publication_lock) do
            lock(tracker.state_lock) do
                elapsed=max(0.0,tracker.clock()-get(tracker.active,"sample_started",tracker.clock()))
                id=get(saved,"benchmark","")
                if haskey(tracker.indices,id)
                    row=tracker.rows[tracker.indices[id]]
                    haskey(row,"finalization_started") && (row["finalization_started"]+=elapsed)
                end
                tracker.active=saved
                tracker.active["sample_outcome"]=outcome
                tracker.paused=false; tracker.work_time=tracker.clock()
            end
            _publish_progress!(tracker;force=true,boundary=true)
        end
    end
end

function _performance_observation!(role;seconds=0.0,sample=0,finished=false,reused=false)
    tracker=_CAMPAIGN_TRACKER[]
    tracker===nothing && return
    lock(tracker.state_lock) do
        id=get(tracker.active,"benchmark","")
        haskey(tracker.indices,id) || return
        row=tracker.rows[tracker.indices[id]]; op=row[string(role)]
        if finished
            op["performance_left"]=0; op["warmup_left"]=false
        elseif sample=="warmup"
            op["warmup_left"]=false
        else
            op["performance_left"]=max(0,get(op,"performance_left",0)-1)
            op["performance_elapsed"]=get(op,"performance_elapsed",0.0)+seconds
            reused || _duration_evidence!(tracker,(op["backend"],op["mode"],"controlled"),seconds)
        end
    end
end
