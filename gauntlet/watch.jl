# The watcher owns terminal rendering. It never restores declarations or results.
function _safe_session(session)
    value=string(session)
    occursin(r"^[A-Za-z0-9_-][A-Za-z0-9_.-]*$",value) ||
        throw(ArgumentError("session must be one safe filename component"))
    return value
end

_clean_text(value) = join(iscntrl(c) ? ' ' : c for c in string(value))
function _fit_text(value,width)
    width=max(0,width); result=IOBuffer(); columns=0
    for c in _clean_text(value)
        next=columns+textwidth(c)
        next>width && break
        print(result,c); columns=next
    end
    return String(take!(result))
end

function _progress_duration(seconds;remaining=false)
    _valid_seconds(seconds) || return "--"
    value=remaining ? (seconds<60 ? ceil(seconds) : seconds<3600 ? 10ceil(seconds/10) : 60ceil(seconds/60)) : floor(seconds)
    value <= typemax(Int) || return "--"
    n=Int(value)
    return lpad(div(n,3600),2,'0')*":"*lpad(div(n%3600,60),2,'0')*":"*lpad(n%60,2,'0')
end

function _watch_lines(snapshot,width;advance=0.0,spinner='|',notice="")
    counts=snapshot["counts"]; active=get(snapshot,"active",Dict())
    total=length(snapshot["benchmarks"])
    finished=sum(get(counts,s,0) for s in _PROGRESS_TERMINAL_STATES)
    label=get(snapshot,"inventory",false) ? "Inventory" : get(snapshot,"selected",false) ? "Selected benchmarks" : "Benchmarks"
    suffix=" $finished/$total finished"
    bars=clamp(width-textwidth(label)-textwidth(suffix)-3,0,20)
    filled=total==0 ? 0 : clamp(floor(Int,bars*finished/total),0,bars)
    first=label*" ["*repeat("=",filled)*repeat(".",bars-filled)*"]"*suffix
    fields=["Complete $(get(counts,"complete",0))","Failed $(get(counts,"failed",0))",
        "Running $(get(counts,"running",0))","Pending $(get(counts,"pending",0))"]
    for state in ("skipped","interrupted")
        get(counts,state,0)>0 && push!(fields,uppercasefirst(state)*" $(counts[state])")
    end
    second=join(fields," | ")
    if textwidth(second)>width
        second=join(["$(label) $(get(counts,state,0))" for (label,state) in
            (("Done","complete"),("F","failed"),("Run","running"),("Pend","pending"),("Skip","skipped"),("Int","interrupted"))
            if state in ("complete","failed","running","pending") || get(counts,state,0)>0]," | ")
    end
    id=get(active,"benchmark","")
    index=findfirst(row->get(row,"id","")==id,snapshot["benchmarks"])
    case=index===nothing ? "No active benchmark" : get(snapshot["benchmarks"][index],"case",id)
    identity="$(something(index,"--"))/$total "
    mode=get(active,"mode","")
    tail=" | "*join(filter(!isempty,[get(active,"role",""),get(active,"backend",""),mode=="ordinary" ? "" : mode])," / ")
    third=identity*_fit_text(case,width-textwidth(identity)-textwidth(tail))*tail
    closed=get(snapshot,"closed",false)
    elapsed=get(snapshot,"elapsed_seconds",-1.0)+(closed ? 0.0 : advance)
    eta=_remaining_budget(snapshot,advance)
    termination=get(snapshot,"termination","")
    eta_text=closed ? (termination=="exhausted" ? "done" : termination=="aborted" ? "stopped" : "--") :
        eta.seconds<0 ? "--" : "~"*_progress_duration(eta.seconds;remaining=true)
    fourth="Elapsed $(_progress_duration(elapsed)) | ETA $eta_text"
    qualifier=eta.revising ? " (revising)" : eta.provisional ? " (provisional)" : ""
    textwidth(fourth*qualifier)<=width && (fourth*=qualifier)
    frame=closed ? ' ' : spinner
    if get(snapshot,"measurement_active",false)
        sample=get(active,"sample","preparation")
        detail=sample=="warmup" ? "warmup; observation suspended" :
            "performance sample $sample/$(get(active,"samples","?")); observation suspended"
        fifth="Current run $frame | $detail"
    else
        done=get(active,"completed",-1); target=get(active,"total",-1)
        fifth="Current run $frame ($(done<0 ? "--" : done)/$(target<0 ? "?" : target) scans)"
        details=String[replace(get(active,"stage","unavailable"),'_'=>' ')]
        haskey(active,"sample_outcome") && push!(details,"sample "*active["sample_outcome"])
        if haskey(active,"point")
            push!(details,"point $(active["point"])")
        elseif haskey(active,"formulation")
            push!(details,"formulation $(active["formulation"])")
        end
        if haskey(active,"frequencies_completed")
            push!(details,"frequencies $(active["frequencies_completed"])/$(get(active,"frequencies_total","?"))")
        elseif get(active,"rejected",0)>0
            push!(details,"rejected $(active["rejected"])")
        elseif get(active,"workers",0)>0
            push!(details,"workers $(active["workers"])")
        end
        fifth*=" | "*join(details[1:min(end,3)],"; ")
    end
    age=max(0.0,get(snapshot,"work_age_seconds",0.0)+(closed ? 0.0 : advance))
    sixth=!isempty(notice) ? notice : closed ? "Invocation "*get(snapshot,"termination","closed") :
        "Last work observation $(floor(Int,age))s ago"*(age>5 ? " (stale; liveness unknown)" : "")
    return [_fit_text(line,width) for line in (first,second,third,fourth,fifth,sixth)]
end

"""Cached observation and viewer-local clock anchors; contains no execution objects."""
Base.@kwdef mutable struct CampaignWatch
    "Explicit or automatically pinned invocation."
    session::Union{Nothing,String} = nothing
    "Last complete valid snapshot."
    snapshot::Dict{String,Any} = Dict{String,Any}()
    "Viewer monotonic attachment to the cached revision in seconds."
    received::Float64 = 0.0
    "Snapshot wall age at attachment in seconds."
    age::Float64 = 0.0
    "Latest poll's viewer monotonic time in seconds."
    polled::Float64 = -Inf
    "Latest poll's wall time in seconds."
    wall::Float64 = NaN
    "Compact availability or identity notice."
    notice::String = ""
    "Undecorated terminal render baseline."
    lines::Vector{String} = String[]
    "Number of terminal rows owned since the last resize."
    painted::Int = 0
    "Last terminal dimensions."
    dimensions::Tuple{Int,Int} = (0,0)
    "Last plain summary time in seconds."
    plain_time::Float64 = -Inf
    "Last plain summary identity."
    plain_key::String = ""
end

function _watch_inventory(root)
    records=try
        campaign_status(root;verify=false)
    catch error
        error isa InterruptException && rethrow()
        NamedTuple[]
    end
    rows=[Dict{String,Any}("id"=>row.id,"state"=>string(row.state)) for row in records]
    counts=Dict(s=>count(row->row["state"]==s,rows) for s in ("pending","running",_PROGRESS_TERMINAL_STATES...))
    snapshot=Dict{String,Any}("schema"=>2,"inventory"=>true,"benchmarks"=>rows,"counts"=>counts,
        "elapsed_seconds"=>-1.0,"active"=>Dict{String,Any}(),"eta"=>Dict("known"=>false))
    return snapshot,records
end

# Validate disposable input before it reaches arithmetic or terminal rendering.
function _valid_watch_snapshot(snapshot,session)
    get(snapshot,"session",nothing)==session && get(snapshot,"schema",0) in (1,2) || return false
    rows=get(snapshot,"benchmarks",nothing); counts=get(snapshot,"counts",nothing)
    active=get(snapshot,"active",nothing)
    rows isa Vector && counts isa AbstractDict && active isa AbstractDict || return false
    all(row->row isa AbstractDict && get(row,"id",nothing) isa String &&
        occursin(r"^[A-Za-z0-9_-][A-Za-z0-9_.-]*$",row["id"]) &&
        get(row,"case","") isa String,rows) || return false
    all(x->x isa Integer && 0<=x<=length(rows),values(counts)) || return false
    bounded(x,low=0) = x isa Real && !(x isa Bool) && isfinite(x) && low<=x<=1e12
    bounded(get(snapshot,"elapsed_seconds",nothing)) || return false
    all(key->get(snapshot,key,false) isa Bool,("closed","selected","measurement_active")) || return false
    get(snapshot,"termination","") isa String || return false
    if get(snapshot,"schema",0)==1
        return bounded(get(snapshot,"updated_unix_seconds",0.0))
    end
    get(snapshot,"revision",nothing) isa Integer && snapshot["revision"]>=0 || return false
    bounded(get(snapshot,"work_age_seconds",nothing)) || return false
    bounded(get(snapshot,"updated_unix_seconds",nothing)) || return false
    all(key->get(active,key,"") isa String,
        ("benchmark","attempt","role","backend","mode","stage","sample_outcome")) || return false
    all(key->get(active,key,0) isa Integer && -1<=get(active,key,0)<=1e12,
        ("completed","total","reused","rejected","workers","frequencies_completed",
            "frequencies_total","point","formulation","samples")) || return false
    get(active,"sample",0) isa Union{String,Integer} || return false
    eta=get(snapshot,"eta",nothing)
    eta isa AbstractDict && get(eta,"known",false) isa Bool && get(eta,"provisional",false) isa Bool || return false
    return bounded(get(eta,"active_seconds",-1.0),-1) &&
        bounded(get(eta,"anchor_age_seconds",0.0)) && bounded(get(eta,"queued_seconds",0.0))
end

function _watch_poll!(watch,root,now,wall)
    previous_poll=watch.polled
    adjusted=isfinite(watch.wall) && abs((wall-watch.wall)-(now-previous_poll))>2.0
    watch.polled=now; watch.wall=wall
    if watch.session===nothing
        inventory,records=_watch_inventory(root)
        sessions=Set{String}()
        for row in records
            row.state===:running || continue
            state=_read_progress_toml(joinpath(root,row.id,"state.toml"))
            session=get(state,"session",nothing)
            session isa String && occursin(r"^[A-Za-z0-9_-][A-Za-z0-9_.-]*$",session) && push!(sessions,session)
        end
        if length(sessions)==1
            watch.session=only(sessions)
        else
            watch.snapshot=inventory; watch.received=now; watch.age=0.0
            watch.notice=length(sessions)>1 ? "Inventory only: multiple sessions; use --session" : "Inventory only: no unambiguous active session"
            return
        end
    end
    path=joinpath(root,"sessions",watch.session*".progress.toml")
    snapshot=_read_progress_toml(path)
    valid=_valid_watch_snapshot(snapshot,watch.session)
    if !valid
        watch.notice=isfile(path) ? "Snapshot malformed or unsupported; retaining last observation" : "Snapshot unavailable; waiting for selected session"
        if isempty(watch.snapshot)
            watch.snapshot=first(_watch_inventory(root)); watch.received=now
        end
        return
    end
    # Older snapshots folded comparison verdicts into job status. Use their
    # recorded execution state for display, without rewriting retained files.
    for row in snapshot["benchmarks"]
        get(row,"state",nothing)=="failed" && get(row,"execution_state",nothing)=="complete" &&
            get(row,"verdict_failed",false)===true || continue
        row["state"]="complete"
    end
    snapshot["counts"]=Dict(state=>count(row->get(row,"state",nothing)==state,snapshot["benchmarks"])
        for state in ("pending","running",_PROGRESS_TERMINAL_STATES...))
    active=snapshot["active"]
    id=get(active,"benchmark",nothing)
    if id!==nothing && !get(snapshot,"closed",false)
        id isa String && any(row->row["id"]==id,snapshot["benchmarks"]) || return
        state=_read_progress_toml(joinpath(root,id,"state.toml"))
        if !isempty(state) && (get(state,"session",watch.session)!=watch.session ||
                get(state,"attempt",nothing)!=get(active,"attempt",nothing))
            watch.notice="Attempt ownership changed; retaining selected invocation"
            return
        end
    end
    if get(snapshot,"schema",0)==1
        # Never reinterpret historical jobs/trials as the new scan unit.
        snapshot["active"]=Dict{String,Any}("stage"=>"legacy snapshot")
        snapshot["eta"]=Dict("known"=>false)
        watch.notice="Legacy snapshot: scan accounting unavailable"
    else
        watch.notice=adjusted ? "Wall clock adjusted; snapshot age approximate" : ""
    end
    if get(watch.snapshot,"session",nothing)==watch.session &&
            get(snapshot,"revision",0)<=get(watch.snapshot,"revision",-1)
        return
    end
    watch.snapshot=snapshot
    watch.received=now
    age=wall-get(snapshot,"updated_unix_seconds",wall)
    watch.age=adjusted || !_valid_seconds(age) ? 0.0 : age
end

_watch_terminal(io) = io isa IOContext ? _watch_terminal(io.io) : io isa Base.TTY
function _watch_color(io)
    return _watch_terminal(io) && get(io,:color,true) && !haskey(ENV,"NO_COLOR") && get(ENV,"TERM","")!="dumb"
end

function _watch_styled(line,index,color)
    color || return line
    if index==2
        return replace(line,r"(?:Complete|Done|Failed|F|Interrupted|Int) [0-9]+"=>
            value->(startswith(value,"Complete") || startswith(value,"Done") ? "\e[32m" : "\e[31m")*value*"\e[0m")
    end
    style=index in (3,5) ? "\e[36m" : index==6 ? "\e[2m" :
        index==4 ? (occursin("done",line) ? "\e[32m" : occursin("stopped",line) ? "\e[31m" : "\e[33m") : ""
    return style*line*"\e[0m"
end

function _watch_render!(watch,io,lines,dimensions;terminal=false,color=false,now=0.0)
    height,width=dimensions
    fits=terminal && height>=7 && width>=62
    buffer=IOBuffer()
    if watch.dimensions!=dimensions
        # Explicit line endings and one spare column keep the panel from wrapping.
        # After a height reduction, touch only owned rows still visible above us.
        rows=min(watch.painted,max(0,height-1))
        if terminal && rows>0
            print(buffer,"\e[$(rows)A")
            for _ in 1:rows
                print(buffer,"\r\e[2K\r\n")
            end
            print(buffer,"\e[$(rows)A")
        end
        empty!(watch.lines); watch.painted=0
    end
    if fits
        if isempty(watch.lines)
            print(buffer,"\e[?25l")
            for (i,line) in enumerate(lines)
                print(buffer,_watch_styled(line,i,color),"\r\n")
            end
            watch.painted=6
        else
            for i in eachindex(lines)
                lines[i]==watch.lines[i] && continue
                distance=7-i
                if i==5 && length(lines[i])>=13 && length(watch.lines[i])>=13 &&
                        lines[i][1:12]==watch.lines[i][1:12] && lines[i][14:end]==watch.lines[i][14:end]
                    print(buffer,"\e[2A\r\e[12C",_watch_styled(string(lines[i][13]),5,color),"\r\e[2B")
                else
                    print(buffer,"\e[$(distance)A\r\e[2K",_watch_styled(lines[i],i,color),"\r\e[$(distance)B")
                end
            end
        end
        watch.lines=copy(lines)
    else
        # Plain output contains neither spinner frames nor ticking-only reports.
        current=replace(lines[5],r"^Current run ."=>"Current run")
        closed=get(watch.snapshot,"closed",false)
        key=join((lines[1],lines[2],lines[3],current,string(closed))," | ")
        if isempty(watch.plain_key) || key!=watch.plain_key && (closed || now-watch.plain_time>=10.0)
            summary=join((lines[1],lines[2],lines[4],current)," | ")
            print(buffer,terminal ? _fit_text(summary,max(0,width-1)) : summary,"\r\n")
            watch.plain_key=key; watch.plain_time=now
            terminal && (watch.painted+=1)
        end
    end
    watch.dimensions=dimensions
    bytes=take!(buffer)
    if !isempty(bytes)
        write(io,bytes); flush(io)
    end
    return nothing
end

"""
    watch_campaign(directory; io=stderr, interval=1.0, session=nothing)

Watch lightweight snapshots in a separate process until interrupted. `session`
pins an invocation, including its final state; omitted selection uses conservative
metadata-only discovery. `interval` is the snapshot polling interval in seconds.
The terminal spinner advances every 0.25 seconds from cached state. Redirected
output receives throttled plain summaries. Closing this viewer never affects
execution. No declarations, numerical results, or worker artifacts are loaded.
"""
function watch_campaign(directory;io=stderr,interval=1.0,session=nothing)
    isfinite(interval) && interval>0 || throw(ArgumentError("watch interval must be positive"))
    watch=CampaignWatch(;session=session===nothing ? nothing : _safe_session(session))
    root=abspath(directory)
    terminal=_watch_terminal(io) && get(ENV,"TERM","")!="dumb"
    frame=0; next_text=-Inf; lines=String[]
    try
        while true
            now=_progress_clock()
            if now-watch.polled>=interval
                _watch_poll!(watch,root,now,time())
            end
            dimensions=displaysize(io)
            if now>=next_text || isempty(lines) || watch.dimensions!=dimensions
                advance=get(watch.snapshot,"closed",false) ? 0.0 : watch.age+max(0.0,now-watch.received)
                # Missing inventory elapsed must stay unavailable.
                get(watch.snapshot,"inventory",false) && (advance=0.0)
                lines=_watch_lines(watch.snapshot,max(0,dimensions[2]-1);advance,notice=watch.notice)
                next_text=now+1.0
            end
            if terminal && !get(watch.snapshot,"closed",false)
                spinner=('|','/','-','\\')[frame%4+1]
                if startswith(lines[5],"Current run ") && ncodeunits(lines[5])>=13
                    lines[5]=lines[5][1:12]*spinner*lines[5][14:end]
                end
                frame+=1
            end
            _watch_render!(watch,io,lines,dimensions;terminal,color=_watch_color(io),now)
            sleep(0.25)
        end
    catch error
        error isa Union{InterruptException,Base.IOError} || rethrow()
    finally
        if terminal
            try
                print(io,"\e[0m\e[?25h"); flush(io)
            catch
                # A closed output pipe cannot affect the calculation process.
            end
        end
    end
    return nothing
end

export watch_campaign
