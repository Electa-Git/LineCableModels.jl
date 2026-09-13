@testitem "Gauntlet / progress / scan ownership, ordering and quiet scopes" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels, TOML, Logging
    const LCM=LineCableModels
    mktempdir() do root
        clock=Ref(0.0)
        tracker=Gauntlet.CampaignProgress(root,[:a,:b],"session";clock=()->clock[],io=IOBuffer())
        emit(event)=Gauntlet._progress_event!(tracker,event)
        LCM.with_progress(emit;scope=(session="session",benchmark=:a,attempt="one")) do
            receiver=LCM.progress_receiver()
            LCM.report_progress(receiver,(kind=:benchmark,state=:running))
            LCM.report_progress(receiver,(kind=:operand,role=:reference,state=:running,backend="Opaque",mode="Monte Carlo"))
            LCM.with_progress_scope(role=:reference) do
                old=Ref{Any}()
                LCM.with_scan_progress(;total=nothing) do scan
                    old[]=scan
                    LCM.report_progress(scan,(kind=:scan,completed=0,total=nothing,stage=:sampling))
                    @test tracker.active["total"] == -1
                    LCM.report_progress(scan,(frequencies_completed=18,frequencies_total=48,workers=3,stage=:solving))
                    @test tracker.active["completed"] == 0
                    @test tracker.active["total"] == -1
                    clock[]=4.0
                    LCM.report_progress(scan,(kind=:scan,completed=0,total=1000,children_completed=2,rejected=2,stage=:sampling))
                    @test !haskey(tracker.active,"workers")
                    clock[]=10.0
                    LCM.report_progress(scan,(kind=:scan,completed=1,total=1000,children_completed=3))
                    @test tracker.active["scan_seconds"] ≈ 10.0 # Includes rejection time.
                    before=tracker.work_time
                    clock[]=12.0
                    LCM.report_progress(scan,(heartbeat_unix_seconds=10.0,))
                    @test tracker.work_time==before
                    @test Gauntlet._progress_snapshot(tracker)["work_age_seconds"]==2.0
                    LCM.with_scan_progress(;total=48) do child
                        LCM.report_progress(child,(kind=:scan,completed=48,total=48))
                    end
                    @test tracker.active["completed"]==1
                    @test tracker.active["total"]==1000
                    for index in 2:1000
                        clock[]+=0.001
                        LCM.report_progress(scan,(kind=:scan,completed=index,total=1000))
                    end
                    @test length(tracker.scopes)==1
                    @test length(tracker.evidence)<=2
                end
                @test isempty(tracker.scopes)
                previous=copy(tracker.active)
                LCM.report_progress(old[],(stage=:obsolete,))
                @test tracker.active==previous
                sequence=tracker.sequence
                revision=tracker.revision
                emit((session="session",benchmark=:a,attempt="one",role=:reference,sequence,stage=:late))
                emit((session="old",benchmark=:a,attempt="one",role=:reference,sequence=sequence+100,stage=:late))
                @test tracker.active==previous
                @test tracker.revision==revision
                emit((kind=:benchmark,session="session",benchmark=:a,attempt="superseded",state=:complete,sequence=sequence+200))
                @test tracker.active==previous
                @test tracker.rows[1]["state"]=="running"
            end
        end
        @test LCM.progress_receiver()===nothing
        events=NamedTuple[]
        LCM.with_progress(e->push!(events,e)) do
            @test_throws ErrorException LCM.with_performance_sample() do
                @test LCM.progress_receiver()===nothing
                LCM.with_scan_progress(;total=1) do receiver
                    @test receiver===nothing
                end
                error("sample failed")
            end
            @test LCM.progress_receiver()!==nothing
        end
        @test isempty(events)
        @test !LCM.performance_sample_active()
        @test_throws ArgumentError Gauntlet._progress_mode(:bad)
    end
end

@testitem "Gauntlet / progress / exhaustion, verdicts and suspended publication" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels, TOML
    mktempdir() do root
        tracker_ref=Ref{Any}()
        Gauntlet._with_campaign_progress(root,[:a,:b],"closed";io=IOBuffer()) do tracker
            tracker_ref[]=tracker
            LineCableModels.with_progress_scope(benchmark=:a,attempt="one") do
                receiver=LineCableModels.progress_receiver()
                LineCableModels.report_progress(receiver,(kind=:benchmark,state=:running))
                LineCableModels.report_progress(receiver,(kind=:operand,role=:reference,state=:running,backend="new",mode="ordinary"))
                LineCableModels.with_progress_scope(role=:reference) do
                    receiver=LineCableModels.progress_receiver()
                    LineCableModels.report_progress(receiver,(kind=:scan_result,completed=1,total=1,reused=0))
                end
                LineCableModels.report_progress(receiver,(kind=:operand,role=:reference,state=:complete,seconds=3.0))
                path=joinpath(root,"sessions","closed.progress.toml")
                recorded=Ref(false)
                Gauntlet._performance_span(;sample=2,samples=3,role=:candidate,record=value->(recorded[]=true)) do
                    before=read(path)
                    snapshot=TOML.parse(String(copy(before)))
                    @test snapshot["measurement_active"]
                    @test snapshot["active"]["sample"]==2
                    @test LineCableModels.progress_receiver()===nothing
                    @test !recorded[]
                    Gauntlet._publish_progress!(tracker;force=true)
                    @test read(path)==before
                    42
                end
                @test recorded[]
                @test !TOML.parsefile(path)["measurement_active"]
                @test TOML.parsefile(path)["active"]["sample_outcome"]=="complete"
                @test_throws InterruptException Gauntlet._performance_span(;sample=3,samples=3,role=:candidate) do
                    throw(InterruptException())
                end
                @test !tracker.paused
                @test TOML.parsefile(path)["active"]["sample_outcome"]=="interrupted"
                LineCableModels.report_progress(receiver,(kind=:benchmark,state=:complete,verdict_failed=true,seconds=4.0))
                @test tracker.evidence[("new","ordinary","operation")]==3.0
            end
            LineCableModels.report_progress(LineCableModels.progress_receiver(),(kind=:benchmark,benchmark=:b,state=:skipped,attempt="two"))
        end
        snapshot=Gauntlet._progress_snapshot(tracker_ref[])
        @test snapshot["termination"]=="exhausted"
        @test snapshot["counts"]["failed"]==1
        @test snapshot["counts"]["skipped"]==1
        @test occursin("ETA done",Gauntlet._watch_lines(snapshot,120)[4])
        elapsed=snapshot["elapsed_seconds"]
        @test Gauntlet._progress_snapshot(tracker_ref[])["elapsed_seconds"]==elapsed
        @test_throws ErrorException Gauntlet._with_campaign_progress(root,[:c,:d],"abort";io=IOBuffer()) do tracker
            LineCableModels.report_progress(LineCableModels.progress_receiver(),(kind=:benchmark,benchmark=:c,attempt="three",state=:running))
            error("required persistence failed")
        end
        aborted=TOML.parsefile(joinpath(root,"sessions","abort.progress.toml"))
        @test aborted["termination"]=="aborted"
        @test aborted["counts"]["failed"]==1
        @test aborted["counts"]["skipped"]==1
        @test occursin("ETA stopped",Gauntlet._watch_lines(aborted,120)[4])
    end
end

@testitem "Gauntlet / progress / scope-correct budgets and genuine overdue extension" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    snapshot=Dict("closed"=>false,"eta"=>Dict("known"=>true,"active_seconds"=>100.0,
        "anchor_age_seconds"=>0.0,"queued_seconds"=>50.0,"provisional"=>false))
    for (age,expected) in ((70.,30.),(80.,20.),(90.,10.),(99.,1.),(100.,1.),(110.,3.5))
        estimate=Gauntlet._remaining_budget(snapshot,age)
        @test estimate.seconds==expected+50
        @test estimate.revising==(age>100)
    end
    @test Gauntlet._remaining_budget(snapshot,100.01).revising
    @test snapshot["eta"]["anchor_age_seconds"]==0.0
    @test Gauntlet._progress_duration(0.01;remaining=true)=="00:00:01"
    @test Gauntlet._progress_duration(61.;remaining=true)=="00:01:10"
    @test Gauntlet._progress_duration(25*3600.)=="25:00:00"
    @test Gauntlet._progress_duration(NaN)=="--"
    mktempdir() do root
        clock=Ref(0.0)
        tracker=Gauntlet.CampaignProgress(root,[:a,:b],"eta";clock=()->clock[],io=IOBuffer())
        Gauntlet._duration_evidence!(tracker,("owned","Monte Carlo","scan"),0.02)
        @test !Gauntlet._progress_snapshot(tracker)["eta"]["known"]
        # A whole benchmark seed supplies unresolved definitions without converting scan units.
        Gauntlet._duration_evidence!(tracker,("","","benchmark"),200.0)
        budget=Gauntlet._progress_snapshot(tracker)["eta"]
        @test budget["known"] && budget["provisional"]
        @test budget["queued_seconds"]==400.0
        @test length(tracker.evidence)==2
        Gauntlet._duration_evidence!(tracker,("owned","Monte Carlo","scan"),NaN)
        Gauntlet._duration_evidence!(tracker,("owned","Monte Carlo","scan"),Inf)
        @test tracker.evidence[("owned","Monte Carlo","scan")]==0.02
        Gauntlet._duration_evidence!(tracker,("owned","Monte Carlo","scan"),0.1)
        @test tracker.evidence[("owned","Monte Carlo","scan")]≈0.04
        row=tracker.rows[1]; row["resolved"]=true
        for role in ("reference","candidate")
            merge!(row[role],Dict("backend"=>"owned","mode"=>"ordinary","total"=>1,"estimate"=>100.0))
        end
        Gauntlet._progress_snapshot(tracker)
        reserve=row["overhead"]
        row["reference"]["state"]="complete"
        clock[]=90.0
        Gauntlet._progress_snapshot(tracker)
        @test row["overhead"]==reserve
        row["state"]="running"
        merge!(row["candidate"],Dict("state"=>"running","scan_seconds"=>1.0,"total"=>100,
            "completed"=>11,"count_updated"=>11.0,"started"=>0.0))
        live=Gauntlet._progress_budget(tracker,11.5)
        @test live["active_seconds"]==89.0
        @test live["anchor_age_seconds"]==0.5 # Counted scans are not subtracted twice.
        # Batch work cannot seed a scalar scan pool by dividing outputs.
        empty!(tracker.evidence)
        a=tracker.active
        merge!(a,Dict("benchmark"=>"a","role"=>"candidate","rate_eligible"=>true,
            "rate_count"=>0,"rate_time"=>0.0,"completed"=>4,"reused"=>0,"batch"=>4))
        Gauntlet._scan_rate!(tracker,8.0;final=true)
        @test a["scan_seconds"]==2.0
        @test isempty(tracker.evidence)
    end
end

@testitem "Gauntlet / progress / recovery, capacity and required work budgets" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    mktempdir() do root
        clock=Ref(0.0)
        tracker=Gauntlet.CampaignProgress(root,[:a],"recovery";clock=()->clock[],io=IOBuffer())
        emit(e)=Gauntlet._progress_event!(tracker,merge((benchmark=:a,attempt="one",role=:reference),e))
        emit((kind=:benchmark,state=:running))
        emit((kind=:operand,state=:running,backend="FutureBackend",mode="Monte Carlo"))
        emit((kind=:scan_start,scan_scope=1,scan_parent=0,total=4))
        clock[]=2.0
        emit((scan_scope=1,partial_recovery=true,recovered=true,capacity=(2,1)))
        emit((kind=:scan,scan_scope=1,completed=1,total=4))
        @test tracker.active["reused"]==1
        @test isempty(tracker.evidence)
        clock[]=6.0
        emit((kind=:scan,scan_scope=1,completed=2,total=4))
        @test tracker.active["scan_seconds"]==4.0
        @test tracker.evidence[("FutureBackend","Monte Carlo","scan")]==4.0
        evidence=copy(tracker.evidence)
        emit((scan_scope=1,frequencies_completed=18,frequencies_total=48,workers=2))
        @test tracker.evidence==evidence
        emit((scan_scope=1,capacity=(4,1)))
        @test !haskey(tracker.active,"scan_seconds")
        @test tracker.rows[1]["reference"]["opaque_interval"]
        clock[]=8.0
        emit((kind=:scan,scan_scope=1,completed=3,total=4))
        @test tracker.active["scan_seconds"]==2.0
        emit((kind=:operand,state=:complete,seconds=8.0))
        @test !haskey(tracker.evidence,("FutureBackend","Monte Carlo","operation"))
        # A delayed duplicate or failed child cannot stretch an accepted scan's duration.
        duplicate=Gauntlet.CampaignProgress(root,[:b],"duplicate";clock=()->clock[],io=IOBuffer())
        clock[]=0.0
        emit_duplicate(e)=Gauntlet._progress_event!(duplicate,merge((benchmark=:b,attempt="two",role=:reference),e))
        emit_duplicate((kind=:benchmark,state=:running))
        emit_duplicate((kind=:operand,state=:running,backend="FutureBackend",mode="ordinary"))
        emit_duplicate((kind=:scan_start,scan_scope=1,scan_parent=0,total=10))
        clock[]=0.2
        emit_duplicate((kind=:scan,scan_scope=1,completed=1,total=10))
        clock[]=3.0
        emit_duplicate((kind=:scan,scan_scope=1,completed=1,total=10))
        @test isempty(duplicate.evidence)
        emit_duplicate((kind=:scan_end,scan_scope=1,state=:failed))
        @test duplicate.evidence[("FutureBackend","ordinary","scan")]==0.2
        clock[]=8.0
        # One fixed allowance plus queued complete calls; pausing cannot double-count a sample.
        row=tracker.rows[1];row["resolved"]=true
        row["performance_limit"]=100.0
        for (role,cost,calls) in (("reference",10.0,2),("candidate",20.0,0))
            merge!(row[role],Dict("state"=>"complete","estimate"=>cost,"total"=>-1,
                "performance_left"=>calls,"warmup_left"=>false,"performance_elapsed"=>0.0))
        end
        empty!(tracker.evidence)
        ordinary=Gauntlet._progress_budget(tracker,8.0)
        @test ordinary["active_seconds"]==5.0
        @test ordinary["queued_seconds"]==20.0
        tracker.paused=true
        tracker.active["sample_started"]=8.0
        frozen=Gauntlet._progress_budget(tracker,8.0)
        @test frozen["active_seconds"]==10.0
        @test frozen["queued_seconds"]==15.0
        @test row["overhead"]==5.0
        row["performance_limit"]=1e300
        @test Gauntlet._progress_budget(tracker,8.0)["queued_seconds"]==15.0
        @test isempty(tracker.evidence) # Forecast assumptions never become observations.
    end
end

@testitem "Gauntlet / progress / publication failure, atomic files and throttling" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using TOML, Logging
    mktempdir() do root
        clock=Ref(0.0)
        tracker=Gauntlet.CampaignProgress(root,[:a],"publish";clock=()->clock[],io=IOBuffer())
        Gauntlet._publish_progress!(tracker;force=true)
        path=joinpath(root,"sessions","publish.progress.toml")
        before=read(path)
        for _ in 1:1000
            Gauntlet._publish_progress!(tracker)
        end
        @test tracker.revision==1
        @test read(path)==before
        clock[]=1.0
        @sync for _ in 1:10
            Threads.@spawn Gauntlet._publish_progress!(tracker;force=true)
        end
        @test TOML.parsefile(path)["revision"]==tracker.revision==11
        @test length(readdir(dirname(path)))==1
        # A failed replace retains the previous complete destination.
        folder=mkpath(joinpath(root,"destination"))
        write(joinpath(folder,"retained"),"complete")
        @test_throws Base.IOError Gauntlet._write_progress(folder,Dict("new"=>1))
        @test read(joinpath(folder,"retained"),String)=="complete"
        blocker=joinpath(root,"blocker"); write(blocker,"file")
        bad=Gauntlet.CampaignProgress(blocker,[:a],"bad";io=IOBuffer())
        @test_logs (:warn,r"snapshot disabled") Gauntlet._publish_progress!(bad;force=true)
        @test bad.snapshot_failed
        @test_logs Gauntlet._publish_progress!(bad;force=true)
        if Sys.isunix()
            required=joinpath(root,"required")
            sessions=mkpath(joinpath(required,"sessions"))
            write(joinpath(sessions,"required.jld2"),"session fixture")
            chmod(sessions,0o500)
            try
                @test_throws Exception Gauntlet._with_campaign_progress(root,[:a],"required";io=IOBuffer()) do tracker
                    Gauntlet._record_campaign_wall(required,(id="required",);progress=:auto) do
                        42
                    end
                end
                @test TOML.parsefile(joinpath(root,"sessions","required.progress.toml"))["termination"]=="aborted"
            finally
                chmod(sessions,0o700)
            end
        end
    end
end

@testitem "Gauntlet / progress / watcher cache, session ownership and six rows" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using TOML
    mktempdir() do root
        clock=Ref(0.0)
        tracker=Gauntlet.CampaignProgress(root,[:a,:b],"watch";clock=()->clock[],wall_clock=()->100.,io=IOBuffer())
        emit(e)=Gauntlet._progress_event!(tracker,e)
        emit((kind=:benchmark,benchmark=:a,attempt="one",state=:running))
        emit((kind=:operand,benchmark=:a,attempt="one",role=:reference,state=:running,backend="opaque",mode="ordinary"))
        tracker.rows[1]["case"]="长名字\n\t\e[31m"*repeat("long",30)
        Gauntlet._publish_progress!(tracker;force=true)
        watch=Gauntlet.CampaignWatch(;session="watch")
        Gauntlet._watch_poll!(watch,root,0.,120.)
        @test watch.age==20.
        received=watch.received
        Gauntlet._watch_poll!(watch,root,1.,121.)
        @test watch.received==received
        lines=Gauntlet._watch_lines(watch.snapshot,79;advance=21.)
        @test length(lines)==6
        @test all(line->textwidth(line)<=79,lines)
        @test all(line->!any(iscntrl,line),lines)
        @test occursin("00:00:21",lines[4])
        @test occursin("21s ago",lines[6])
        output=IOBuffer()
        Gauntlet._watch_render!(watch,output,lines,(24,80);terminal=true)
        take!(output)
        next=Gauntlet._watch_lines(watch.snapshot,79;advance=21.,spinner='/')
        Gauntlet._watch_render!(watch,output,next,(24,80);terminal=true)
        update=String(take!(output))
        @test !occursin("Benchmarks",update)
        @test ncodeunits(update)<30
        @test !occursin("\e[2K",update)
        # The pinned invocation legitimately advances to another benchmark/attempt.
        emit((kind=:benchmark,benchmark=:a,attempt="one",state=:complete))
        emit((kind=:benchmark,benchmark=:b,attempt="two",state=:running))
        Gauntlet._watch_poll!(watch,root,2.,122.)
        @test watch.snapshot["active"]["benchmark"]=="b"
        @test watch.session=="watch"
        path=joinpath(root,"sessions","watch.progress.toml")
        last_good=watch.snapshot
        write(path,"broken = [")
        Gauntlet._watch_poll!(watch,root,3.,123.)
        @test watch.snapshot===last_good
        @test occursin("malformed",watch.notice)
        for mutate in (s->(s["eta"]["queued_seconds"]=NaN),
                s->(s["active"]["completed"]="bad"),s->(s["work_age_seconds"]=Inf),
                s->(s["closed"]="yes"),s->(s["updated_unix_seconds"]="yesterday"))
            invalid=deepcopy(last_good); mutate(invalid)
            @test !Gauntlet._valid_watch_snapshot(invalid,"watch")
        end
        @test_throws ArgumentError Gauntlet._safe_session("../old")
        @test occursin("'\"'\"'",Gauntlet._shell_quote("a'b"))
        # Plain output cannot receive terminal control or spinner frames.
        plain=Gauntlet.CampaignWatch(;snapshot=watch.snapshot)
        output=IOBuffer()
        for time in (0.,0.25,0.5,1.)
            Gauntlet._watch_render!(plain,output,lines,(24,80);now=time)
        end
        data=String(take!(output))
        @test count(==('\n'),data)==1
        @test !occursin('\e',data)
        @test !occursin("Current run |",data)
        short=Gauntlet.CampaignWatch(;snapshot=watch.snapshot)
        Gauntlet._watch_render!(short,output,lines,(4,80);terminal=true)
        @test !occursin('\e',String(take!(output)))
        @test occursin("\e[32m",Gauntlet._watch_styled("OK 2 | Failed 1",2,true))
        @test occursin("\e[31m",Gauntlet._watch_styled("OK 2 | Failed 1",2,true))
        @test Gauntlet._watch_styled("ETA --",4,false)=="ETA --"
    end
end

@testitem "Gauntlet / progress / conservative discovery and legacy data" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using TOML
    mktempdir() do root
        watch=Gauntlet.CampaignWatch(;session="legacy")
        legacy=Dict{String,Any}("schema"=>1,"session"=>"legacy","elapsed_seconds"=>10.0,
            "updated_unix_seconds"=>100.0,"closed"=>true,"measurement_active"=>false,
            "benchmarks"=>[Dict("id"=>"a","state"=>"complete")],"counts"=>Dict("complete"=>1),
            "active"=>Dict("benchmark"=>"a","unit"=>"trials","completed"=>999),
            "campaign_eta_seconds"=>1.0)
        Gauntlet._write_progress(joinpath(root,"sessions","legacy.progress.toml"),legacy)
        Gauntlet._watch_poll!(watch,root,0.0,120.0)
        lines=Gauntlet._watch_lines(watch.snapshot,100;advance=100.,notice=watch.notice)
        @test occursin("ETA --",lines[4]) # Legacy closure has no execution-exhaustion fact.
        @test occursin("00:00:10",lines[4])
        @test !occursin("999",lines[5])
        @test occursin("Legacy",lines[6])
        legacy["updated_unix_seconds"]="bad"
        @test !Gauntlet._valid_watch_snapshot(legacy,"legacy")
        @test !Gauntlet._watch_color(IOContext(IOBuffer(),:color=>false))
        if Sys.isunix()
            Gauntlet._write_progress(joinpath(root,"campaign.toml"),Dict("schema"=>3,"benchmarks"=>["a","b"]))
            leases=IO[]
            try
                for (id,session) in (("a","first"),("b","second"))
                    directory=mkpath(joinpath(root,id))
                    Gauntlet._write_progress(joinpath(directory,"state.toml"),Dict("state"=>"running","session"=>session,"attempt"=>"one"))
                    lease=open(joinpath(directory,"execution.lock"),"a+");push!(leases,lease)
                    @test ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
                end
                auto=Gauntlet.CampaignWatch()
                Gauntlet._watch_poll!(auto,root,0.0,120.0)
                @test auto.session===nothing
                @test auto.snapshot["inventory"]
                @test occursin("multiple sessions",auto.notice)
                @test !auto.snapshot["eta"]["known"]
                Gauntlet._write_progress(joinpath(root,"b","state.toml"),Dict("state"=>"complete","session"=>"second"))
                Gauntlet._watch_poll!(auto,root,1.0,121.0)
                @test auto.session=="first"
                @test occursin("unavailable",auto.notice)
            finally
                foreach(close,leases)
            end
        end
    end
end

@testitem "Gauntlet / performance / compute scope, callbacks and legacy observations" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels, Logging, JLD2
    const calls=NamedTuple[]
    const delivered=Ref(0)
    struct TimedProbe<:LineCableModels.Grammar.AbstractFormulation end
    function LineCableModels.compute(problem::LineParametersProblem, ::TimedProbe; options = (;))
        push!(calls,
            (quiet = LineCableModels.performance_sample_active(),
                receiver = LineCableModels.progress_receiver(), callback = haskey(options, :on_result)))
        value=LineParameters(PhaseDomain, fill(1.0+0im, 2, 2, 2), fill(0.0+1im, 2, 2, 2), [
            1.0, 10.0])
        haskey(options, :on_result)&&options.on_result(problem, 1, value)
        return value
    end
    model=Gauntlet.load_case(:two_insulated_wires;
        variation = Gauntlet.ExactOverrides(frequencies = [1.0, 10.0]))
    options=(on_result = (args...)->(delivered[]+=1),)
    calculation=Gauntlet.BenchmarkCalculation(:probe, model.problem, TimedProbe(); options)
    direct=Gauntlet._execute(calculation)
    @test delivered[]==1
    @test direct.timing.scope === :compute_call_wall
    @test direct.timing.seconds >= 0
    empty!(calls)
    measured=Gauntlet._benchmark_owned(calculation, (samples = 3, seconds = 10.0))
    @test measured.samples==3
    @test length(calls)==4 # warmup plus three measured calls
    @test all(call->call.quiet && call.receiver===nothing && !call.callback, calls)
    @test delivered[]==1
    @test measured.scope === :compute_call_wall
    @test measured.policy.callbacks === false
    @test !LineCableModels.performance_sample_active()
    mktempdir() do root
        saved=Gauntlet._execute(calculation; directory = root, model)
        retained=Gauntlet.read_calculation(joinpath(root, "calculation.jld2"))
        @test retained.metadata.timing.scope === :compute_call_wall
        @test retained.metadata.timing.seconds == saved.timing.seconds
        @test isfile(joinpath(root, "timing.toml"))
        reused=Gauntlet._execute(calculation; directory = root, model)
        @test reused.reused
        @test reused.timing == saved.timing
    end
end
