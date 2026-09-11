@testitem "Gauntlet / progress / absolute counters, quiet spans and bounded ETA" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels, TOML, Logging
    const LCM=LineCableModels
    events=NamedTuple[]
    @test LCM.progress_receiver() === nothing
    LCM.with_progress(e->push!(events, e); scope = (benchmark = :a,)) do
        LCM.with_progress_scope(role = :reference, point = 2) do
            LCM.report_progress(LCM.progress_receiver(), (
                stage = :sampling, completed = 3, total = 10))
            @test_throws ErrorException LCM.with_performance_sample() do
                @test LCM.performance_sample_active()
                @test LCM.progress_receiver() === nothing
                @info "suppressed"
                @warn "suppressed"
                error("restore context")
            end
            @test !LCM.performance_sample_active()
            @test LCM.progress_receiver() !== nothing
        end
    end
    @test only(events).benchmark === :a
    @test only(events).role === :reference
    @test only(events).point == 2
    @test LCM.progress_receiver() === nothing
    @test_throws ArgumentError Gauntlet._progress_mode(:verbose)
    diagnostics=IOBuffer()
    LCM.with_performance_sample() do
        Logging.with_logger(LCM.Engine.ConsoleVerbosityLogger(
            ConsoleLogger(diagnostics), (default = 2,))) do
            @warn "owned logger must honor quiet samples"
        end
    end
    @test isempty(take!(diagnostics))
    mktempdir() do root
        clock=Ref(0.0)
        output=IOBuffer()
        tracker=Gauntlet.CampaignProgress(root, [:a, :b], "test";
            mode = :plain, io = output, clock = ()->clock[])
        @test all(row[role]["total"]==0 for row in tracker.rows
        for role in ("reference", "candidate"))
        for row in tracker.rows, role in ("reference", "candidate")

            row[role]["total"]=1
            row[role]["known"]=true
        end
        emit(event)=Gauntlet._progress_event!(tracker, event)
        emit((
            kind = :benchmark, benchmark = :a, state = :running, attempt = "attempts/one"))
        emit((kind = :operand, role = :reference, state = :running))
        emit((stage = :sampling, unit = :trials, completed = 0,
            total = 100, attempts = 0, rejected = 0))
        for index in 1:10
            clock[]=Float64(index)
            emit((stage = :sampling, unit = :trials, completed = index, total = 100,
                attempts = index+2, rejected = 2))
        end
        snapshot=Gauntlet._progress_snapshot(tracker)
        @test snapshot["stage_eta_seconds"] ≈ 90.0
        @test snapshot["campaign_eta_seconds"] == -1
        @test snapshot["unestimated_benchmarks"] == 2
        emit((stage = :sampling, unit = :trials, completed = 10,
            total = 100, attempts = 12, rejected = 2))
        @test length(tracker.rates)==11
        for index in 11:1000
            clock[]=Float64(index)
            emit((stage = :sampling, unit = :trials, completed = index, total = 2000))
        end
        @test length(tracker.rates)==32
        emit((point = 2, stage = :sampling, unit = :trials, completed = 0, total = 100))
        @test length(tracker.rates)==1
        @test Gauntlet._progress_snapshot(tracker)["stage_eta_seconds"] == -1
        tracker.estimates["a"]=2000.0
        tracker.estimates["b"]=3000.0
        @test Gauntlet._progress_snapshot(tracker)["campaign_eta_seconds"] ≈ 4000.0
        clock[]=2500.0
        @test Gauntlet._progress_snapshot(tracker)["campaign_eta_seconds"] == -1
        Gauntlet._progress_paint!(tracker; force = true)
        path=joinpath(root, "sessions", "test.progress.toml")
        @test isfile(path)
        @test TOML.parsefile(path)["active"]["benchmark"] == "a"
        before=read(path)
        Base.ScopedValues.with(Gauntlet._CAMPAIGN_TRACKER=>tracker) do
            LCM.with_progress(emit; scope = (benchmark = :a,)) do
                Gauntlet._performance_span(sample = 1, samples = 2, role = :reference) do
                    @test tracker.paused
                    @test LCM.progress_receiver() === nothing
                    paused=read(path)
                    Gauntlet._progress_paint!(tracker)
                    @test read(path)==paused
                    emit((stage = :should_not_publish,))
                    @test tracker.active["stage"] == "performance"
                end
            end
        end
        @test !tracker.paused
        emit((kind = :operand, role = :reference, state = :complete, reused = true))
        emit((kind = :operand, role = :candidate, state = :running))
        tracker.rows[1]["candidate"]["estimate_seconds"]=20.0
        clock[]+=1.0
        @test Gauntlet._progress_snapshot(tracker)["active_eta_scope"]=="operand"
        @test Gauntlet._progress_snapshot(tracker)["stage_eta_seconds"]≈19.0
        emit((kind = :benchmark, benchmark = :a, state = :failed))
        final=Gauntlet._progress_snapshot(tracker)
        @test final["counts"]["failed"]==1
        @test final["counts"]["pending"]==1
        @test final["benchmarks"][1]["reference"]["reused"]==1
        @test final["benchmarks"][1]["candidate"]["state"]=="failed"
        @test !occursin('\e', String(take!(output)))
        @test Gauntlet._progress_duration(NaN)=="estimating"
        @test Gauntlet._progress_duration(-1)=="estimating"
        @test all(length(line)<=25 for line in Gauntlet._progress_lines(final, 25))
        withenv("TERM"=>"xterm") do
            terminal_output=IOBuffer()
            display=Gauntlet.CampaignProgress(root,[:a],"terminal";
                io=terminal_output,terminal=true)
            Gauntlet._progress_paint!(display;force=true)
            Gauntlet._progress_paint!(display;force=true)
            @test occursin('\e',String(take!(terminal_output)))
            LCM.with_progress(event->Gauntlet._progress_event!(display,event)) do
                LCM.with_progress_output() do
                    @test display.lines==0
                    println(terminal_output,"visible warning")
                end
            end
            @test occursin("visible warning",String(take!(terminal_output)))
            Gauntlet._progress_event!(display,(kind=:native_console,))
            @test !display.redraw
            @test display.lines==0
        end
    end
end

@testitem "Gauntlet / progress / consistent persistent completion summaries" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels, TOML
    core = Formulation()
    @test Gauntlet._progress_backend(core) == "Owned"
    @test Gauntlet._progress_backend(Formulation(:pscad)) == "PSCAD"
    @test Gauntlet._progress_backend(Formulation(:LineCableModelsFEM)) == "FEM"
    @test Gauntlet._progress_backend(MonteCarlo(core; trials=3)) == "Owned / Monte Carlo"
    @test Gauntlet._progress_backend(LinearError(core)) == "Owned / LEP"
    @test Gauntlet._progress_backend([core, core]) == "Owned"
    @test Gauntlet._progress_backend(Formulation(earth_impedance=Grid((:default, :Xue2018)))) == "Owned"
    unbuilt = Gridspace{typeof(core)}((args...)->error("display must not call a formulation factory"),
        (Grid((1, 2)),))
    @test Gauntlet._progress_backend(unbuilt) == "Owned"
    for terminal in (false, true)
        withenv("TERM"=>"xterm") do
            mktempdir() do root
                output = IOBuffer()
                clock = Ref(0.0)
                tracker = Gauntlet.CampaignProgress(root, [:pscad_case, :fem_case, :uq_case], "summary";
                    io=IOContext(output, :displaysize=>(30, 120)), terminal, clock=()->clock[])
                for row in tracker.rows, role in ("reference", "candidate")
                    row[role]["known"] = true
                    row[role]["total"] = 1
                end
                emit(event) = Gauntlet._progress_event!(tracker, event)
                emit((kind=:benchmark, benchmark=:pscad_case, state=:running))
                emit((kind=:operand, role=:reference, state=:running, backend="PSCAD"))
                emit((stage=:compiling, backend=:pscad))
                @test tracker.active["backend"] == "PSCAD"
                Gauntlet._progress_paint!(tracker; force=true)
                take!(output)
                clock[] = 20.0
                done = (kind=:operand, role=:reference, state=:complete, seconds=20.0,
                    compute_seconds=18.0, timing_reused=false)
                emit(done)
                # Events only update state; the renderer owns all output.
                @test isempty(take!(output))
                emit((kind=:operand, role=:candidate, state=:running, backend="Owned"))
                emit((kind=:operand, role=:candidate, state=:complete, seconds=0.2,
                    compute_seconds=nothing, reused=true, timing_reused=true))
                emit((kind=:benchmark, benchmark=:pscad_case, state=:complete))
                emit((kind=:benchmark, benchmark=:fem_case, state=:running))
                emit((kind=:operand, role=:reference, state=:running, backend="FEM"))
                emit((stage=:solving, backend=:fem, unit=:jobs, completed=2, total=101))
                @test tracker.active["backend"] == "FEM"
                Gauntlet._progress_paint!(tracker) # Completions bypass ordinary throttling.
                printed = String(take!(output))
                @test count("Completed | reference | PSCAD", printed) == 1
                @test count("Completed | candidate | Owned", printed) == 1
                @test occursin("Execution wall 20.0s | Compute-call wall 18.0s", printed)
                @test occursin("not run (saved result)", printed)
                @test occursin("Reused 1 | Recovery yes", printed)
                @test !terminal || tracker.lines == length(Gauntlet._progress_lines(Gauntlet._progress_snapshot(tracker), 120))
                @test terminal || !occursin('\e', printed)
                # A repeated completion notification cannot duplicate a summary.
                emit(merge(done, (benchmark=:pscad_case,)))
                Gauntlet._progress_paint!(tracker; force=true)
                @test !occursin("Completed |", String(take!(output)))
                snapshot = TOML.parsefile(joinpath(root, "sessions", "summary.progress.toml"))
                @test snapshot["benchmarks"][1]["reference"]["execution_seconds"] == 20.0
                @test snapshot["benchmarks"][1]["reference"]["compute_seconds"] == 18.0
                # In-process native recovery still has a current compute-call wall
                # duration, unlike loading an already saved Gauntlet operand.
                emit((kind=:operand, benchmark=:fem_case, role=:reference, state=:complete,
                    seconds=2.0, compute_seconds=1.5, timing_reused=true))
                tracker.paused = true
                Gauntlet._progress_paint!(tracker)
                @test isempty(take!(output))
                tracker.paused = false
                Gauntlet._progress_paint!(tracker)
                printed = String(take!(output))
                @test occursin("Completed | reference | FEM", printed)
                @test occursin("Compute-call wall 1.5s", printed)
                @test occursin("Reused 0 | Recovery yes", printed)
                emit((kind=:operand, role=:candidate, state=:running, backend="Owned"))
                clock[] += 3.0
                emit((kind=:benchmark, benchmark=:fem_case, state=:failed))
                Gauntlet._progress_paint!(tracker)
                printed = String(take!(output))
                @test occursin("Failed | candidate | Owned", printed)
                @test occursin("Execution wall 3.0s | Compute-call wall unavailable", printed)
                emit((kind=:benchmark, benchmark=:uq_case, state=:running))
                emit((kind=:operand, role=:reference, state=:running, backend="Owned / Monte Carlo"))
                emit((kind=:benchmark, benchmark=:uq_case, state=:interrupted))
                Gauntlet._progress_paint!(tracker)
                printed = String(take!(output))
                @test occursin("Interrupted | reference | Owned / Monte Carlo", printed)
                @test !occursin("candidate |", printed) # Pending work was skipped, not executed.
                @test all(textwidth(line) <= 25 for line in
                    Gauntlet._progress_completion_lines(tracker.rows[1], "reference", 25))
            end
        end
    end
end

@testitem "Gauntlet / progress / fresh campaign learns and exposes partial coverage" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels, TOML
    definitions = [
        Gauntlet._catalogue_fem_benchmark(:two_insulated_wires, @__FILE__; frequencies=[1.0, 10.0]),
        Gauntlet._catalogue_fem_benchmark(:two_bare_wires, @__FILE__; frequencies=[1.0, 10.0]),
        Gauntlet._catalogue_pscad_benchmark(:two_insulated_wires, @__FILE__; frequencies=[1.0, 10.0]),
    ]
    mktempdir() do root
        clock = Ref(0.0)
        tracker = Gauntlet.CampaignProgress(root, getproperty.(definitions, :id), "learning";
            clock=()->clock[], io=IOBuffer())
        Gauntlet._progress_declarations!(tracker, definitions)
        emit(event) = Gauntlet._progress_event!(tracker, event)
        @test Gauntlet._progress_snapshot(tracker)["unestimated_benchmarks"] == 3
        emit((kind=:benchmark, benchmark=definitions[1].id, state=:running))
        emit((kind=:operand, role=:reference, state=:running))
        clock[] = 100.0
        emit((kind=:operand, role=:reference, state=:complete, seconds=100.0))
        @test tracker.rows[2]["reference"]["estimate_seconds"] > 0
        @test !haskey(tracker.rows[3]["reference"], "estimate_seconds") # PSCAD is still unobserved.
        emit((kind=:operand, role=:candidate, state=:running))
        clock[] = 120.0
        emit((kind=:operand, role=:candidate, state=:complete, seconds=20.0))
        clock[] = 125.0
        emit((kind=:benchmark, benchmark=definitions[1].id, state=:complete,
            finalization_seconds=5.0, reused=false))
        snapshot = Gauntlet._progress_snapshot(tracker)
        @test snapshot["campaign_eta_seconds"] == -1
        @test snapshot["estimated_remaining_seconds"] > 0
        @test snapshot["unestimated_benchmarks"] == 1
        @test any(contains("Estimated portion ~"), Gauntlet._progress_lines(snapshot, 200))
        @test !haskey(snapshot["benchmarks"][1]["reference"], "group")
        emit((kind=:benchmark, benchmark=definitions[2].id, state=:running))
        emit((kind=:operand, role=:reference, state=:running))
        expected = tracker.rows[2]["reference"]["estimate_seconds"]
        clock[] += 1.0
        @test Gauntlet._progress_snapshot(tracker)["stage_eta_seconds"] ≈ expected-1
        # A native recovered execution must not replace its cost with a cache hit.
        emit((kind=:operand, role=:reference, state=:complete,
            seconds=0.01, timing_reused=true))
        @test !haskey(tracker.rows[2]["reference"], "observed_seconds")
        @test tracker.rows[2]["reference"]["estimate_seconds"] == expected
        emit((kind=:operand, role=:candidate, state=:running))
        emit((kind=:operand, role=:candidate, state=:complete, seconds=0.01, reused=true))
        @test !haskey(tracker.rows[2]["candidate"], "observed_seconds")
        emit((kind=:benchmark, benchmark=definitions[2].id, state=:complete,
            finalization_seconds=0.01, reused=true))
        @test !haskey(tracker.rows[2]["finalization"], "observed_seconds")

        # A subsequent fresh invocation can seed its cohorts from exact compatible history.
        Gauntlet._write_toml(joinpath(root, string(definitions[1].id), "state.toml"), Dict(
            "fresh_timing_key"=>Gauntlet._progress_history_key(definitions[1]),
            "fresh_wall_seconds"=>125.0, "fresh_reference_seconds"=>100.0,
            "fresh_candidate_seconds"=>20.0))
        again = Gauntlet.CampaignProgress(root, getproperty.(definitions[1:2], :id), "history";
            clock=()->clock[], io=IOBuffer())
        Gauntlet._progress_declarations!(again, definitions[1:2])
        @test again.rows[1]["reference"]["estimate_source"] == "history"
        @test again.rows[2]["reference"]["estimate_source"] == "comparable calculations"
        @test Gauntlet._progress_snapshot(again)["campaign_eta_seconds"] > 125

        # Propagation families and worker settings must not borrow one another's timings.
        problem = definitions[1].model.nominal_problem
        core = Formulation()
        descriptor(form; options=(;)) = Gauntlet._progress_workload(
            Gauntlet.BenchmarkCalculation(:probe, problem, form; options), problem)
        mc = descriptor(MonteCarlo(core; trials=100))
        @test descriptor(MonteCarlo(core; trials=200)).work == 2mc.work
        @test descriptor(MonteCarlo(core; trials=200)).group == mc.group
        @test descriptor(LinearError(core)).group != mc.group
        @test descriptor(core).group != mc.group
        @test descriptor(Formulation(:LineCableModelsFEM); options=(frequency_workers=1,)).group !=
            descriptor(Formulation(:LineCableModelsFEM); options=(frequency_workers=2,)).group
        function nested(options; overrides=(;))
            calculation = Gauntlet.BenchmarkCalculation(:probe,
                ParametricProblem(problem, options), Combinatorial(LineCableModelsFEM());
                options=overrides)
            Gauntlet._progress_workload(calculation, problem)
        end
        @test nested((frequency_workers=1,)).group != nested((frequency_workers=2,)).group
        @test nested((frequency_workers=1,); overrides=(frequency_workers=2,)).group ==
            nested((frequency_workers=2,)).group
    end
end

@testitem "Gauntlet / progress / heterogeneous FEM ETA, recovery and concurrent tail" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using LineCableModels
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    mktempdir() do root
        small = joinpath(root, "small.msh")
        large = joinpath(root, "large.msh")
        write(small, "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Nodes\n1 100 1 100\n")
        write(large, "\$MeshFormat\n4.1 0 8\n\$EndMeshFormat\n\$Nodes\n1 400 1 400\n")
        @test FEM._mesh_work_size(small) == 100.0
        @test FEM._mesh_work_size(large) == 400.0
        write(joinpath(root, "binary.msh"), "\$MeshFormat\n4.1 1 8\n")
        @test FEM._mesh_work_size(joinpath(root, "binary.msh")) === nothing
        @test FEM._mesh_work_size(joinpath(root, "missing.msh")) === nothing
        # Three measurements of one second per unit; include a large queued mesh.
        rates = [1.0, 1.0, 1.0]
        work = [10.0, 40.0, 20.0]
        now = UInt64(5_000_000_000)
        active = [(1, UInt64(0)), (2, UInt64(0))]
        @test FEM._fem_remaining_seconds(work, rates[1:2], active, 3, 2, now) == -1
        @test FEM._fem_remaining_seconds(work, rates, active, 3, 2, now) == 35.0
        @test FEM._fem_remaining_seconds(work, rates, [(1, UInt64(0))], 2, 1, now) == 65.0
        @test FEM._fem_remaining_seconds([40.0], rates, [(1, UInt64(0))], 2, 2, now) == 35.0
        @test FEM._fem_remaining_seconds(work, rates, active, 3, 2, UInt64(50_000_000_000)) == -1
        # Recovered columns are absent from work and cannot inflate either ETA or the rate sample count.
        @test FEM._fem_remaining_seconds([20.0], rates, Tuple{Int,UInt64}[], 1, 2, now) == 20.0
        @test FEM._fem_remaining_seconds(Float64[], rates, Tuple{Int,UInt64}[], 1, 2, now) == 0.0
        clock = Ref(0.0)
        tracker = Gauntlet.CampaignProgress(root, [:a], "mesh"; clock=()->clock[], io=IOBuffer())
        emit(event) = Gauntlet._progress_event!(tracker, event)
        emit((kind=:benchmark, benchmark=:a, state=:running))
        emit((kind=:operand, role=:reference, state=:running))
        emit((stage=:solving, unit=:jobs, completed=3, total=101,
            estimate_throughput=false, remaining_seconds=35.0))
        clock[] = 5.0
        @test Gauntlet._progress_snapshot(tracker)["stage_eta_seconds"] == 30.0
        clock[] = 40.0
        @test Gauntlet._progress_snapshot(tracker)["stage_eta_seconds"] == -1
        emit((stage=:completed,))
        @test Gauntlet._progress_snapshot(tracker)["stage_eta_seconds"] == -1
    end
end

@testitem "Gauntlet / progress / read-only watch and optional IO failures" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using .GauntletSupport: Gauntlet
    using TOML, Logging
    mktempdir() do root
        Gauntlet._write_toml(joinpath(root, "campaign.toml"), Dict("schema"=>3, "benchmarks"=>["a"]))
        mkpath(joinpath(root, "a"))
        state_path=joinpath(root, "a", "state.toml")
        state=Dict("state"=>"running", "session"=>"watch", "attempt"=>"one")
        Gauntlet._write_toml(state_path, state)
        # An invalid numerical file would fail immediately if watch read results.
        write(joinpath(root, "a", "calculation.jld2"), "must not deserialize")
        tracker=Gauntlet.CampaignProgress(root, [:a], "watch"; io = IOBuffer())
        Gauntlet._progress_event!(tracker,
            (kind = :benchmark, benchmark = :a,
                state = :running, attempt = "one", stage = :solving))
        Gauntlet._progress_paint!(tracker; force = true)
        snapshot_path=joinpath(root, "sessions", "watch.progress.toml")
        snapshot=TOML.parsefile(snapshot_path)
        snapshot["updated_unix_seconds"]=0.0
        Gauntlet._write_toml(snapshot_path, snapshot)
        open(joinpath(root, "a", "execution.lock"), "a+") do lease
            held=Sys.iswindows() ?
                 ccall(:_locking, Cint, (Cint, Cint, Clong), fd(lease), 2, 1)==0 :
                 ccall(:flock, Cint, (Cint, Cint), fd(lease), 6)==0
            @test held
            before=read(state_path), read(snapshot_path)
            lines=Gauntlet._campaign_watch_lines(root, 160; now = 20.0)
            @test any(line->occursin("solving", line), lines)
            @test any(line->occursin("liveness not inferred", line), lines)
            @test (read(state_path), read(snapshot_path))==before
            snapshot["measurement_active"]=true
            Gauntlet._write_toml(snapshot_path, snapshot)
            @test any(line->occursin("paused during performance", line),
                Gauntlet._campaign_watch_lines(root, 160; now = 20.0))
            snapshot["active"]["attempt"]="obsolete"
            Gauntlet._write_toml(snapshot_path, snapshot)
            @test length(Gauntlet._campaign_watch_lines(root, 160))==1
        end
        @test length(Gauntlet._campaign_watch_lines(root, 160))==1
        @test only(Gauntlet.campaign_status(root; verify = false)).state === :interrupted
        broken=IOBuffer()
        close(broken)
        tracker.io=broken
        @test_logs (:warn, r"progress output disabled") Gauntlet._progress_paint!(tracker; force = true)
        @test tracker.output_failed
        @test_logs Gauntlet._progress_paint!(tracker; force = true)
        blocker=joinpath(root, "not-a-directory")
        write(blocker, "fixture")
        tracker=Gauntlet.CampaignProgress(blocker, [:a], "broken"; io = IOBuffer())
        @test_logs (:warn, r"snapshot disabled") Gauntlet._progress_paint!(tracker; force = true)
        @test tracker.snapshot_failed
        @test_logs Gauntlet._progress_paint!(tracker; force = true)
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
