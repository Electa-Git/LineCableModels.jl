@testitem "PSCAD / transport failure retains diagnostics and interruption cancels" tags=[:integration] begin
    const P = LineCableModels.PSCAD
    using Logging
    function P.remote_command(::Val{:local_failure_probe}, config::P.RemoteConfig, command::AbstractString)
        script = occursin("-SharedCase", command) ?
            "println(stderr, \"synthetic transport failure\"); exit(9)" :
            "write(" * repr(joinpath(config.local_root, "cancelled")) * ", \"cancelled\")"
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end
    function P.remote_command(::Val{:local_interrupt_probe}, config::P.RemoteConfig, command::AbstractString)
        # Model a transport command with default OS termination. Julia's native
        # SIGTERM diagnostic handler can deadlock while printing its own stack.
        script = "write(" * repr(joinpath(config.local_root, "started")) *
            ", \"started\"); println(\"running\"); println(stderr, \"diagnostic\"); flush(stdout); flush(stderr); sleep(30)"
        return `$(Base.julia_cmd()) --startup-file=no --handle-signals=no --project=@stdlib -e $script`
    end
    mktempdir() do root
        config = P.RemoteConfig("local-fixture", root, "scratch", "julia", "python";
            local_root = root, transport = :local_failure_probe)
        directory = mkpath(joinpath(root, "nested", "case"))
        project = joinpath(directory, "generated.pscx")
        write(project, "protocol fixture")
        output = joinpath(directory, "outputs")
        caught = try
            P.run_remote_pscad(config, project, output, Formulation(:pscad),
                10.0.^range(-1, 5; length = 101); output_stem = "fixture")
            nothing
        catch error
            error
        end
        @test caught isa ErrorException
        @test occursin("Full PSCAD diagnostics:", sprint(showerror, caught))
        @test occursin("Transport stderr:", sprint(showerror, caught))
        @test read(joinpath(output, "transport-stderr.txt"), String) == "synthetic transport failure\n"
        @test isfile(joinpath(root, "cancelled"))
        @test occursin("before the remote failure", read(joinpath(output, "pscad-console.txt"), String))
        @test_throws ArgumentError P.run_remote_pscad(config, project, output,
            Formulation(:pscad), 10.0.^range(-1, 5; length = 101); output_stem = "fixture")
        missing = joinpath(root, "missing.log")
        @test P._diagnostic_tail(missing) == "PSCAD produced no diagnostic log."
        write(missing, "\n   \n")
        @test P._diagnostic_tail(missing) == "PSCAD diagnostic log is empty."
        write(missing, "one\n\ntwo\nthree\n")
        @test P._diagnostic_tail(missing; count = 2) == "two\nthree"
    end
    for callback_fails in (false, true)
        mktempdir() do root
            config = P.RemoteConfig("local-fixture", root, "scratch", "julia", "python";
                local_root = root, transport = :local_interrupt_probe)
            cancelled = Ref(false)
            callback = () -> begin
                cancelled[] = true
                callback_fails && error("synthetic cancellation diagnostic")
            end
            log = Test.TestLogger()
            task = with_logger(log) do
                @async try
                    P._run_remote(config, "ignored";
                        stdout_path = joinpath(root, "out.log"),
                        stderr_path = joinpath(root, "err.log"),
                        stream = true, on_interrupt = callback)
                catch error
                    error
                end
            end
            ready = timedwait(() -> isfile(joinpath(root, "started")) || istaskdone(task), 15)
            @test ready === :ok
            @test !istaskdone(task)
            schedule(task, InterruptException(); error = true)
            @test fetch(task) isa InterruptException
            @test cancelled[]
            @test callback_fails == any(record -> occursin("cancellation could not be confirmed",
                string(record.message)), log.logs)
        end
    end
end

@testitem "PSCAD / progress protocol is independent of diagnostic verbosity" tags=[:integration] begin
    const P = LineCableModels.PSCAD
    events=NamedTuple[]
    function P.remote_command(::Val{:local_progress_probe}, config::P.RemoteConfig, command::AbstractString)
        script="println(\"LCM_PROGRESS_V1\\tstage\\tcompiling\"); println(\"LCM_PROGRESS_V1\\theartbeat\\t1\"); println(\"human diagnostics\")"
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end
    mktempdir() do root
        config=P.RemoteConfig("fixture",root,"scratch","julia","python";
            local_root=root,transport=:local_progress_probe)
        command()=P._supervisor_command(config,root,"scratch","fixture",
            Formulation(:pscad),10.0.^range(-1,7;length=101);output_stem="fixture",verbosity=0)
        @test !occursin("-TrackProgress",command())
        LineCableModels.with_progress(event->push!(events,event)) do
            @test occursin("-TrackProgress",command())
            P._run_remote(config,"ignored";stream=false,
                stdout_path=joinpath(root,"transport.txt"))
            @test only(filter(e->haskey(e,:stage),events)).stage === :compiling
            @test any(e->haskey(e,:heartbeat_unix_seconds),events)
            @test all(e->!haskey(e,:completed),events)
            count=length(events)
            @test P._remote_progress(LineCableModels.progress_receiver(),"LCM_PROGRESS_V1\tstage\tunknown")
            @test !P._remote_progress(LineCableModels.progress_receiver(),"compiling 50 Hz")
            @test length(events)==count
            LineCableModels.with_performance_sample() do
                @test !occursin("-TrackProgress",command())
                P._remote_progress(LineCableModels.progress_receiver(),"LCM_PROGRESS_V1\tstage\tcompiling")
            end
            @test length(events)==count
        end
        @test occursin("human diagnostics",read(joinpath(root,"transport.txt"),String))
        @test !occursin("LCM_PROGRESS_V1",read(joinpath(root,"transport.txt"),String))
    end
end

@testitem "PSCAD / completion verbosity and compile-call timing scope" tags=[:integration] begin
    using Logging
    const P = LineCableModels.PSCAD
    function P.remote_command(::Val{:local_completion_probe}, config::P.RemoteConfig, command::AbstractString)
        output = joinpath(config.local_root, "case", "outputs")
        script = "for name in (\"pscad-console.txt\", \"result_zm.out\", \"result_zp.out\", \"result_ym.out\", \"result_yp.out\"); write(joinpath(" *
            repr(output) * ", name), \"fixture\"); end; write(joinpath(" * repr(output) *
            ", \"timing.txt\"), \"0.0328471\")"
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end
    for level in 0:2
        mktempdir() do root
            config = P.RemoteConfig("fixture", root, "scratch", "julia", "python";
                local_root=root, transport=:local_completion_probe)
            directory = mkpath(joinpath(root, "case"))
            project = joinpath(directory, "generated.pscx")
            write(project, "protocol fixture")
            log = Test.TestLogger()
            events = NamedTuple[]
            result = with_logger(log) do
                LineCableModels.with_progress(e->push!(events, e)) do
                    P.run_remote_pscad(config, project, joinpath(directory, "outputs"),
                        Formulation(:pscad), 10.0.^range(-1, 7; length=101);
                        output_stem="fixture", verbosity=level)
                end
            end
            @test result.elapsed_seconds ≈ 0.0328471
            @test result.elapsed_scope == P.PSCAD_TIMING_SCOPE
            @test any(e->get(e, :stage, nothing) === :validating, events)
            completions = filter(record->record.message == "PSCAD frequency scan completed", log.logs)
            if level == 0
                @test isempty(log.logs)
            else
                record = only(completions)
                @test record.kwargs[:compile_call_seconds] ≈ 0.0328471
                @test record.kwargs[:timing_scope] == P.PSCAD_TIMING_SCOPE
                @test !haskey(record.kwargs, :elapsed_seconds)
            end
        end
    end
end

@testitem "PSCAD / unsupported indexed equations fail without fallback" tags=[:integration] begin
    const P = LineCableModels.PSCAD
    selected = Formulation(:pscad).methods
    for (equation, selection) in ((P.earth_impedance, selected.earth_impedance.air),
            (P.earth_potential_coefficient, selected.earth_admittance)),
            (kind, s, t) in ((:self, 1, 2), (:mutual, 2, 3), (:self, 3, 3))
        caught = try
            equation(selection, Val(kind), Val(s), Val(t), Val(:pscad))
            nothing
        catch error
            error
        end
        @test caught isa ArgumentError
        @test occursin("source in layer $s and target in layer $t", sprint(showerror, caught))
    end
    @test_throws ArgumentError P.internal_impedance(selected.internal_impedance, Val(:invalid), Val(:pscad))
    @test_throws ArgumentError Formulation(:pscad; insulation_impedance=:not_registered)
    for (s, t) in ((1, 2), (2, 1))
        @test P.earth_impedance(Formulation(:pscad; earth_impedance=:ametani2009).methods.earth_impedance, Val(:mutual), Val(s), Val(t), Val(:pscad)) ==
            (EarthForm3 = (value = 0, readback = "AMETANIL"),)
        @test P.earth_impedance(Formulation(:pscad; earth_impedance=:lucca1994).methods.earth_impedance, Val(:mutual), Val(s), Val(t), Val(:pscad)) ==
            (EarthForm3 = (value = 2, readback = "LUCCA"),)
    end
    const pipe = LineCableModels.Engine.PipeImpedance.Formula(:default)
    @test_throws ArgumentError Formulation(Val(:pscad), pipe, Val(:pipe))
end

@testitem "PSCAD / mixed native equations survive coaxial author withdrawal" tags=[:integration] begin
    const E = LineCableModels.Engine
    const P = LineCableModels.PSCAD
    selected = E.EarthImpedance.Formula(:ametani2009)
    @test occursin("mixed", description(selected))
    native = Formulation(:pscad; earth_impedance=:ametani2009).methods.earth_impedance
    for (s,t) in ((1,2), (2,1))
        @test_throws r"not yet implemented" E.EarthImpedance.earth_impedance(
            selected, Val(:mutual), Val(s), Val(t), nothing, nothing, nothing)
        @test P.earth_impedance(native, Val(:mutual), Val(s), Val(t), Val(:pscad)).EarthForm3.readback == "AMETANIL"
    end
    copper = Material(:conductor, 1.72e-8, 1.0)
    design = build(CableDesign, "mixed-order", terminal(:core,
        solid(copper, Disk(0.004)), insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
    choices = (air = :carson1926, earth = :pollaczek1926, mixed = :ametani2009)
    configuration = Formulation(earth_impedance = choices)
    function problem(poses)
        system = build(LineCableSystem, [design, design], poses;
            connections = [Dict(:core => 1), Dict(:core => 2)])
        LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
    end
    for poses in ([Pose2(0, 2), Pose2(0.75, -1)], [Pose2(0.75, -1), Pose2(0, 2)])
        @test_throws r"not yet implemented" compute(problem(poses),configuration)
    end
end
