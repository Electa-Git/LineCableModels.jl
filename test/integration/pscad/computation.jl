@testitem "PSCAD / compile failure retains diagnostics and an interruption unloads the project" tags=[:integration, :pscad] setup=[PSCADStandIn] begin
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    const S = PSCADStandIn
    frequencies = 10.0 .^ range(-1, 5; length = 101)
    # The exported project of one insulated wire, staged under `root`.
    function staged(root)
        design = build(CableDesign, "wire", terminal(:core,
            core(Material(:conductor, 1.72e-8, 1.0); r = 0.004),
            insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
        system = build(LineCableSystem, [design], [Pose2(0, -1)];
            connections = [Dict(:core => 1)], system_id = "stand-in")
        problem = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies)
        project = joinpath(mkpath(joinpath(root, "case")), "stand_in.pscx")
        write(project, P._pscad_inputs(problem, Formulation(:pscad),
            LineCableModels.DataModel.flatten.(system.designs)).project)
        return project
    end
    settings = Dict(row => Dict{String, Any}() for row in ("ground", "frequency", "configuration"))
    mktempdir() do root
        application = S.install(root)
        application.compile_error = "synthetic compile failure"
        project = staged(root)
        output = joinpath(dirname(project), "outputs")
        caught = try
            P.run_remote_pscad(P.RemoteConfig(S.remote(root)), project, output, settings,
                frequencies; line_name = "lcm_failed")
            nothing
        catch error
            error
        end
        @test caught isa ErrorException
        message = sprint(showerror, caught)
        @test occursin("during line-constants computation", message)
        @test occursin("synthetic compile failure", message)
        @test occursin("PSCAD stand-in message", message)
        @test occursin("Full PSCAD log: $(joinpath(output, "pscad-console.txt"))", message)
        console = read(joinpath(output, "pscad-console.txt"), String)
        @test occursin("PSCAD stand-in message", console)
        @test occursin("stand-in project output", console)
        @test occursin("synthetic compile failure", console)
        @test isempty(S.loaded(application))
        @test !S.silence(application)
        @test all(S.ended, S.processes("forward"))
        @test_throws ArgumentError P.run_remote_pscad(P.RemoteConfig(S.remote(root)), project,
            output, settings, frequencies; line_name = "lcm_failed")
        missing = joinpath(root, "missing.log")
        @test P._diagnostic_tail(missing) == "PSCAD produced no diagnostic log."
        write(missing, "\n   \n")
        @test P._diagnostic_tail(missing) == "PSCAD diagnostic log is empty."
        write(missing, "one\n\ntwo\nthree\n")
        @test P._diagnostic_tail(missing; count = 2) == "two\nthree"
    end
    mktempdir() do root
        application = S.install(root)
        application.solves = false
        project = staged(root)
        console = joinpath(dirname(project), "outputs", "pscad-console.txt")
        task = @async try
            P.run_remote_pscad(P.RemoteConfig(S.remote(root)), project, dirname(console),
                settings, frequencies; line_name = "lcm_interrupted")
        catch error
            error
        end
        waiting = timedwait(() -> istaskdone(task) ||
            (isfile(console) && occursin("Waiting for", read(console, String))), 60)
        @test waiting === :ok
        @test !istaskdone(task)
        @test S.loaded(application) == ["stand_in"]
        @test S.silence(application)
        schedule(task, InterruptException(); error = true)
        @test fetch(task) isa InterruptException
        log = read(console, String)
        @test occursin("during detailed-output completion", log)
        @test occursin("PSCAD stand-in message", log)
        @test isempty(S.loaded(application))
        @test !S.silence(application)
        @test all(S.ended, S.processes("forward"))
    end
end

@testitem "PSCAD / step lines, completion records and compile-call timing scope" tags=[:integration, :pscad] setup=[PSCADStandIn] begin
    using Logging, TOML
    using PythonCall: pyconvert
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    const S = PSCADStandIn
    frequencies = 10.0 .^ range(-1, 7; length = 101)
    function staged(root)
        design = build(CableDesign, "wire", terminal(:core,
            core(Material(:conductor, 1.72e-8, 1.0); r = 0.004),
            insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
        system = build(LineCableSystem, [design], [Pose2(0, -1)];
            connections = [Dict(:core => 1)], system_id = "stand-in")
        problem = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies)
        project = joinpath(mkpath(joinpath(root, "case")), "stand_in.pscx")
        write(project, P._pscad_inputs(problem, Formulation(:pscad),
            LineCableModels.DataModel.flatten.(system.designs)).project)
        return project
    end
    settings = Dict(row => Dict{String, Any}() for row in ("ground", "frequency", "configuration"))
    # Each step line of a run is a line of its log and a debug record of the extension, as
    # is the completion record. A logger that admits debug records shows them. A run's
    # verbosity narrows the caller's logger and never widens it.
    narrowed(levels) = log -> LineCableModels.VerbosityLogger(log, levels)
    for (min_level, scope, shown) in ((Logging.Info, identity, false),
            (Logging.Debug, identity, true),
            (Logging.Debug, narrowed((default = 0,)), false),
            (Logging.Debug, narrowed((default = 0, LineCableModelsPSCADExt = 1)), false),
            (Logging.Debug, narrowed((default = 0, LineCableModelsPSCADExt = 2)), true),
            (Logging.Info, narrowed((default = 0, LineCableModelsPSCADExt = 2)), false))
        mktempdir() do root
            application = S.install(root)
            project = staged(root)
            output = joinpath(dirname(project), "outputs")
            log = Test.TestLogger(min_level = min_level)
            result = with_logger(scope(log)) do
                P.run_remote_pscad(P.RemoteConfig(S.remote(root)), project, output, settings,
                    frequencies; line_name = "lcm_steps")
            end
            @test result.elapsed_seconds ==
                  parse(Float64, read(joinpath(output, "timing.txt"), String)) >= 0
            @test result.elapsed_scope == P.PSCAD_TIMING_SCOPE
            @test result.console_path == joinpath(output, "pscad-console.txt")
            @test result.output_dir == output
            steps = filter(record -> record.message ==
                                     "Starting PSCAD line-constants computation", log.logs)
            @test length(steps) == (shown ? 1 : 0)
            @test all(record -> record.level == Logging.Debug && record._module === P, steps)
            @test occursin("Starting PSCAD line-constants computation",
                read(result.console_path, String))
            completions = filter(record -> record.message == "PSCAD frequency scan completed",
                log.logs)
            if !shown
                @test isempty(completions)
            else
                record = only(completions)
                @test record._module === P
                @test record.kwargs[:compile_call_seconds] == result.elapsed_seconds
                @test record.kwargs[:timing_scope] == P.PSCAD_TIMING_SCOPE
                @test !haskey(record.kwargs, :elapsed_seconds)
            end
            # PSCAD loads the project from the station's side of the shared directory, and
            # every file of the run's line moves from that directory into the outputs.
            @test pyconvert(Vector{String}, application.loads) ==
                  [raw"Z:\exchange\case\stand_in.pscx"]
            @test sort(readdir(output)) == ["lcm_steps.log", "lcm_steps_zs.out",
                "native-settings.toml", "pscad-console.txt", "result_ym.out", "result_yp.out",
                "result_zm.out", "result_zp.out", "solver.toml", "timing.txt"]
            @test !any(startswith("lcm_steps"), readdir(root))
            @test TOML.parsefile(joinpath(output, "solver.toml"))["version"] == "5.1.0"
            @test isempty(S.loaded(application))
            @test !S.silence(application)
        end
    end
end


@testitem "PSCAD / borrowed and started applications" tags=[:integration, :pscad] setup=[PSCADStandIn] begin
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    const S = PSCADStandIn
    frequencies = 10.0 .^ range(-1, 5; length = 101)
    function staged(root)
        design = build(CableDesign, "wire", terminal(:core,
            core(Material(:conductor, 1.72e-8, 1.0); r = 0.004),
            insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
        system = build(LineCableSystem, [design], [Pose2(0, -1)];
            connections = [Dict(:core => 1)], system_id = "stand-in")
        problem = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies)
        project = joinpath(mkpath(joinpath(root, "case")), "stand_in.pscx")
        write(project, P._pscad_inputs(problem, Formulation(:pscad),
            LineCableModels.DataModel.flatten.(system.designs)).project)
        return project
    end
    settings = Dict(row => Dict{String, Any}() for row in ("ground", "frequency", "configuration"))
    solve(root) = P.run_remote_pscad(P.RemoteConfig(S.remote(root)), staged(root),
        joinpath(root, "case", "outputs"), settings, frequencies; line_name = "lcm_application")
    # A PSCAD that listens on the station is borrowed: its silence returns, it keeps running,
    # and the run ends its forward.
    mktempdir() do root
        application = S.install(root)
        @test isfile(solve(root).console_path)
        @test !isfile(joinpath(S.STATION[], "started.txt"))
        @test isempty(S.processes("start"))
        @test !S.silence(application)
        @test S.quits(application) == 0
        @test length(S.processes("forward")) == 1
        @test all(S.ended, S.processes("forward"))
    end
    # Otherwise the run starts PSCAD in the shared directory and quits it.
    mktempdir() do root
        application = S.install(root; listening = false)
        @test isfile(solve(root).console_path)
        @test read(joinpath(S.STATION[], "started.txt"), String) ==
              "cmd /c \"cd /d Z:\\exchange && \"C:\\Program Files (x86)\\PSCAD\\5.1.0\\bin\\win64\\" *
              "Pscad.exe\" /startup:au /port:64339 /silence:true /splash:false\""
        @test S.quits(application) == 1
        @test length(S.processes("start")) == 1
        @test length(S.processes("forward")) == 1
        @test all(S.ended, [S.processes("forward"); S.processes("start")])
        @test isempty(S.loaded(application))
    end
end

@testitem "PSCAD / unsupported indexed equations fail without fallback" tags=[:integration, :pscad] begin
    using PythonCall
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    pscad = Formulation(:pscad)
    selected = pscad.formulas
    for (equation, selection) in ((P.earth_impedance, selected.earth_impedance.air),
            (P.earth_potential_coefficient, selected.earth_admittance)),
            (kind, s, t) in ((:self, 1, 2), (:mutual, 2, 3), (:self, 3, 3))
        caught = try
            equation(selection, Val(kind), Val(s), Val(t), pscad)
            nothing
        catch error
            error
        end
        @test caught isa ArgumentError
        @test occursin("source in layer $s and target in layer $t", sprint(showerror, caught))
    end
    @test_throws ArgumentError P.internal_impedance(selected.internal_impedance, Val(:invalid), pscad)
    @test_throws ArgumentError Formulation(:pscad; insulation_impedance=:not_registered)
    for (s, t) in ((1, 2), (2, 1))
        @test P.earth_impedance(Formulation(:pscad; earth_impedance=:ametani2009).formulas.earth_impedance, Val(:mutual), Val(s), Val(t), pscad) ==
            (EarthForm3 = (value = 0, readback = "AMETANIL"),)
        @test P.earth_impedance(Formulation(:pscad; earth_impedance=:lucca1994).formulas.earth_impedance, Val(:mutual), Val(s), Val(t), pscad) ==
            (EarthForm3 = (value = 2, readback = "LUCCA"),)
    end
    # PSCAD's coaxial cable admits coaxial topology only. The topology rule of flatten
    # rejects the pipe.
    copper = Material(:conductor, 1.72e-8)
    dielectric = Material(:insulator, 1e14, 2.3)
    conductors = [terminal(name, solid(copper, Disk(0.005)), insulation(dielectric; t = 0.002))
                  for name in (:a, :b)]
    wall = terminal(:pipe, sheath(copper; t = 0.001), insulation(dielectric; t = 0.002))
    enclosed = build(CableDesign, "pscad-pipe", pipe(at(conductors[1], -0.01, 0),
        at(conductors[2], 0.01, 0); shape = Disk(0.025), fill = Material(:insulator, Inf, 1.0), wall))
    coaxial = build(CableDesign, "pscad-coaxial", first(conductors))
    @test length(LineCableModels.DataModel.flatten(coaxial)) == 1
    @test_throws "pipe-type cables are not supported: the conductor :pipe encloses the " *
                 "coaxial groups (:a), (:b)" LineCableModels.DataModel.flatten(enclosed)
end

@testitem "PSCAD / mixed native equations survive coaxial author withdrawal" tags=[:integration, :pscad, :slow] begin
    using PythonCall
    const E = LineCableModels.Engine
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    selected = E.EarthImpedance.Formula(:ametani2009)
    @test occursin("mixed", description(selected))
    pscad = Formulation(:pscad; earth_impedance=:ametani2009)
    native = pscad.formulas.earth_impedance
    for (s,t) in ((1,2), (2,1))
        @test_throws r"not yet implemented" E.EarthImpedance.earth_impedance(
            selected, Val(:mutual), Val(s), Val(t), nothing, nothing)
        @test P.earth_impedance(native, Val(:mutual), Val(s), Val(t), pscad).EarthForm3.readback == "AMETANIL"
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
