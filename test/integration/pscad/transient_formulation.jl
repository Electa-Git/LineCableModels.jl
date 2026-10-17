@testitem "PSCAD / transient compute stages and resumes offline" tags=[:integration, :pscad] setup=[TestFixtures] begin
    using TOML, SHA, EzXML
    using LineCableModels.PSCAD: PSCADTransientFormulation
    const P = LineCableModels.PSCAD
    const launches = Ref(0)
    const identifications = Ref(0)
    const changed_settings = Ref(false)
    const missing_channel = Ref(false)
    const identity = Dict("schema" => "1", "version" => "5.1.0",
        "pscad_sha256" => repeat("1", 64), "line_constants_sha256" => repeat("2", 64),
        "master_library_sha256" => repeat("3", 64))

    function P.remote_command(::Val{:transient_compute_fixture}, config::P.RemoteConfig,
            command::AbstractString)
        if occursin("-SharedCase", command)
            launches[] += 1
            @test occursin("-Mode 'transient'", command)
            directory = replace(only(match(r"-SharedCase '([^']+)'", command).captures), '\\' => '/')
            script = """
                using TOML
                directory = $(repr(directory))
                input = TOML.parsefile(joinpath(directory, "computation.toml"))
                input["mode"] == "transient" || error("incorrect execution mode")
                identity = TOML.parse($(repr(sprint(TOML.print, identity))))
                get(input, "expected_solver", identity) == identity || error("unexpected solver identity")
                output = joinpath(directory, "outputs")
                stem = splitext(input["project_settings"]["output_filename"])[1]
                names = ["sending:cable:1:core:v", "sending:cable:1:core:i",
                    "sending:cable:1:sheath:v", "sending:cable:1:sheath:i",
                    "sending:cable:2:sheath:v", "sending:cable:2:sheath:i"]
                $(missing_channel[]) && (names[2] = "sending:cable:2:core:i")
                open(joinpath(output, stem * ".inf"), "w") do io
                    for (index, name) in enumerate(names)
                        println(io, "PGB(", index, ") Output Desc=", repr(name),
                            " Group=", repr("Main"), " Max=1000 Min=0 Units=", repr("V"))
                    end
                end
                open(joinpath(output, stem * "_01.out"), "w") do io
                    for t in range(0.0, input["project_settings"]["time_duration"]; length=201)
                        println(io, t, " 1000 2 30 3 40 4")
                    end
                end
                readback = Dict{String, Any}(component => Dict(field => control["readback"]
                    for (field, control) in controls) for (component, controls) in input["native_settings"])
                readback["project"] = input["project_settings"]
                $(changed_settings[]) && (readback["project"]["time_step"] = 99.0)
                for (name, value) in (("native-settings.toml", readback), ("solver.toml", identity),
                        ("execution.toml", Dict("mode" => "transient", "stem" => stem)))
                    open(io -> TOML.print(io, value), joinpath(output, name), "w")
                end
                write(joinpath(output, "timing.txt"), "1.25")
                write(joinpath(output, "pscad-console.txt"), "offline transient transport fixture")
                """
        else
            @test occursin("identify.py", command)
            identifications[] += 1
            script = "print(" * repr(sprint(TOML.print, identity)) * ")"
        end
        return `$(Base.julia_cmd()) --startup-file=no --project=@stdlib -e $script`
    end

    design = TestFixtures.coaxial_design()
    system = build(LineCableSystem, [design, design], [(-0.08, -1.0), (0.08, -1.0)];
        connections = [Dict(:core => 1, :sheath => 2), Dict(:core => 3, :sheath => 4)],
        system_id = "transient-compute", line_length = 600.0)
    line = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0),
        frequencies = [1.0, 1e6])
    sending = (Step("cable:1:core", 1000.0; t0 = 2e-6),
        TerminalImpedance("cable:1:sheath"; R = 10),
        TerminalImpedance("cable:2:sheath"; R = 10))
    problem = TransientAnalysisProblem(line; sending, receiving = (),
        window = (tmax = 20e-6, step = 1e-7))
    formulation = PSCADTransientFormulation(fit = (FS = 1.0, FE = 2e6, Numf = 200))
    parent = mkpath(joinpath(@__DIR__, "..", "..", "..", ".linecablemodels", "test", "transient-compute"))
    mktempdir(parent) do directory
        remote = P.RemoteConfig("fixture", directory, "unused", "julia", "python";
            local_root = directory, transport = :transient_compute_fixture)
        options = (; remote, timing = true)
        for invalid in ((;), (remote = 1,), (; remote, unknown = true),
                (; remote, output_stem = "bad/path"), (; remote, timing = 1),
                (; remote, output_basis = :invalid), (; remote, resume_run_directory = ""),
                (; remote, solver_identity = 1), (; remote, work_root = dirname(directory)))
            @test_throws ArgumentError compute(problem, formulation; options = invalid)
        end
        shifted_line = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0),
            frequencies = [1.0, 1e6], sigma = 1.0)
        shifted = TransientAnalysisProblem(shifted_line; sending, receiving = (), window = problem.window)
        @test_throws ArgumentError compute(shifted, formulation; options)
        @test launches[] == identifications[] == 0
        coarse = TransientAnalysisProblem(line; sending, receiving = (),
            window = (tmax = 20e-6, step = 2e-6))
        @test_logs (:warn, r"one tenth of the shortest possible travel time") P.validate(coarse, PSCADTransientFormulation)
        @test_logs P.validate(problem, PSCADTransientFormulation)
        result = compute(problem, formulation; options)
        @test launches[] == 1
        @test identifications[] == 0
        @test size(result.voltage) == size(result.current) == (2, 4, 201)
        @test result.coordinates == ["cable:1:core", "cable:1:sheath", "cable:2:core", "cable:2:sheath"]
        @test result.voltage[1, 1, :] == fill(1000.0, 201)
        @test result.current[1, 1, :] == fill(2.0, 201)
        @test result.voltage[1, 2, :] == fill(30.0, 201)
        @test result.current[1, 4, :] == fill(4.0, 201)
        @test all(iszero, result.current[2, :, :])
        @test result.t ≈ collect(range(0.0, 20e-6; length = 201))
        retained = details(result).data
        source = retained.execution.source_run
        @test !retained.execution.reused
        @test retained.execution.solver_identity == identity
        @test retained.formulations.model === :fd_phase
        @test retained.formulations.fit == formulation.fit
        @test retained.timing.run_call_seconds == 1.25
        @test retained.timing.wall_seconds >= 0
        @test retained.window == problem.window
        @test retained.native_readback["project"]["time_step"] == 0.1
        @test retained.native_readback["frequency"]["FE"] == 2e6
        @test P.computation_details(typeof(formulation), result) === details(result)
        expected = export_data(:pscad, problem, formulation; file_name = joinpath(directory, "expected.pscx"))
        @test read(joinpath(source, "generated.pscx"), String) == read(expected, String)
        document = readxml(joinpath(source, "generated.pscx"))
        @test length(findall("//User[@defn='master:cable_interface']", document)) == 2
        @test length(findall("//User[@defn='master:pgb']", document)) == 6
        @test length(findall("//User[@defn='master:Line_FrePhase_Options']", document)) == 1
        @test TOML.parsefile(joinpath(source, "computation.toml"))["mode"] == "transient"
        @test Set(keys(TOML.parsefile(joinpath(source, "complete.toml"))["outputs"])) ==
            Set(("solver.toml", "native-settings.toml", "timing.txt", "execution.toml", "transient.inf", "transient_01.out"))
        original = Dict(joinpath(root, name) => (bytes2hex(open(sha256, joinpath(root, name))),
            mtime(joinpath(root, name))) for (root, _, names) in walkdir(source) for name in names)
        callbacks = Int[]
        reused = compute(problem, formulation; options = (; options..., resume_run_directory = source,
            on_result = (p, i, r) -> push!(callbacks, i)))
        @test launches[] == 1
        @test identifications[] == 1
        @test callbacks == [1]
        @test details(reused).data.execution.reused
        @test isempty(details(reused).data.timing)
        @test reused.voltage == result.voltage && reused.current == result.current
        @test original == Dict(path => (bytes2hex(open(sha256, path)), mtime(path)) for path in keys(original))
        latest = compute(problem, formulation; options = (; options..., timing = false, resume_run_directory = :latest))
        @test details(latest).data.execution.reused
        @test !haskey(details(latest).data, :timing)
        @test_throws ArgumentError compute(problem, PSCADTransientFormulation(fit = (FE = 3e6,));
            options = (; options..., resume_run_directory = source))
        @test_throws ArgumentError compute(problem, formulation;
            options = (; options..., output_stem = "other", resume_run_directory = source))
        mode = compute(problem, PSCADTransientFormulation(model = :fd_mode, fit = (FT = 500.0,));
            options = (; remote, output_stem = "mode"))
        @test details(mode).data.formulations.model === :fd_mode
        @test details(mode).data.native_readback["frequency"]["FT"] == 500
        @test isfile(joinpath(details(mode).data.execution.source_run, "outputs", "mode.inf"))
        fresh = compute(problem, formulation; options)
        @test launches[] == 3
        @test details(fresh).data.execution.source_run != source
        changed_settings[] = true
        @test_throws r"project setting time_step differs" compute(problem, formulation; options)
        changed_settings[] = false
        missing_channel[] = true
        @test_throws r"missing PSCAD channel sending:cable:1:core:i" compute(problem, formulation; options)
        missing_channel[] = false
        @test count(path -> isfile(joinpath(path, "complete.toml")),
            readdir(directory; join = true)) == 3
        open(joinpath(source, "outputs", "transient_01.out"), "a") do io
            println(io, "damaged")
        end
        @test_throws r"output integrity check failed" compute(problem, formulation;
            options = (; options..., resume_run_directory = source))
        @test launches[] == 5
    end
end
