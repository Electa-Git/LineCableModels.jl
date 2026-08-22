@testitem "Gauntlet / mode and snapshot behavior" tags=[:unit] setup=[
    GauntletSupport,
    TestFixtures
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport
    using .TestFixtures

    @test !isdefined(GauntletSupport, :PSCADHarness)
    @test all(package.name != "PythonCall" for package in keys(Base.loaded_modules))

    mktempdir() do directory
        source=joinpath(directory, "case.jl")
        write(source, "# deliberate case source\n")
        parameters=LineParameters(
            PhaseDomain,
            reshape(ComplexF64[1 + 2im, 2 + 3im], 1, 1, :),
            reshape(ComplexF64[2im, 4im], 1, 1, :),
            [1.0, 10.0]
        )
        case=GauntletCase(
            :fixture,
            source,
            nothing,
            nothing,
            :underground_wedepohl_pollaczek_lossless,
            ["core"],
            false,
            (1, 1, 2),
            (;),
            false
        )
        snapshot=joinpath(directory, "fixture.jld2")
        write_snapshot(
            snapshot,
            case,
            parameters;
            pscad_version = "5.1.0",
            pscad_elapsed_seconds = 1.25
        )
        loaded=load_snapshot(case; path = snapshot)
        @test loaded.reference isa LineParameters
        @test Z(loaded.reference) == Z(parameters)
        @test loaded.metadata["port_order"] == ["core"]
        @test loaded.metadata["kron_reduced"] === false

        write(source, "# changed source\n")
        @test_throws ArgumentError load_snapshot(case; path = snapshot)
        @test_throws ArgumentError load_snapshot(
            case; path = joinpath(directory, "missing.jld2")
        )

        write(source, "# executable snapshot case\n")
        problem=TestFixtures.line_parameters_problem(; frequencies = [50.0])
        formulation=Formulation()
        reference=compute!(problem, formulation)
        executable=GauntletCase(
            :executable,
            source,
            problem,
            formulation,
            :underground_wedepohl_pollaczek_lossless,
            ["phase:1", "phase:2", "phase:3"],
            true,
            size(Z(reference)),
            (;),
            false
        )
        executable_snapshot=joinpath(directory, "executable.jld2")
        write_snapshot(
            executable_snapshot,
            executable,
            reference;
            pscad_version = "5.1.0",
            pscad_elapsed_seconds = 1.25
        )
        outcome=run_snapshot(executable; path = executable_snapshot)
        @test outcome.mode === :snapshot
        @test outcome.comparison.Z == RMSError(0.0, 0.0)
        @test outcome.comparison.Y == RMSError(0.0, 0.0)
        @test !isdefined(GauntletSupport, :PSCADHarness)
    end

    previous=get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", nothing)
    previous_ci=get(ENV, "CI", nothing)
    try
        delete!(ENV, "CI")
        for mode in ("snapshot", "live", "record")
            ENV["LINECABLEMODELS_GAUNTLET_MODE"]=mode
            @test gauntlet_mode() == Symbol(mode)
        end
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="automatic"
        @test_throws ArgumentError gauntlet_mode()
        ENV["CI"]="true"
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="live"
        @test_throws ArgumentError gauntlet_mode()
    finally
        previous===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_MODE") :
        (ENV["LINECABLEMODELS_GAUNTLET_MODE"]=previous)
        previous_ci===nothing ? delete!(ENV, "CI") : (ENV["CI"]=previous_ci)
    end
end

@testitem "Gauntlet / PSCAD parser and fixed mappings" tags=[:unit] setup=[
    GauntletSupport
] begin
    using Base64: base64decode
    using Test
    using LineCableModels.Engine
    using .GauntletSupport

    harness=GauntletSupport._load_live_support!()
    @test harness.formulation_spec(:overhead_deri_carson_lossless).pscad_parameter ===
          :EarthForm2
    @test_throws ArgumentError harness.formulation_spec(:unknown)

    mktempdir() do directory
        frequency=[1.0, 10.0]
        magnitude=[1.0 2.0 3.0 4.0
                   5.0 6.0 7.0 8.0]
        phase=zeros(2, 4)
        function table(path, values)
            open(path, "w") do io
                println(io, " LOG10(FN) FN V11 V12 V21 V22")
                for index in eachindex(frequency)
                    println(
                        io,
                        join(
                            (log10(frequency[index]),
                                frequency[index], values[index, :]...), ' ')
                    )
                end
            end
        end
        table(joinpath(directory, "result_zm.out"), magnitude)
        table(joinpath(directory, "result_zp.out"), phase)
        table(joinpath(directory, "result_ym.out"), magnitude .* 1.0e-6)
        table(joinpath(directory, "result_yp.out"), phase)
        result=harness.read_pscad_result(directory, frequency, (2, 2, 2))
        @test result isa LineParameters
        @test size(Z(result)) == (2, 2, 2)
        @test Z(result)[:, :, 1] == ComplexF64[1 2; 3 4]
        @test Y(result)[:, :, 2] == ComplexF64[5 6; 7 8] .* 1.0e-6
        rounded_frequency=[1.0, 10.0 * (1 + 4.0e-8)]
        rounded_result=harness.read_pscad_result(
            directory, rounded_frequency, (2, 2, 2)
        )
        @test frequencies(rounded_result) == rounded_frequency
        @test_throws ArgumentError harness.read_pscad_result(
            directory, [1.0, 11.0], (2, 2, 2)
        )
        mv(joinpath(directory, "result_zm.out"), joinpath(directory, "raw_zm.out"))
        @test_throws ArgumentError harness.read_pscad_result(directory, frequency, (2, 2, 2))
    end

    config=harness.RemoteConfig(
        "host",
        raw"Z:\gauntlet\cases\.work",
        raw"C:\gauntlet",
        "julia",
        "python";
        transport = :ssh,
        verbosity = 2,
        timeout_seconds = 60
    )
    powershell="[IO.Directory]::CreateDirectory('C:\\gauntlet') | Out-Null"
    command=harness.remote_command(config, powershell)
    @test command.exec[1] == "ssh"
    @test "host" in command.exec
    @test "-EncodedCommand" in command.exec
    encoded=command.exec[end]
    bytes=base64decode(encoded)
    decoded=transcode(String, ltoh.(collect(reinterpret(UInt16, bytes))))
    @test decoded == powershell
    @test !occursin("Out-Null", encoded)
    @test harness._remote_project_name(raw"C:\gauntlet\case\generated.pscx") ==
          "generated"
    @test_throws ArgumentError harness.RemoteConfig(
        "host", "shared", "remote", "julia", "python"; verbosity = 3
    )
    @test_throws ArgumentError harness.RemoteConfig(
        "host", "shared", "remote", "julia", "python"; timeout_seconds = 0
    )
    @test_throws ArgumentError harness._supervisor_command(
        config,
        raw"Z:\gauntlet\cases\.work\case\current",
        raw"C:\gauntlet\case\current",
        "case",
        :overhead_deri_carson_lossless,
        [1.0, 3.0, 10.0]
    )
    @test_throws ArgumentError harness._validate_frequencies(
        collect(10.0 .^ range(0, stop = 6, length = 61))
    )
    frequency_probe=collect(10.0 .^ range(0, stop = 6, length = 101))
    @test harness._validate_frequencies(frequency_probe) === frequency_probe
    supervisor_command=harness._supervisor_command(
        config,
        raw"Z:\gauntlet\cases\.work\case\current",
        raw"C:\gauntlet\case\current",
        "generated",
        :overhead_deri_carson_lossless,
        frequency_probe
    )
    @test occursin("supervisor.ps1", supervisor_command)
    @test occursin(raw"Z:\gauntlet\cases\.work\case\current", supervisor_command)
    @test occursin(raw"C:\gauntlet\case\current", supervisor_command)
    @test occursin("-Verbosity '2'", supervisor_command)
    @test occursin("-TimeoutSeconds '60'", supervisor_command)
    @test !occursin("OpenStandardInput", supervisor_command)
    cancel_command=harness._cancel_command(raw"C:\gauntlet\case\current")
    @test occursin("owner.txt", cancel_command)
    @test occursin("taskkill.exe /PID", cancel_command)
    @test occursin("CommandLine.IndexOf", cancel_command)

    mktempdir() do directory
        variant=joinpath(directory, "case", "current")
        output=joinpath(variant, "outputs")
        project=joinpath(variant, "generated.pscx")
        mkpath(variant)
        write(project, "fixture")
        staged=harness._stage_toolkit(project, output)
        @test staged == joinpath(variant, "toolkit")
        @test sort(readdir(staged)) == [
            "Manifest.toml",
            "Project.toml",
            "runner.jl",
            "supervisor.ps1"
        ]
    end

    Core.eval(
        harness,
        :(function remote_command(
                ::Val{:fixture_fail}, config::RemoteConfig, powershell::AbstractString
        )
            Cmd(["julia", "--startup-file=no", "-e", "exit(7)"])
        end)
    )
    Core.eval(
        harness,
        quote
            function remote_command(
                    ::Val{:fixture_logs},
                    config::RemoteConfig,
                    powershell::AbstractString
            )
                executable=Base.julia_cmd().exec[1]
                return Cmd([
                    executable,
                    "--startup-file=no",
                    "-e",
                    "println(\"runner output\"); println(stderr, \"runner error\")"
                ])
            end
        end
    )
    failed_config=harness.RemoteConfig(
        "host", "shared", raw"C:\gauntlet", "julia", "python";
        transport = :fixture_fail
    )
    mktempdir() do directory
        stdout_path=joinpath(directory, "stdout.txt")
        stderr_path=joinpath(directory, "stderr.txt")
        @test_throws ErrorException harness._run_remote(
            failed_config, "ignored"; stdout_path, stderr_path
        )
        @test isfile(stdout_path)
        @test isfile(stderr_path)
    end
    mktempdir() do directory
        logs_config=harness.RemoteConfig(
            "host", "shared", raw"C:\gauntlet", "julia", "python";
            transport = :fixture_logs
        )
        stdout_path=joinpath(directory, "stdout.txt")
        stderr_path=joinpath(directory, "stderr.txt")
        output=harness._run_remote(
            logs_config, "ignored"; stdout_path, stderr_path
        )
        @test occursin("runner output", output)
        @test read(stdout_path, String) == "runner output\n"
        @test read(stderr_path, String) == "runner error\n"
    end

    supervisor=read(
        joinpath(
            GauntletSupport.GAUNTLET_ROOT,
            "pscad",
            "remote",
            "supervisor.ps1"
        ),
        String
    )
    @test occursin("Stop-OwnedRunner", supervisor)
    @test occursin("taskkill.exe /PID", supervisor)
    @test occursin("TimeoutSeconds", supervisor)
    @test occursin("Copy-Result", supervisor)
    @test occursin("exit-code.txt", supervisor)
    @test occursin("%ERRORLEVEL%", supervisor)
    @test !occursin("process.ExitCode", supervisor)
    @test !occursin("taskkill.exe /IM", supervisor)
    runner=read(
        joinpath(GauntletSupport.GAUNTLET_ROOT, "pscad", "remote", "runner.jl"),
        String
    )
    @test occursin("Starting PSCAD line-constants calculation", runner)
    @test occursin("line.compile()", runner)
    @test !occursin("line._command", runner)
end

@testitem "Gauntlet / private legacy PSCAD exporter" tags=[:unit] setup=[
    GauntletSupport,
    TestFixtures
] begin
    using Test
    using LineCableModels.ImportExport
    using .GauntletSupport
    using .TestFixtures

    harness=GauntletSupport._load_live_support!()
    system=TestFixtures.three_phase_system()
    problem=TestFixtures.line_parameters_problem()
    earth=problem.earth_props
    mktempdir() do directory
        current=export_data(
            :pscad,
            system,
            earth;
            file_name = joinpath(directory, "current.pscx")
        )
        legacy=harness.export_legacy_pscad(
            system,
            earth;
            file_name = joinpath(directory, "legacy.pscx")
        )
        @test isfile(current)
        @test isfile(legacy)
        for path in (current, legacy)
            content=read(path, String)
            @test occursin("master:Line_FrePhase_Options", content)
            @test occursin("master:Cable_Coax", content)
            @test occursin("master:Line_Ground", content)
            @test occursin("name=\"Main\"", content)
            @test occursin("name=\"CableSystem\"", content)
        end
    end

    timing=harness.benchmark_local(
        (; problem, formulation = LineCableModels.Formulation()); samples = 1, seconds = 1
    )
    @test timing.samples == 1
    @test timing.minimum_seconds >= 0
    @test timing.median_seconds >= timing.minimum_seconds
    @test timing.bytes >= 0
    @test timing.allocations >= 0
end
