@testitem "Gauntlet / mode and snapshot behavior" tags=[:gauntlet_toolkit] setup=[
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

    function fixture_case(source)
        problem=TestFixtures.line_parameters_problem(; frequencies = [50.0])
        next_phase=1
        ports=String[]
        for (cable_index, position) in enumerate(problem.system.cables)
            for component_index in eachindex(position.conn)
                position.conn[component_index]=next_phase
                component=position.design_data.components[component_index]
                push!(ports, "cable:$cable_index:$(component.id)")
                next_phase+=1
            end
        end
        formulation=Formulation(
            earth_impedance = EarthImpedance.Pollaczek(),
            earth_admittance = EarthAdmittance.IdealGround(),
            insulation_admittance = InsulationAdmittance.Lossless(),
            options = (kron_reduction = false, reduce_bundle = false)
        )
        tolerances=(
            reference = (
                Z = (absolute = 1.0e-6, relative = 1.0e-3),
                Y = (absolute = 1.0e-9, relative = 1.0e-3)
            ),
            regression = (
                Z = (absolute = 1.0e-12, relative = 1.0e-9),
                Y = (absolute = 1.0e-15, relative = 1.0e-9)
            ),
            performance = (
                median_time_ratio = 1.20,
                bytes_ratio = 1.05,
                allocations_ratio = 1.05
            )
        )
        count=length(ports)
        return GauntletCase(
            :fixture,
            source,
            problem,
            formulation,
            :underground_wedepohl_pollaczek_lossless,
            ports,
            (count, count, 1),
            tolerances
        )
    end

    mktempdir() do directory
        source=joinpath(directory, "case.jl")
        write(source, "# deliberate case source\n")
        case=fixture_case(source)
        parameters=compute!(case.problem, case.formulation)
        comparison=compare(parameters, parameters)
        benchmark=(
            minimum_seconds = 1.0,
            median_seconds = 1.1,
            bytes = 100,
            allocations = 10,
            samples = 10,
            environment = (
                julia_version = string(VERSION), kernel = string(Sys.KERNEL),
                architecture = string(Sys.ARCH), threads = Threads.nthreads(),
                blas = "fixture"
            )
        )
        artifact_root=joinpath(directory, ".artifacts")
        artifacts_toml=joinpath(directory, "Artifacts.toml")
        persisted=persist_snapshot(
            case,
            parameters,
            parameters,
            comparison,
            benchmark;
            pscad_version = "5.1.0",
            pscad_elapsed_seconds = 1.25,
            artifact_root,
            artifacts_toml
        )
        @test isfile(joinpath(persisted.path, "snapshot.jld2"))
        @test isfile(joinpath(persisted.path, "snapshot.sha256"))
        @test isfile(joinpath(persisted.path, "fixture.tar.gz"))
        @test isfile(snapshot_path(case; artifacts_toml))
        published=publish_artifact(
            :fixture,
            "file://$(abspath(joinpath(persisted.path, "fixture.tar.gz")))";
            artifact_root,
            artifacts_toml
        )
        @test published.artifact == "pscad_gauntlet_fixture"
        @test occursin("file://", read(artifacts_toml, String))

        loaded=load_snapshot(case; artifacts_toml)
        @test loaded.reference isa LineParameters
        @test loaded.accepted isa LineParameters
        @test Z(loaded.reference) == Z(parameters)
        @test loaded.metadata["port_order"] == case.port_order
        @test loaded.metadata["formulation"] isa NamedTuple
        @test loaded.metadata["formulation"].pscad ==
              "underground_wedepohl_pollaczek_lossless"

        outcome=run_snapshot(
            case;
            artifacts_toml,
            benchmark_samples = 1,
            benchmark_seconds = 1
        )
        @test outcome.mode === :snapshot
        @test outcome.comparison.Z.absolute ==
              zeros(size(Z(parameters), 1), size(Z(parameters), 2))
        @test outcome.regression.Y.relative ==
              zeros(size(Y(parameters), 1), size(Y(parameters), 2))
        @test !isdefined(GauntletSupport, :PSCADHarness)

        write(source, "# changed source\n")
        @test_throws ArgumentError load_snapshot(case; artifacts_toml)
        empty_toml=joinpath(directory, "empty.toml")
        write(empty_toml, "")
        @test_throws ArgumentError snapshot_path(case; artifacts_toml = empty_toml)

        snapshot=joinpath(persisted.path, "snapshot.jld2")
        open(snapshot, "a") do io
            write(io, UInt8(0))
        end
        @test_throws ArgumentError load_snapshot(case; path = snapshot)
    end

    previous=get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", nothing)
    previous_ci=get(ENV, "CI", nothing)
    previous_persist=get(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST", nothing)
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
        delete!(ENV, "CI")
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="record"
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST")
        source=tempname()
        write(source, "# persist guard\n")
        @test_throws ArgumentError run_case(fixture_case(source))
    finally
        previous===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_MODE") :
        (ENV["LINECABLEMODELS_GAUNTLET_MODE"]=previous)
        previous_ci===nothing ? delete!(ENV, "CI") : (ENV["CI"]=previous_ci)
        previous_persist===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST") :
        (ENV["LINECABLEMODELS_GAUNTLET_PERSIST"]=previous_persist)
    end
end

@testitem "Gauntlet / explicit ports and phase-domain results" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels.Engine
    using .GauntletSupport

    tolerances=(
        reference = (
            Z = (absolute = 1.0, relative = 1.0),
            Y = (absolute = 1.0, relative = 1.0)
        ),
        regression = (
            Z = (absolute = 1.0, relative = 1.0),
            Y = (absolute = 1.0, relative = 1.0)
        ),
        performance = (
            median_time_ratio = 1.2,
            bytes_ratio = 1.05,
            allocations_ratio = 1.05
        )
    )
    function make_case(assignments; kron = false, bundle = false)
        problem=(
            system = (cables = [(conn = collect(assignments),)],),
            frequencies = [1.0, 10.0]
        )
        formulation=(options = (kron_reduction = kron, reduce_bundle = bundle),)
        count=length(assignments)
        return GauntletCase(
            :ports,
            @__FILE__,
            problem,
            formulation,
            :underground_wedepohl_pollaczek_lossless,
            ["port:$index" for index in 1:count],
            (count, count, 2),
            tolerances
        )
    end

    valid=make_case([1, 2, 3])
    @test validate_case(valid) === valid
    @test_throws ArgumentError make_case([1, 0, 3])
    @test_throws ArgumentError make_case([1, 1, 2])
    @test_throws ArgumentError make_case([1, 2, 4])
    @test_throws ArgumentError make_case([1, 2, 3]; kron = true)
    @test_throws ArgumentError make_case([1, 2, 3]; bundle = true)

    modal=LineParameters(
        ModalDomain,
        zeros(ComplexF64, 3, 3, 2),
        zeros(ComplexF64, 3, 3, 2),
        [1.0, 10.0]
    )
    error=try
        validate_structure(valid, modal)
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test occursin("does not provide modal Z and Y", sprint(showerror, error))
end

@testitem "Gauntlet / PSCAD parser and fixed mappings" tags=[:gauntlet_toolkit] setup=[
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
        rounded_frequency=[1.0, 10.0*(1+4.0e-8)]
        rounded_result=harness.read_pscad_result(
            directory, rounded_frequency, (2, 2, 2)
        )
        @test frequencies(rounded_result) == rounded_frequency
        @test_throws ArgumentError harness.read_pscad_result(
            directory, [1.0, 11.0], (2, 2, 2)
        )
        mv(joinpath(directory, "result_zm.out"), joinpath(directory, "raw_zm.out"))
        @test_throws ArgumentError harness.read_pscad_result(directory, frequency, (
            2, 2, 2))
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
    @test occursin("_verify_retained_ports(components)", runner)
    @test occursin("r\"^elim\\d+\$\"", runner)
    @test !occursin("parameters[\"LC\"]", runner)
    @test !occursin("parameters[\"LL\"]", runner)
    @test occursin("requests conductor elimination", runner)
    @test occursin("line.compile()", runner)
    @test !occursin("line._command", runner)
end

@testitem "Gauntlet / local performance comparison" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport,
    TestFixtures
] begin
    using Test
    using LineCableModels
    using .GauntletSupport
    using .TestFixtures

    problem=TestFixtures.line_parameters_problem()
    timing=benchmark_local(
        (; problem, formulation = LineCableModels.Formulation()); samples = 1, seconds = 1
    )
    @test timing.samples == 1
    @test timing.minimum_seconds >= 0
    @test timing.median_seconds >= timing.minimum_seconds
    @test timing.bytes >= 0
    @test timing.allocations >= 0
    @test timing.environment.julia_version == string(VERSION)

    tolerance=(median_time_ratio = 1.2, bytes_ratio = 1.05, allocations_ratio = 1.05)
    accepted=(;
        timing...,
        median_seconds = timing.median_seconds/1.1,
        bytes = max(timing.bytes, 1),
        allocations = max(timing.allocations, 1)
    )
    compared=performance_comparison(accepted, timing, tolerance)
    @test compared.comparable
    @test compared.passes isa Bool
    other_environment=(;
        accepted...,
        environment = (; accepted.environment..., julia_version = "different")
    )
    diagnostic=performance_comparison(other_environment, timing, tolerance)
    @test !diagnostic.comparable
    @test diagnostic.passes === nothing
end
