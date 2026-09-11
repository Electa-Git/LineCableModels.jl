@testitem "PSCAD / parser and fixed mappings" tags=[:integration] begin
    using Base64: base64decode
    using Test
    import LineCableModels
    using LineCableModels: description
    using LineCableModels.Engine

    harness=LineCableModels.PSCAD
    overhead=Formulation(:pscad; earth_impedance = :Gary1976)
    underground=Formulation(:pscad; earth_impedance = :WedepohlWilcox1973)
    @test overhead isa harness.PSCADFormulation
    @test underground isa harness.PSCADFormulation
    @test hasmethod(compute, Tuple{LineParametersProblem, harness.PSCADFormulation})
    @test !isdefined(harness,:Gauntlet)
    @test harness.earth_impedance(Val(:Gary1976), Val(:mutual), Val(1), Val(1), Val(:pscad)) ==
          (EarthForm2 = (value = 0, readback = "DERISEMLYEN"),)
    @test harness.earth_impedance(Val(:WedepohlWilcox1973), Val(:mutual), Val(2), Val(2), Val(:pscad)) ==
          (EarthForm = (value = 0, readback = "WEDEPOHL"),)
    @test EarthImpedance.formula_id(overhead.methods.earth_impedance) === :Gary1976
    @test EarthImpedance.formula_id(underground.methods.earth_impedance) ===
          :WedepohlWilcox1973
    @test LineCableModels.formula_id(overhead.methods.earth_admittance) === :default
    @test LineCableModels.formula_id(overhead.methods.insulation_admittance) === :default
    @test occursin("lossless", lowercase(description(overhead.methods.insulation_admittance)))
    @test occursin("PSCAD native earth admittance", harness._formulation_label(overhead))
    @test overhead.options == (reduce_bundle = false, kron_reduction = false,
        ideal_transposition = false,base_frequency=50.0)
    @test_throws ArgumentError Formulation(:pscad; options = (output_stem = "invalid",))
    @test !isdefined(harness, :NativeEarthAdmittance)
    @test !isdefined(harness, :NativeInsulationAdmittance)
    @test !isdefined(harness, :DirectNumericalIntegration)

    identifiers=(:Gary1976, :Carson1926, :Pollaczek1926, :WedepohlWilcox1973,
        :Saad1996, :Ametani2009, :Lucca1994)
    @test all(
        id -> Formulation(:pscad; earth_impedance = LineCableModels.formula(id)) isa
              harness.PSCADFormulation,
        identifiers)
    @test haskey(harness.earth_impedance(Val(:Saad1996), Val(:mutual), Val(2), Val(2), Val(:pscad)), :EarthForm)
    @test harness.earth_impedance(Val(:Lucca1994), Val(:mutual), Val(1), Val(2), Val(:pscad)).EarthForm3.readback == "LUCCA"
    @test harness.earth_impedance(Val(:Carson1926), Val(:mutual), Val(1), Val(1), Val(:pscad)).EarthForm2.value == 2
    @test harness.earth_impedance(Val(:Pollaczek1926), Val(:mutual), Val(2), Val(2), Val(:pscad)).EarthForm.value == 2
    @test_throws ArgumentError harness.earth_impedance(Val(:Pollaczek1926), Val(:mutual), Val(1), Val(1), Val(:pscad))
    @test_throws ArgumentError harness.earth_impedance(Val(:DirectNumericalIntegration), Val(:mutual), Val(1), Val(2), Val(:pscad))

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
        raw"Z:\gauntlet\benchmarks\.work",
        raw"C:\gauntlet",
        "julia",
        "python"; local_root=mktempdir(),
        transport = :ssh,
        timeout_seconds = 60
    )
    owner=harness.PSCADFormulation
    @test LineCableModels.formulation_options(owner, (;)) == overhead.options
    @test_throws ArgumentError LineCableModels.formulation_options(
        owner,
        (output_stem = "case",)
    )
    @test_throws ArgumentError LineCableModels.computation_options(owner, (;))
    execution=LineCableModels.computation_options(owner,
        (
            output_stem = "case",
            remote = config,
            verbosity = (default = 0, PSCAD = 2),
            output_basis = :total
        ))
    @test execution == (
        work_root=config.local_root,
        output_stem = "case",
        remote = config,
        verbosity = (default = 0, PSCAD = 2),
        output_basis = Val(:total),
        on_result = nothing,
        resume_run_directory = nothing,
        solver_identity = nothing
    )
    @test_throws ArgumentError LineCableModels.computation_options(
        owner, (
            output_stem = "benchmark_525kV_1600mm2_bipole_pscad",
            remote = config
        ))
    @test_throws MethodError LineCableModels.computation_options(
        Val(:pscad),
        (remote = config,)
    )
    powershell="[IO.Directory]::CreateDirectory('C:\\gauntlet') | Out-Null"
    command=harness.remote_command(config, powershell)
    @test command.exec[1] == "ssh"
    @test "host" in command.exec
    @test "-EncodedCommand" in command.exec
    encoded=command.exec[end]
    bytes=base64decode(encoded)
    decoded=transcode(String, ltoh.(collect(reinterpret(UInt16, bytes))))
    @test decoded == "\$ProgressPreference='SilentlyContinue'; $powershell"
    @test !occursin("Out-Null", encoded)
    @test harness._remote_project_name(raw"C:\gauntlet\case\generated.pscx") ==
          "generated"
    root=config.local_root
    @test LineCableModels.computation_options(owner,(remote=config,work_root=joinpath(root,"any","depth"))).work_root == joinpath(root,"any","depth")
    @test_throws ArgumentError LineCableModels.computation_options(owner,(remote=config,work_root=dirname(root)))
    @test_throws ArgumentError harness.RemoteConfig(
        "host", "shared", "remote", "julia", "python"; local_root=mktempdir(), timeout_seconds = 0
    )
    verbosity_error=try
        harness.RemoteConfig(
            "host", "shared", "remote", "julia", "python"; local_root=mktempdir(), verbosity = 2
        )
        nothing
    catch error
        error
    end
    @test verbosity_error isa ArgumentError
    @test occursin("options=(verbosity=", sprint(showerror, verbosity_error))
    @test_throws ArgumentError harness._supervisor_command(
        config,
        raw"Z:\gauntlet\benchmarks\.work\case\current",
        raw"C:\gauntlet\case\current",
        "case",
        overhead,
        [1.0, 3.0, 10.0];
        output_stem = "gauntlet",
        verbosity = 2
    )
    @test_throws ArgumentError harness._validate_frequencies(
        collect(10.0 .^ range(0, stop = 6, length = 61))
    )
    frequency_probe=collect(10.0 .^ range(0, stop = 6, length = 101))
    @test harness._validate_frequencies(frequency_probe) === frequency_probe
    minimum_frequency_probe=collect(10.0 .^ range(-1, stop = 6, length = 101))
    @test harness._validate_frequencies(minimum_frequency_probe) ===
          minimum_frequency_probe
    @test_throws DomainError harness._validate_frequencies(
        collect(10.0 .^ range(-2, stop = 6, length = 101))
    )
    supervisor_command=harness._supervisor_command(
        config,
        raw"Z:\gauntlet\benchmarks\.work\case\current",
        raw"C:\gauntlet\case\current",
        "generated",
        overhead,
        frequency_probe;
        output_stem = "gauntlet",
        verbosity = 2
    )
    @test occursin("supervisor.ps1", supervisor_command)
    @test occursin(raw"Z:\gauntlet\benchmarks\.work\case\current", supervisor_command)
    @test occursin(raw"C:\gauntlet\case\current", supervisor_command)
    @test occursin("-OutputStem 'gauntlet'", supervisor_command)
    @test occursin("-Verbosity '2'", supervisor_command)
    @test occursin("-TimeoutSeconds '60'", supervisor_command)
    @test !occursin("-EarthField", supervisor_command)
    @test !occursin("OpenStandardInput", supervisor_command)
    saad_command=harness._supervisor_command(
        config,
        raw"Z:\gauntlet\benchmarks\.work\case\current",
        raw"C:\gauntlet\case\current",
        "generated",
        Formulation(:pscad; earth_impedance = :Saad1996),
        frequency_probe;
        output_stem = "saad",
        verbosity = 0
    )
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
            "files.jl",
            "identity.py",
            "runner.jl",
            "supervisor.ps1"
        ]
        for (name, source) in harness.PSCAD_REMOTE_SOURCES
            @test read(joinpath(staged, name), String) == source
        end
    end

    mktempdir() do directory
        partial=joinpath(directory, "fixture_zm.out")
        write(partial, "header\n")
        writer=@async begin
            sleep(0.05)
            open(partial, "a") do io
                println(io, "0.0 1.0 2.0")
                println(io, "1.0 10.0 3.0")
            end
        end
        completed=harness._wait_output(
            [directory], "_zm.out", 2;
            timeout_seconds = 1,
            poll_seconds = 0.01
        )
        wait(writer)
        @test completed == realpath(partial)
        @test harness._data_rows(completed) == 2
    end
    mktempdir() do directory
        partial=joinpath(directory, "fixture_zm.out")
        write(partial, "header\n")
        error=try
            harness._wait_output(
                [directory], "_zm.out", 2;
                timeout_seconds = 0.05,
                poll_seconds = 0.01
            )
            nothing
        catch caught
            caught
        end
        @test error isa ArgumentError
        @test occursin("0 of 2 rows", sprint(showerror, error))
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
        "host", "shared", raw"C:\gauntlet", "julia", "python"; local_root=mktempdir(),
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
            "host", "shared", raw"C:\gauntlet", "julia", "python"; local_root=mktempdir(),
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

    for count in (101,201,501,1001), T in (Float32,Float64)
        f=T.(10.0.^range(-1,6;length=count))
        @test harness._validate_frequencies(f) === f
    end
    @test_throws ArgumentError harness._validate_frequencies(collect(range(.1,1e6;length=201)))
end
