@testitem "Gauntlet / mode and snapshot behavior" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport,
    TestFixtures
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using Pkg.Artifacts
    using .GauntletSupport
    using .TestFixtures

    @test isdefined(GauntletSupport, :PSCADBenchmarks)
    run_owner=Val(GauntletCase)
    run_defaults=LineCableModels.computation_options(run_owner, (;))
    @test run_defaults == (
        output_basis = :per_length,
        reference = (;),
        candidate = (;),
        benchmark = (samples = 10, seconds = 10.0)
    )
    routed=LineCableModels.computation_options(run_owner,
        (
            output_basis = :total,
            reference = (solver = :external,),
            candidate = (solver = :native,),
            benchmark = (samples = 2, seconds = 0.5)
        ))
    @test routed.reference == (solver = :external,)
    @test routed.candidate == (solver = :native,)
    @test routed.benchmark == (samples = 2, seconds = 0.5)
    @test_throws ArgumentError LineCableModels.computation_options(
        run_owner,
        (verbosity = (default = 0,),)
    )
    @test_throws ArgumentError LineCableModels.computation_options(
        run_owner,
        (reference = Dict{Symbol, Any}(),)
    )
    @test_throws ArgumentError LineCableModels.computation_options(
        run_owner,
        (benchmark = (samples = 0,),)
    )
    @test gauntlet_instrumented() == (
        !iszero(Base.JLOptions().code_coverage) ||
        !iszero(Base.JLOptions().malloc_log)
    )
    @test all(package.name != "PythonCall" for package in keys(Base.loaded_modules))

    function fixture_case(source; name = :fixture, backend = :fixture)
        problem=TestFixtures.line_parameters_problem(; frequencies = [50.0])
        problem.system.system_id=string(name)
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
        reference_problem=LineParametersProblem(
            problem.system;
            temperature = problem.temperature,
            earth_props = problem.earth_props,
            frequencies = problem.frequencies
        )
        reference_formulation=Formulation(
            :pscad;
            earth_impedance = EarthImpedance.Wedepohl()
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
            name,
            backend,
            source,
            reference_problem,
            reference_formulation,
            problem,
            formulation,
            ports,
            (count, count, 1),
            tolerances
        )
    end

    mktempdir() do directory
        source=joinpath(directory, "case.jl")
        write(source, "# deliberate case source\n")
        case=fixture_case(source)
        @test case.backend === :fixture
        @test_throws ArgumentError fixture_case(source; backend = :PSCAD)
        parameters=compute(case.problem, case.formulation)
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
            reference_execution = (
                backend = :fixture,
                version = "1.0",
                elapsed_seconds = 1.25,
                exit_code = 0
            ),
            artifact_root
        )
        @test isfile(joinpath(persisted.path, "snapshot.jld2"))
        @test isfile(joinpath(persisted.path, "snapshot.sha256"))
        @test persisted.backend === :fixture
        @test persisted.gauntlet_version == GAUNTLET_VERSION
        aggregate=report(:fixture; artifact_root)
        @test size(aggregate) == (1, 39)
        @test only(aggregate.case) == "fixture"
        @test only(aggregate.basis) == "per_length"
        @test only(aggregate.domain) == "PhaseDomain"
        @test only(aggregate.Z_zero_atol) == 1.0e-6
        @test only(aggregate.Z_rms_absolute) == 0.0
        @test ismissing(only(aggregate.Z_rms_relative))
        @test only(aggregate.Z_rms_relative_raw) == 0.0
        @test only(aggregate.Y_rms_absolute) == 0.0
        @test ismissing(only(aggregate.Y_rms_relative))
        @test only(aggregate.Y_rms_relative_raw) == 0.0
        @test only(aggregate.Y_zero_atol) == 1.0e-9
        @test only(aggregate.reference_seconds) == 1.25
        @test only(aggregate.snapshot_sha256) == persisted.snapshot_sha256
        @test_throws ArgumentError report(
            :fixture;
            artifact_root,
            zero_atol = (Z = -1.0, Y = 0.0)
        )
        collections=finalize_artifacts(; artifact_root, artifacts_toml)
        @test only(collections).backend === :fixture
        @test basename(only(collections).archive) ==
              "benchmarks-fixture-v$(GAUNTLET_VERSION).tar.gz"
        @test isfile(only(collections).archive)
        @test isfile(joinpath(only(collections).path, "report.jld2"))
        @test isfile(joinpath(only(collections).path, "report.tsv"))
        @test isfile(joinpath(only(collections).path, "report.sha256"))
        @test isequal(only(collections).report, aggregate)
        installed=artifact_path(Base.SHA1(only(collections).tree_hash))
        @test isfile(joinpath(installed, "report.jld2"))
        @test isfile(joinpath(installed, "report.tsv"))
        @test isfile(joinpath(installed, "report.sha256"))
        @test isfile(snapshot_path(case; artifacts_toml))
        @test_throws ArgumentError finalize_artifacts(; artifact_root, artifacts_toml)
        forced=finalize_artifacts(; artifact_root, artifacts_toml, force = true)
        @test only(forced).tree_hash == only(collections).tree_hash
        published=publish_artifact(
            :fixture,
            "file://$(abspath(only(collections).archive))";
            artifact_root,
            artifacts_toml
        )
        @test published.artifact == "gauntlet_fixture_v1_0_0"
        @test release_tag() == "gauntlet-v1.0.0"
        @test backend_archive_name(:fixture) == "benchmarks-fixture-v1.0.0.tar.gz"
        @test occursin("file://", read(artifacts_toml, String))

        loaded=load_snapshot(case; artifacts_toml)
        @test loaded.reference isa LineParameters
        @test loaded.accepted isa LineParameters
        @test Z(loaded.reference) == Z(parameters)
        @test loaded.metadata["port_order"] == case.port_order
        @test loaded.metadata["formulation"] isa NamedTuple
        @test loaded.metadata["formulation"] == formulation_record(case)
        @test loaded.metadata["reference_execution"] == (
            backend = :fixture,
            version = "1.0",
            elapsed_seconds = 1.25,
            exit_code = 0
        )
        @test load_prior_snapshot(case; artifacts_toml) !== nothing
        @test load_prior_snapshot(case; artifacts_toml, force = true) === nothing

        new_source=joinpath(directory, "new_case.jl")
        write(new_source, "# new case with no accepted baseline\n")
        new_case=fixture_case(new_source; name = :new_fixture)
        @test load_prior_snapshot(new_case; artifacts_toml) === nothing

        outcome=run_snapshot(
            case;
            artifacts_toml,
            options = (
                output_basis = :per_length,
                candidate = (output_basis = :total,),
                benchmark = (samples = 1, seconds = 1)
            )
        )
        @test outcome.mode === :snapshot
        @test outcome.comparison.Z.absolute ==
              zeros(size(Z(parameters), 1), size(Z(parameters), 2))
        @test outcome.regression.Y.relative ==
              zeros(size(Y(parameters), 1), size(Y(parameters), 2))
        @test all(package.name != "PythonCall" for package in keys(Base.loaded_modules))

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
        @test_throws ArgumentError report(:fixture; artifact_root)

        @test_throws ArgumentError prepare_artifacts(; artifact_root, artifacts_toml)
        prepare_artifacts(; artifact_root, artifacts_toml, force = true)
        @test !ispath(backend_stage(:fixture; artifact_root))
    end

    previous=get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", nothing)
    previous_ci=get(ENV, "CI", nothing)
    previous_persist=get(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST", nothing)
    previous_runner=get(ENV, "LINECABLEMODELS_GAUNTLET_RUNNER", nothing)
    previous_cleanup=get(ENV, "LINECABLEMODELS_GAUNTLET_CLEANUP", nothing)
    previous_force=get(ENV, "LINECABLEMODELS_GAUNTLET_FORCE", nothing)
    try
        delete!(ENV, "CI")
        for mode in ("snapshot", "live", "record")
            ENV["LINECABLEMODELS_GAUNTLET_MODE"]=mode
            @test gauntlet_mode() == Symbol(mode)
        end
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="automatic"
        @test_throws ArgumentError gauntlet_mode()
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_CLEANUP")
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_FORCE")
        @test !gauntlet_cleanup()
        @test !gauntlet_force()
        ENV["LINECABLEMODELS_GAUNTLET_CLEANUP"]="true"
        ENV["LINECABLEMODELS_GAUNTLET_FORCE"]="TRUE"
        @test gauntlet_cleanup()
        @test gauntlet_force()
        ENV["LINECABLEMODELS_GAUNTLET_FORCE"]="yes"
        @test_throws ArgumentError gauntlet_force()
        ENV["CI"]="true"
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="live"
        @test_throws ArgumentError gauntlet_mode()
        delete!(ENV, "CI")
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="record"
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST")
        source=tempname()
        write(source, "# persist guard\n")
        @test_throws ArgumentError run_case(fixture_case(source))
        ENV["LINECABLEMODELS_GAUNTLET_PERSIST"]="true"
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_RUNNER")
        @test_throws ArgumentError run_case(fixture_case(source))
    finally
        previous===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_MODE") :
        (ENV["LINECABLEMODELS_GAUNTLET_MODE"]=previous)
        previous_ci===nothing ? delete!(ENV, "CI") : (ENV["CI"]=previous_ci)
        previous_persist===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST") :
        (ENV["LINECABLEMODELS_GAUNTLET_PERSIST"]=previous_persist)
        previous_runner===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_RUNNER") :
        (ENV["LINECABLEMODELS_GAUNTLET_RUNNER"]=previous_runner)
        previous_cleanup===nothing ?
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_CLEANUP") :
        (ENV["LINECABLEMODELS_GAUNTLET_CLEANUP"]=previous_cleanup)
        previous_force===nothing ? delete!(ENV, "LINECABLEMODELS_GAUNTLET_FORCE") :
        (ENV["LINECABLEMODELS_GAUNTLET_FORCE"]=previous_force)
    end

    mktempdir() do directory
        work_root=joinpath(directory, ".work")
        mkpath(work_root)
        write(joinpath(work_root, "output.txt"), "disposable\n")
        @test cleanup_work(; work_root) == work_root
        @test !ispath(work_root)
    end
end

@testitem "Gauntlet / explicit ports and phase-domain results" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels: ModalDomain
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
            system = (system_id = "ports", cables = [(conn = collect(assignments),)]),
            frequencies = [1.0, 10.0]
        )
        formulation=(options = (kron_reduction = kron, reduce_bundle = bundle),)
        count=length(assignments)
        return GauntletCase(
            :ports,
            :fixture,
            @__FILE__,
            problem,
            (kind = :reference,),
            problem,
            formulation,
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
    @test occursin("modal Gauntlet comparison", sprint(showerror, error))
end

@testitem "Gauntlet / PSCAD parser and fixed mappings" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Base64: base64decode
    using Test
    using LineCableModels: description
    using LineCableModels.Engine
    using .GauntletSupport

    harness=GauntletSupport.PSCADBenchmarks
    overhead=Formulation(:pscad; earth_impedance = EarthImpedance.Deri())
    underground=Formulation(:pscad; earth_impedance = EarthImpedance.Wedepohl())
    @test overhead isa harness.PSCADFormulation
    @test underground isa harness.PSCADFormulation
    @test hasmethod(
        compute,
        Tuple{LineParametersProblem, harness.PSCADFormulation}
    )
    @test hasmethod(
        benchmark_metadata,
        Tuple{LineParametersProblem, harness.PSCADFormulation}
    )
    @test harness.pscad_field(overhead.earth_impedance) === :EarthForm2
    @test harness.pscad_value(overhead.earth_impedance) == 0
    @test harness.pscad_readback(overhead.earth_impedance) == "DERISEMLYEN"
    @test harness.pscad_field(underground.earth_impedance) === :EarthForm
    @test harness.pscad_readback(underground.earth_impedance) == "WEDEPOHL"
    @test description(overhead.earth_admittance) == "PSCAD native earth admittance"
    @test description(overhead.insulation_admittance) ==
          "PSCAD native insulation admittance"
    @test harness._formulation_label(overhead) ==
          "Deri-Semlyen/PSCAD native earth admittance/PSCAD native insulation admittance"
    @test overhead.options == (;)
    @test_throws ArgumentError Formulation(
        :pscad;
        earth_impedance = EarthImpedance.Wedepohl(),
        options = (output_stem = "525kV_bipole",)
    )
    @test !isdefined(EarthImpedance, :ReferenceEarthImpedance)
    @test !isdefined(EarthImpedance, :DirectNumericalIntegration)

    methods=(
        EarthImpedance.Deri(),
        harness.DirectNumericalIntegration(:overhead),
        EarthImpedance.Wedepohl(),
        harness.DirectNumericalIntegration(:underground),
        EarthImpedance.Saad(),
        EarthImpedance.Ametani(),
        EarthImpedance.Lucca()
    )
    @test all(
        method -> Formulation(:pscad; earth_impedance = method) isa
                  harness.PSCADFormulation, methods)
    @test all(
        method -> supertype(typeof(method)) ===
                  LineCableModels.Engine.EarthImpedanceFormulation,
        methods
    )
    @test harness.pscad_field(EarthImpedance.Saad()) === :EarthForm
    @test harness.pscad_readback(EarthImpedance.Lucca()) == "LUCCA"
    @test_throws MethodError EarthImpedance.Wedepohl()(:self)
    unsupported=Formulation(
        :pscad;
        earth_impedance = EarthImpedance.Carson()
    )
    @test unsupported isa harness.PSCADFormulation
    @test_throws MethodError harness.pscad_field(unsupported.earth_impedance)
    @test_throws ArgumentError harness.DirectNumericalIntegration(:mutual)

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
        timeout_seconds = 60
    )
    owner=Val(harness.PSCADFormulation)
    @test LineCableModels.formulation_options(owner, (;)) == (;)
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
        output_stem = "case",
        remote = config,
        verbosity = (default = 0, PSCAD = 2),
        output_basis = Val(:total)
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
    backend_output=joinpath(
        GauntletSupport.WORK_ROOT,
        "pscad",
        "benchmark_525kV_1600mm2_bipole_pscad",
        "reference",
        "outputs"
    )
    work_parts=harness._work_parts(backend_output)
    @test work_parts == [
        "pscad", "benchmark_525kV_1600mm2_bipole_pscad", "reference"
    ]
    @test harness._remote_path(config.shared_root, work_parts...) ==
          raw"Z:\gauntlet\cases\.work\pscad\benchmark_525kV_1600mm2_bipole_pscad\reference"
    @test_throws ArgumentError harness._work_parts(
        joinpath(GauntletSupport.WORK_ROOT, "case", "reference", "outputs")
    )
    @test_throws ArgumentError harness.RemoteConfig(
        "host", "shared", "remote", "julia", "python"; timeout_seconds = 0
    )
    verbosity_error=try
        harness.RemoteConfig(
            "host", "shared", "remote", "julia", "python"; verbosity = 2
        )
        nothing
    catch error
        error
    end
    @test verbosity_error isa ArgumentError
    @test occursin("options=(reference=(verbosity=", sprint(showerror, verbosity_error))
    @test_throws ArgumentError harness._supervisor_command(
        config,
        raw"Z:\gauntlet\cases\.work\case\current",
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
    supervisor_command=harness._supervisor_command(
        config,
        raw"Z:\gauntlet\cases\.work\case\current",
        raw"C:\gauntlet\case\current",
        "generated",
        overhead,
        frequency_probe;
        output_stem = "gauntlet",
        verbosity = 2
    )
    @test occursin("supervisor.ps1", supervisor_command)
    @test occursin(raw"Z:\gauntlet\cases\.work\case\current", supervisor_command)
    @test occursin(raw"C:\gauntlet\case\current", supervisor_command)
    @test occursin("-OutputStem 'gauntlet'", supervisor_command)
    @test occursin("-Verbosity '2'", supervisor_command)
    @test occursin("-TimeoutSeconds '60'", supervisor_command)
    @test occursin("-EarthField 'EarthForm2'", supervisor_command)
    @test occursin("-EarthValue '0'", supervisor_command)
    @test occursin("-EarthReadback 'DERISEMLYEN'", supervisor_command)
    @test !occursin("OpenStandardInput", supervisor_command)
    saad_command=harness._supervisor_command(
        config,
        raw"Z:\gauntlet\cases\.work\case\current",
        raw"C:\gauntlet\case\current",
        "generated",
        Formulation(:pscad; earth_impedance = EarthImpedance.Saad()),
        frequency_probe;
        output_stem = "saad",
        verbosity = 0
    )
    @test occursin("-EarthField 'EarthForm'", saad_command)
    @test occursin("-EarthValue '3'", saad_command)
    @test occursin("-EarthReadback 'SAAD'", saad_command)
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
            "runner.jl",
            "supervisor.ps1"
        ]
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
    @test occursin("\$ProgressPreference = \"SilentlyContinue\"", supervisor)
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
    @test occursin("_set!(line, \"Name\", output_stem)", runner)
    @test occursin("_wait_output(roots, suffix, expected_rows)", runner)
    @test occursin("PSCAD diagnostics at failure", runner)
    @test !occursin("FORMULATIONS", runner)
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

    instrumented=performance_comparison(
        accepted,
        timing,
        tolerance;
        instrumented = true
    )
    @test !instrumented.comparable
    @test instrumented.passes === nothing
end
