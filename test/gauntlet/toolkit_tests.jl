@testitem "Gauntlet / mode and snapshot behavior" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport,
    TestFixtures
] begin
    using Test
    using SHA
    using LineCableModels
    using LineCableModels.Engine
    using Pkg.Artifacts
    using .GauntletSupport
    using .TestFixtures

    @test isdefined(GauntletSupport, :PSCADBenchmarks)
    run_owner=Val(GauntletCase)
    run_defaults=LineCableModels.computation_options(run_owner, (;))
    @test run_defaults == (
        output_basis = :pul,
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
        next_phase=1
        ports=String[]
        connections=Vector{Int}[]
        for (cable_index, design) in enumerate(problem.system.designs)
            cable_connections=Int[]
            for terminal in design.terminal_order
                push!(cable_connections, next_phase)
                push!(ports, "cable:$cable_index:$terminal")
                next_phase+=1
            end
            push!(connections, cable_connections)
        end
        system=LineCableModels.build(
            LineCableModels.LineCableSystem,
            problem.system.designs,
            problem.system.positions;
            connections,
            environment = problem.system.environment,
            system_id = string(name),
            line_length = problem.system.line_length
        )
        problem=LineParametersProblem(
            system;
            temperature = problem.temperature,
            earth_props = problem.earth_props,
            frequencies = problem.frequencies
        )
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
        @test persisted.schema_version == SNAPSHOT_SCHEMA_VERSION
        aggregate=GauntletSupport.report(:fixture; artifact_root)
        @test size(aggregate) == (1, 40)
        @test only(aggregate.benchmark) == "fixture"
        @test only(aggregate.case) == "fixture"
        @test only(aggregate.basis) == "pul"
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
        @test_throws ArgumentError GauntletSupport.report(
            :fixture;
            artifact_root,
            zero_atol = (Z = -1.0, Y = 0.0)
        )
        collections=finalize_staging(; artifact_root)
        @test only(collections).collection === :fixture
        @test only(collections).schema_version == SNAPSHOT_SCHEMA_VERSION
        @test isfile(joinpath(only(collections).path, "report.jld2"))
        @test isfile(joinpath(only(collections).path, "report.tsv"))
        @test isfile(joinpath(only(collections).path, "report.sha256"))
        @test isequal(only(collections).report, aggregate)
        @test !isfile(artifacts_toml)

        package=package_collection(
            :fixture,
            v"1.0.0";
            reason = "Initial fixture baseline",
            git_commit = repeat("a", 40),
            artifact_root
        )
        @test package.collection === :fixture
        @test package.version == v"1.0.0"
        @test package.tag == "gauntlet-fixture-v1.0.0"
        @test package.artifact == "gauntlet_fixture"
        @test basename(package.archive) == "benchmarks-fixture-v1.0.0.tar.gz"
        @test isfile(package.archive)
        @test isfile(package.package_path)
        installed=artifact_path(Base.SHA1(package.tree_hash))
        @test isfile(joinpath(installed, "report.jld2"))
        @test isfile(joinpath(installed, "report.tsv"))
        @test isfile(joinpath(installed, "report.sha256"))
        @test isfile(joinpath(installed, "release.toml"))
        @test_throws ArgumentError package_collection(
            :fixture,
            v"1.0.0";
            reason = "Initial fixture baseline",
            git_commit = repeat("a", 40),
            artifact_root
        )
        forced=package_collection(
            :fixture,
            v"1.0.0";
            reason = "Initial fixture baseline",
            git_commit = repeat("a", 40),
            artifact_root,
            force = true
        )
        @test forced.tree_hash == package.tree_hash
        published=bind_published_artifact(
            :fixture,
            v"1.0.0",
            "file://$(abspath(forced.archive))";
            artifact_root,
            artifacts_toml
        )
        @test published.artifact == "gauntlet_fixture"
        @test release_tag(:fixture, v"1.0.0") == "gauntlet-fixture-v1.0.0"
        @test collection_archive_name(:fixture, v"1.0.0") ==
              "benchmarks-fixture-v1.0.0.tar.gz"
        @test_throws ArgumentError artifact_name(:Fixture)
        @test_throws ArgumentError release_tag(:fixture, v"0.1.0")
        @test occursin("file://", read(artifacts_toml, String))
        @test isfile(snapshot_path(case; artifacts_toml))

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

        new_source=joinpath(directory, "new_case.jl")
        write(new_source, "# new case with no accepted baseline\n")
        new_case=fixture_case(new_source; name = :new_fixture)
        @test load_prior_snapshot(new_case; artifacts_toml) === nothing

        outcome=run_snapshot(
            case;
            artifacts_toml,
            options = (
                output_basis = :pul,
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
        @test_throws ArgumentError GauntletSupport.report(:fixture; artifact_root)

        @test_throws ArgumentError prepare_staging(; artifact_root)
        prepare_staging(; artifact_root, force = true)
        @test !ispath(collection_stage(:fixture; artifact_root))
        @test isfile(artifacts_toml)
    end

    previous=get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", nothing)
    previous_ci=get(ENV, "CI", nothing)
    previous_persist=get(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST", nothing)
    previous_runner=get(ENV, "LINECABLEMODELS_GAUNTLET_RUNNER", nothing)
    previous_cleanup=get(ENV, "LINECABLEMODELS_GAUNTLET_CLEANUP", nothing)
    previous_stage_force=get(
        ENV,
        "LINECABLEMODELS_GAUNTLET_STAGE_FORCE",
        nothing
    )
    try
        delete!(ENV, "CI")
        for mode in ("snapshot", "live", "record")
            ENV["LINECABLEMODELS_GAUNTLET_MODE"]=mode
            @test gauntlet_mode() == Symbol(mode)
        end
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="automatic"
        @test_throws ArgumentError gauntlet_mode()
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_CLEANUP")
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_STAGE_FORCE")
        @test !gauntlet_cleanup()
        @test !gauntlet_stage_force()
        ENV["LINECABLEMODELS_GAUNTLET_CLEANUP"]="true"
        ENV["LINECABLEMODELS_GAUNTLET_STAGE_FORCE"]="TRUE"
        @test gauntlet_cleanup()
        @test gauntlet_stage_force()
        ENV["LINECABLEMODELS_GAUNTLET_STAGE_FORCE"]="yes"
        @test_throws ArgumentError gauntlet_stage_force()
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
        previous_stage_force===nothing ?
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_STAGE_FORCE") :
        (ENV["LINECABLEMODELS_GAUNTLET_STAGE_FORCE"]=previous_stage_force)
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
            system = (
                system_id = "ports",
                connection_order = collect(assignments)
            ),
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

@testitem "Gauntlet / case catalog and benchmark-file invariants" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using SHA
    using LineCableModels
    using .GauntletSupport

    index=case_index()
    expected_sizes=Dict(
        :cable_18kv_1000mm2_trefoil=>(9, 9, 101),
        :cable_132kv_630mm2_flathor=>(9, 9, 101),
        :cable_380kv_2000mm2_flatver=>(9, 9, 101),
        :cable_525kv_1600mm2_bipole=>(6, 6, 101),
        :cable_640kv_2000mm2_bipole=>(6, 6, 101),
        :solid_1000mm2_single=>(1, 1, 101),
        :two_bare_wires=>(2, 2, 101)
    )
    @test Set(keys(index)) == Set(keys(expected_sizes))
    case_files=Set(realpath(path)
    for path in filter(
        path->endswith(path, ".jl"),
        readdir(GauntletSupport.CASE_ROOT; join = true)
    ))
    @test Set(values(index)) == case_files

    mktempdir() do directory
        root=joinpath(directory, "cases")
        mkpath(root)
        outside=joinpath(directory, "outside.jl")
        write(outside, "nothing\n")

        escaping_index=joinpath(root, "escaping.toml")
        write(escaping_index, "[cases]\nescape = \"../outside.jl\"\n")
        @test_throws ArgumentError case_index(
            index_path = escaping_index,
            case_root = root
        )

        shared=joinpath(root, "shared.jl")
        write(shared, "nothing\n")
        duplicate_index=joinpath(root, "duplicate.toml")
        write(
            duplicate_index,
            "[cases]\nfirst = \"shared.jl\"\nsecond = \"shared.jl\"\n"
        )
        @test_throws ArgumentError case_index(
            index_path = duplicate_index,
            case_root = root
        )

        mismatched=joinpath(root, "requested.jl")
        write(
            mismatched,
            "case_definition(\n"*
            "    :different,\n"*
            "    (frequency = case_parameter(:frequency, [50.0]; "*
            "tags = (:operation, :frequency)),),\n"*
            "    [\"terminal\"]\n"*
            ") do _\n"*
            "    error(\"an ID mismatch must fail before construction\")\n"*
            "end\n"
        )
        mismatch_index=joinpath(root, "mismatch.toml")
        write(mismatch_index, "[cases]\nrequested = \"requested.jl\"\n")
        @test_throws ArgumentError load_case(
            :requested;
            index_path = mismatch_index,
            case_root = root
        )
    end

    expected_wire_counts=Dict(
        :cable_132kv_630mm2_flathor=>(
            (1, 6, 12, 18, 24), (19,), ()
        ),
        :cable_18kv_1000mm2_trefoil=>(
            (1, 6, 12, 18, 24), (49,), ()
        ),
        :cable_380kv_2000mm2_flatver=>(
            (1, 6, 12, 18, 24), (14,), ()
        ),
        :cable_525kv_1600mm2_bipole=>(
            (1, 6, 12, 18, 24, 30, 36), (), (68,)
        ),
        :cable_640kv_2000mm2_bipole=>(
            (6, 12, 18, 24), (), ()
        ),
        :solid_1000mm2_single=>((),),
        :two_bare_wires=>((),)
    )
    for (id, expected_size) in expected_sizes
        model=load_case(id)
        @test model.id === id
        @test model.expected_size == expected_size
        @test model.source_file == index[id]
        @test model.source_sha256 == bytes2hex(SHA.sha256(read(index[id])))
        source=read(index[id], String)
        @test length(collect(eachmatch(r"\bcase_definition\s*\(", source))) == 1
        @test !occursin("@testitem", source)
        @test !occursin(r"\bFormulation\s*\(", source)
        @test !occursin(
            r"\b(?:maxfill|make_screened|make_stranded)\s*\(",
            source
        )
        for parameter in values(model.definition.parameters)
            :topology in parameter.tags||continue
            topology=parameter.nominal
            @test topology isa Integer ||
                  (topology isa Tuple && all(value -> value isa Integer, topology))
        end
        design=first(model.nominal_problem.system.designs)
        groups=filter(item->item isa Group, design.root.items)
        component_wire_counts=Tuple(
            Tuple(
                group.pattern.n
            for group in groups
            if group.name===terminal&&group.pattern isa Ring
            ) for terminal in design.terminal_order
        )
        @test component_wire_counts == expected_wire_counts[id]
        @test basename(index[id]) == string(id, ".jl")
    end

    benchmark_root=joinpath(GauntletSupport.GAUNTLET_ROOT, "benchmarks")
    benchmark_files=sort!(filter(
        path->endswith(path, ".jl")&&
        !startswith(path, joinpath(benchmark_root, ".work")*Base.Filesystem.path_separator),
        collect(Iterators.flatten(
            (joinpath(root, file) for file in files)
        for (root, _, files) in walkdir(benchmark_root)
        ))
    ))
    contains_for(value) = value isa Expr&&
    (value.head===:for||any(contains_for, value.args))
    benchmark_ids=Symbol[]
    benchmark_cases=Dict(:pscad=>Symbol[], :uq=>Symbol[])
    for path in benchmark_files
        source=read(path, String)
        parsed=Meta.parseall(source)
        expressions=parsed.head===:toplevel ?
                    filter(expression->!(expression isa LineNumberNode), parsed.args) :
                    Any[parsed]
        testitems=filter(expressions) do expression
            expression isa Expr&&expression.head===:macrocall&&
                expression.args[1]===Symbol("@testitem")
        end
        @test length(testitems) == 1
        @test length(collect(eachmatch(r"\bload_case\s*\(", source))) == 1
        case_matches=collect(eachmatch(
            r"\bload_case\s*\(\s*:([a-z][a-z0-9_]*)"s,
            source
        ))
        @test length(case_matches) == 1
        if length(case_matches)==1
            case_id=Symbol(only(only(case_matches).captures))
            @test case_id in keys(index)
            collection=Symbol(basename(dirname(path)))
            haskey(benchmark_cases, collection)&&
            push!(benchmark_cases[collection], case_id)
        end
        @test !contains_for(only(testitems))
        benchmark_matches=collect(eachmatch(
            r"(?:GauntletCase|benchmark_definition)\s*\(\s*:(benchmark_[a-z0-9_]+)"s,
            source
        ))
        @test length(benchmark_matches) == 1
        length(benchmark_matches)==1||continue
        benchmark_id=Symbol(only(only(benchmark_matches).captures))
        push!(benchmark_ids, benchmark_id)
        @test splitext(basename(path))[1] == string(benchmark_id)
    end
    @test length(unique(benchmark_ids)) == length(benchmark_ids)
    for collection in (:pscad, :uq)
        @test sort(benchmark_cases[collection]; by = string) ==
              sort!(collect(keys(index)); by = string)
    end
end

@testitem "Gauntlet / case variations and correlation" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using Measurements
    using LineCableModels
    using .GauntletSupport

    first_load=load_case(:cable_18kv_1000mm2_trefoil)
    first_load.nominal_problem.system.connections[1][1]=99
    second_load=load_case(:cable_18kv_1000mm2_trefoil)
    @test second_load.nominal_problem.system.connections[1][1] == 1

    single_frequency=load_case(
        :cable_18kv_1000mm2_trefoil;
        variation = ExactOverrides(frequencies = [50.0])
    )
    @test single_frequency.expected_size == (9, 9, 1)
    @test single_frequency.problem.frequencies == [50.0]
    @test single_frequency.selected_parameters == [:frequencies]

    frequency_grid=load_case(
        :two_bare_wires;
        variation = ParameterGrids(frequencies = Grid(([50.0], [60.0])))
    )
    @test length(frequency_grid.problem) == 2
    @test [only(problem.frequencies) for problem in frequency_grid.problem] == [50.0, 60.0]

    @test_throws ArgumentError load_case(
        :two_bare_wires;
        variation = ExactOverrides(unknown_parameter = 1.0)
    )
    @test_throws ArgumentError load_case(
        :two_bare_wires;
        variation = RelativeStandardUncertainty(10.0; tags = (:not_a_tag,))
    )
    @test_throws ArgumentError load_case(
        :two_bare_wires;
        variation = ExactOverrides(frequencies = Grid(([50.0],)))
    )
    @test_throws ArgumentError load_case(:unknown_case)

    uncertain=load_case(
        :cable_18kv_1000mm2_trefoil;
        variation = RelativeStandardUncertainty(
            10.0; tags = (:geometry, :cable_layer)
        )
    )
    expected_selected=[
        :core_strand_diameter,
        :core_ring_lay_ratio_1,
        :core_ring_lay_ratio_2,
        :core_ring_lay_ratio_3,
        :core_ring_lay_ratio_4,
        :screen_wire_diameter,
        :screen_wire_lay_ratio,
        :semicon_tape_thickness,
        :inner_semicon_thickness,
        :insulation_thickness,
        :outer_semicon_thickness,
        :copper_tape_thickness,
        :copper_tape_width,
        :copper_tape_lay_ratio,
        :water_blocking_thickness,
        :aluminum_tape_thickness,
        :pe_face_thickness,
        :jacket_thickness
    ]
    @test uncertain.selected_parameters == expected_selected
    @test length(uncertain.problem) == 1
    for id in expected_selected
        descriptor=only(getproperty(uncertain.sources, id))
        nominal=getproperty(uncertain.definition.parameters, id).nominal
        @test LineCableModels.nominal(descriptor) == nominal
        @test LineCableModels.uncertainty(descriptor) ≈ 0.1abs(nominal)
    end
    @test uncertain.sources.formation_clearance_ratio == 0.2
    @test uncertain.sources.screen_wires == 49
    @test uncertain.sources.earth_rho == 100.0

    nominal_system=uncertain.nominal_problem.system
    nominal_radius=outer_radius(first(nominal_system.designs))
    nominal_spacing=hypot(
        nominal_system.positions[1].x-nominal_system.positions[2].x,
        nominal_system.positions[1].y-nominal_system.positions[2].y
    )
    @test nominal_radius ≈ 34.575e-3
    @test nominal_spacing ≈ 2.2nominal_radius
    @test nominal_spacing - 2nominal_radius ≈ 0.2nominal_radius

    expanded_armor=load_case(
        :cable_525kv_1600mm2_bipole;
        variation = ExactOverrides(armor_wire_diameter = 1.2*5.827e-3)
    )
    expanded_design=first(expanded_armor.problem.system.designs)
    armor_group=only(filter(
        item->item isa Group&&item.name===:armor,
        expanded_design.root.items
    ))
    armor=only(filter(
        value->value.name===:armor,
        LineCableModels.Engine.analytical_components(
            LineCableModels.Formulation(), expanded_design, 50.0
        )
    )).conductor
    armor_wire_radius=armor_group.item.primitive.r
    @test armor.num_wires == 68
    @test armor.r_in >=
          armor_wire_radius / sinpi(1 / armor.num_wires) - armor_wire_radius
    expanded_bedding=only(filter(
        source->source.source.tag===:sheath_bedding,
        expanded_design.geometry.regions
    ))
    @test thickness(expanded_bedding.shape) > 3.0e-3

    nominal_armor=load_case(:cable_525kv_1600mm2_bipole)
    nominal_design=first(nominal_armor.problem.system.designs)
    nominal_bedding=only(filter(
        source->source.source.tag===:sheath_bedding,
        nominal_design.geometry.regions
    ))
    nominal_unbuffered_outer_radius=63.22185e-3+5.827e-3+10.0e-3
    @test thickness(nominal_bedding.shape) ≈
          3.0e-3 + 0.2nominal_unbuffered_outer_radius

    materialized=only(uncertain.problem)
    uncertain_system=materialized.system
    uncertain_radius=outer_radius(first(uncertain_system.designs))
    uncertain_spacing=hypot(
        uncertain_system.positions[1].x-uncertain_system.positions[2].x,
        uncertain_system.positions[1].y-uncertain_system.positions[2].y
    )
    @test uncertain_spacing - 2uncertain_radius ≈ 0.2uncertain_radius

    geometry=first(materialized.system.designs).geometry.regions
    first_tape=thickness(only(filter(
        source->source.source.tag===:core_semicon_tape_inner,
        geometry
    )).shape)
    second_tape=thickness(only(filter(
        source->source.source.tag===:core_semicon_tape_outer,
        geometry
    )).shape)
    @test Measurements.cov(first_tape, second_tape) > 0
    @test Measurements.cov(first_tape, second_tape) ≈
          Measurements.uncertainty(first_tape)^2

    uq_model=load_case(
        :cable_132kv_630mm2_flathor;
        variation = RelativeStandardUncertainty(
            10.0; tags = (:geometry, :cable_layer)
        )
    )
    uq_selected=[
        :core_strand_diameter,
        :core_lay_ratio,
        :screen_wire_diameter,
        :screen_lay_ratio,
        :semicon_tape_thickness,
        :inner_semicon_thickness,
        :insulation_thickness,
        :outer_semicon_thickness,
        :copper_tape_thickness,
        :copper_tape_width,
        :copper_tape_lay_ratio,
        :water_blocking_thickness,
        :aluminum_tape_thickness,
        :jacket_thickness
    ]
    @test uq_model.selected_parameters == uq_selected
    @test length(uq_model.problem) == 1
    for parameter in values(uq_model.definition.parameters)
        source=getproperty(uq_model.sources, parameter.id)
        if parameter.id in uq_selected
            descriptor=only(source)
            @test LineCableModels.nominal(descriptor) == parameter.nominal
            @test LineCableModels.uncertainty(descriptor) ≈
                  0.1abs(parameter.nominal)
        else
            @test isequal(source, parameter.nominal)
        end
    end
end

@testitem "Gauntlet / UQ moment comparison contract" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport

    frequencies_value=[1.0, 10.0]
    ports=["a", "b"]
    values=map((R = 1.0, L = 2.0, C = 3.0, G = 4.0)) do scale
        (
            mean = fill(scale, 2, 2, 2),
            std = fill(scale/10, 2, 2, 2)
        )
    end
    reference=MomentResult(values, frequencies_value, :pul, PhaseDomain, ports)
    equal_comparison=compare(reference, reference)
    tolerance=(
        mean = map(_->(absolute = 0.0, relative = 0.0), values),
        std = map(_->(absolute = 0.0, relative = 0.0), values)
    )
    @test reference ≈ reference
    @test moment_comparison_passes(equal_comparison, tolerance)
    @test all(iszero, equal_comparison.errors.R.mean.absolute)

    changed_values=merge(values, (
        R = (mean = copy(values.R.mean), std = copy(values.R.std)),
    ))
    changed_values.R.mean[1, 2, :].=1.5
    changed=MomentResult(changed_values, frequencies_value, :pul, PhaseDomain, ports)
    changed_comparison=compare(reference, changed)
    @test changed_comparison.errors.R.mean.absolute[1, 2] == 0.5
    @test !moment_comparison_passes(changed_comparison, tolerance)

    zero_values=map(values) do product
        (mean = zeros(size(product.mean)), std = zeros(size(product.std)))
    end
    small_values=map(zero_values) do product
        (mean = fill(1.0e-12, size(product.mean)), std = copy(product.std))
    end
    zero_reference=MomentResult(zero_values, frequencies_value, :pul, PhaseDomain, ports)
    small_candidate=MomentResult(small_values, frequencies_value, :pul, PhaseDomain, ports)
    floor_tolerance=(
        mean = map(_->(absolute = 1.0e-11, relative = 0.0), values),
        std = map(_->(absolute = 0.0, relative = 0.0), values)
    )
    @test all(isinf, compare(zero_reference, small_candidate).errors.R.mean.relative)
    @test moment_comparison_passes(
        compare(zero_reference, small_candidate), floor_tolerance
    )

    @test_throws ArgumentError compare(
        reference,
        MomentResult(values, [1.0, 11.0], :pul, PhaseDomain, ports)
    )
    @test_throws ArgumentError compare(
        reference,
        MomentResult(values, frequencies_value, :total, PhaseDomain, ports)
    )
    @test_throws ArgumentError compare(
        reference,
        MomentResult(values, frequencies_value, :pul, PhaseDomain, reverse(ports))
    )
    wrong_shape=merge(values, (
        R = (mean = zeros(1, 1, 2), std = zeros(1, 1, 2)),
    ))
    @test_throws DimensionMismatch compare(
        reference,
        MomentResult(wrong_shape, frequencies_value, :pul, PhaseDomain, ports)
    )
end

@testitem "Gauntlet / fixed-seed Monte Carlo reproducibility" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport

    model=load_case(
        :cable_132kv_630mm2_flathor;
        variation = compose_variations(
            ExactOverrides(frequencies = [50.0]),
            RelativeStandardUncertainty(
                10.0; tags = (:geometry, :cable_layer)
            )
        )
    )
    inner=Formulation(
        earth_impedance = EarthImpedance.Pollaczek(),
        earth_admittance = EarthAdmittance.IdealGround(),
        insulation_admittance = InsulationAdmittance.Lossless(),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
    formulation=MonteCarlo(
        inner;
        trials = 8,
        seed = 0x51eed,
        distribution = :normal,
        return_samples = false,
        return_histograms = false
    )
    problem=ParametricProblem(model.problem)
    first_result=compute(problem, formulation)
    second_result=compute(problem, formulation)
    @test first_result.root_seed == second_result.root_seed == UInt64(0x51eed)
    @test first_result.trial_counts == second_result.trial_counts == [8]
    @test statistics(first_result) == statistics(second_result)
    @test samples(first_result) === nothing
    @test histograms(first_result) === nothing
end

@testitem "Gauntlet / owned formula comparison dispatch" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Test
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport

    previous_mode=get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", nothing)
    try
        ENV["LINECABLEMODELS_GAUNTLET_MODE"]="live"
        model=load_case(
            :two_bare_wires;
            variation = ExactOverrides(frequencies = [50.0])
        )
        pollaczek=Formulation(
            earth_impedance = EarthImpedance.Pollaczek(),
            earth_admittance = EarthAdmittance.IdealGround(),
            insulation_admittance = InsulationAdmittance.Lossless(),
            options = (kron_reduction = false, reduce_bundle = false)
        )
        papadopoulos=Formulation(
            earth_impedance = EarthImpedance.Papadopoulos(),
            earth_admittance = EarthAdmittance.IdealGround(),
            insulation_admittance = InsulationAdmittance.Lossless(),
            options = (kron_reduction = false, reduce_bundle = false)
        )
        limits=(
            Z = (absolute = 1.0, relative = 1.0),
            Y = (absolute = 1.0, relative = 1.0)
        )
        benchmark=benchmark_definition(
            :benchmark_owned_formula_fixture,
            model.id,
            :owned,
            @__FILE__,
            model,
            benchmark_calculation(
                :pollaczek, :engine, model.problem, pollaczek
            ),
            benchmark_calculation(
                :papadopoulos, :engine, model.problem, papadopoulos
            ),
            LineParametersPolicy(),
            (reference = limits, regression = limits)
        )
        outcome=run_benchmark(benchmark)
        @test outcome.reference isa LineParameters
        @test outcome.candidate isa LineParameters
        @test outcome.passes
        @test outcome.metadata.calculations.reference.id === :pollaczek
        @test outcome.metadata.calculations.candidate.id === :papadopoulos
        @test size(Z(outcome.reference)) == (2, 2, 1)
    finally
        previous_mode===nothing ?
        delete!(ENV, "LINECABLEMODELS_GAUNTLET_MODE") :
        (ENV["LINECABLEMODELS_GAUNTLET_MODE"]=previous_mode)
    end
end

@testitem "Gauntlet / PSCAD parser and fixed mappings" tags=[:gauntlet_toolkit] setup=[
    GauntletSupport
] begin
    using Base64: base64decode
    using Test
    import LineCableModels
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
        raw"Z:\gauntlet\benchmarks\.work",
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
          raw"Z:\gauntlet\benchmarks\.work\pscad\benchmark_525kV_1600mm2_bipole_pscad\reference"
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
    @test occursin("-EarthField 'EarthForm2'", supervisor_command)
    @test occursin("-EarthValue '0'", supervisor_command)
    @test occursin("-EarthReadback 'DERISEMLYEN'", supervisor_command)
    @test !occursin("OpenStandardInput", supervisor_command)
    saad_command=harness._supervisor_command(
        config,
        raw"Z:\gauntlet\benchmarks\.work\case\current",
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
