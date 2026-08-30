@testitem "Gmsh FEM / public API, strict parsing, and UI transitions" tags=[:extension] begin
    import LineCableModels
    using Gmsh
    using LinearAlgebra

    extension_module = Base.get_extension(
        LineCableModels, :LineCableModelsGmshExt
    )
    @test extension_module !== nothing
    @test LineCableModels.LineCableModelsFEM <:
          LineCableModels.Engine.AbstractFormulationBackend
    @test LineCableModels.LineCableModelsFEMOptions <:
          LineCableModels.Engine.AbstractFormulationOptions

    options = LineCableModels.LineCableModelsFEMOptions(
        mesh_policy = :remesh,
        keep_run_directory = true,
        gmsh_verbosity = 0,
        getdp_verbosity = 5
    )
    @test options.mesh_policy === :remesh
    @test options.keep_run_directory
    @test options.getdp_verbosity == 5
    @test_throws ArgumentError LineCableModels.LineCableModelsFEMOptions(
        mesh_policy = :invalid
    )
    @test_throws ArgumentError LineCableModels.LineCableModelsFEMOptions(
        getdp_verbosity = 6
    )

    formulation = LineCableModels.Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,),
        fem_options = options
    )
    @test formulation isa LineCableModels.LineCableModelsFEM
    @test !formulation.options.ideal_transposition
    @test formulation.execution === options

    transition = extension_module._ui_transition
    @test transition(:geometry_ready, :run_model) === :mesh_required
    @test transition(:geometry_ready, :generate_mesh) === :mesh_requested
    @test transition(:mesh_ready, :run_model) === :solve_requested
    @test transition(:geometry_ready, Symbol("")) === :geometry_ready
    @test transition(:geometry_ready, Symbol(""), false) === :closed_before_mesh
    @test transition(:mesh_ready, Symbol(""), false) === :closed_before_solve
    @test extension_module._ui_state_labels(:geometry_ready) == (
        mesh_state = "not generated", solve_state = "not run"
    )
    @test extension_module._ui_state_labels(:mesh_ready) == (
        mesh_state = "ready", solve_state = "not run"
    )
    @test extension_module._ui_state_labels(:results_ready) == (
        mesh_state = "ready", solve_state = "completed"
    )

    mktempdir() do directory
        run = extension_module.FEMRun(
            directory, extension_module.created, "fixture", :none, ""
        )
        path = joinpath(directory, "Z.tsv")
        header = join(extension_module.FEM_RAW_HEADER, '\t')
        rows = [
            "1\t50\t1\t1\t1\t2",
            "1\t50\t2\t1\t3\t4",
            "1\t50\t1\t2\t5\t6",
            "1\t50\t2\t2\t7\t8"
        ]
        open(path, "w") do io
            println(io, header)
            println.(Ref(io), rows)
        end
        matrix = extension_module._parse_raw_matrix(
            Float64, path, [50.0], 2, run
        )
        @test size(matrix) == (2, 2, 1)
        @test matrix[1, 1, 1] == 1 + 2im
        @test matrix[2, 1, 1] == 3 + 4im
        @test matrix[1, 2, 1] == 5 + 6im
        @test matrix[2, 2, 1] == 7 + 8im

        open(path, "w") do io
            println(io, header)
            println.(Ref(io), [rows[1], rows[1], rows[3], rows[4]])
        end
        @test_throws LineCableModels.LineCableModelsFEMError extension_module._parse_raw_matrix(
            Float64, path, [50.0], 2, run
        )

        extension_module._transition!(run, extension_module.running, "fixture running")
        @test run.state === extension_module.running
        @test occursin("\"state\": \"running\"", read(
            joinpath(directory, "run.json"), String
        ))
        missing_completion = try
            extension_module._validate_completion(
                joinpath(directory, "missing.tsv"), 1, 2, run
            )
            nothing
        catch exception
            exception
        end
        @test missing_completion isa LineCableModels.LineCableModelsFEMError
        @test missing_completion.run_directory == directory

        open(path, "w") do io
            println(io, header)
            println.(Ref(io), ["1\t50\t1\t1\tNaN\t0", rows[2], rows[3], rows[4]])
        end
        @test_throws LineCableModels.LineCableModelsFEMError extension_module._parse_raw_matrix(
            Float64, path, [50.0], 2, run
        )
    end

    potential = Array{ComplexF64, 3}(undef, 2, 2, 2)
    potential[:, :, 1] = [2.0 0.25; 0.25 1.5]
    potential[:, :, 2] = [3.0 + 0.1im 0.5; 0.5 2.0 + 0.2im]
    inversion = LineCableModels.Engine.potential_to_admittance(
        potential; diagnostics = true
    )
    for frequency in axes(potential, 3)
        @test potential[:, :, frequency] * inversion.Y[:, :, frequency] ≈ I
        @test inversion.residuals[frequency] ≤ 100eps(Float64)
    end
end

@testitem "Gmsh FEM / deterministic geometry and mesh lifecycle" tags=[:extension] begin
    using LineCableModels
    using Gmsh

    copper = Material(kind = :conductor, rho = 1 / 5.8e7)
    dielectric = Material(kind = :insulator, rho = 1.0e15, eps_r = 2.3)
    design = build(CableDesign,
        "fem-mesh",
        Stack(
            Group(:core, Region(:core_metal, Disk(0.005), copper)),
            Region(:dielectric, Annulus(0.005, 0.01), dielectric),
            Group(:sheath, Region(:sheath_metal, Annulus(0.01, 0.011), copper)),
            Region(:jacket, Annulus(0.011, 0.013), dielectric)
        ))
    system = build(
        LineCableSystem,
        design,
        (0.0, -0.1);
        connections = Dict(:core => 1, :sheath => 0),
        system_id = "fem-mesh",
        line_length = 1.0
    )
    problem = LineParametersProblem(
        system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0]
    )
    formulation = Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,),
        fem_options = (gmsh_verbosity = 0,)
    )
    extension_module = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    model = extension_module._resolved_fem_model(problem, formulation)
    @test model.terminal_ids == [
        "cable_0001/fem-mesh/core",
        "cable_0001/fem-mesh/sheath"
    ]
    @test getproperty.(model.material_plans, :physical_tag) == 10_001:10_004
    @test model.tags.terminal_base == 3_000
    @test model.tags.infinite_domain == 1_005
    @test model.shell_outer_radius == 1.25model.domain_radius

    gamma_problem = LineParametersProblem(
        system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0],
        Γ = [1.0e-12im]
    )
    gamma_error = try
        extension_module._resolved_fem_model(gamma_problem, formulation)
        nothing
    catch exception
        exception
    end
    @test gamma_error isa LineCableModelsFEMError
    @test gamma_error.category === :unsupported
    @test gamma_error.field === :Γ

    vertical_problem = LineParametersProblem(
        system;
        earth_props = LineCableModels.EarthProps.EarthModel(
            100.0, 10.0, 1.0; vertical_layers = true
        ),
        frequencies = [50.0]
    )
    vertical_error = try
        extension_module._resolved_fem_model(vertical_problem, formulation)
        nothing
    catch exception
        exception
    end
    @test vertical_error isa LineCableModelsFEMError
    @test vertical_error.category === :unsupported
    @test vertical_error.field === :vertical_layers

    Gmsh.initialize(String[]; finalize_atexit = false)
    gmsh.model.add("caller-owned")
    caller_view = gmsh.view.add("caller-view")
    gmsh.option.set_number("General.Terminal", 1)
    gmsh.option.set_number("General.Verbosity", 4)
    caller_parameter = extension_module._onelab_name("caller_parameter")
    gmsh.onelab.set_string(caller_parameter, ["preserve-me"])
    caller_models = Set(String.(gmsh.model.list()))
    try
        mktempdir() do runtime_root
            session = extension_module._start_gmsh(0)
            try
                @test !session.owned
                gmsh.view.add("temporary-view")
                gmsh.model.add("fem-exact-curves")
                registry = extension_module.FEMLoopRegistry(1.0e-3)
                ellipse_loop = extension_module._boundary_loop!(
                    registry,
                    Ellipse(0.02, 0.01, Pose2(-0.05, 0.0, π / 7))
                )
                sector_loop = extension_module._boundary_loop!(
                    registry,
                    Sector(0.01, 0.02, -π / 4, π / 2, Pose2(0.0, 0.0, π / 9))
                )
                rounded_loop = extension_module._boundary_loop!(
                    registry,
                    LineCableModels.DataModel.RoundedSectorShape(
                        RoundedSector(π / 3, 0.005, 0.02, 0.001),
                        Pose2(0.05, 0.0, -π / 11)
                    )
                )
                gmsh.model.geo.synchronize()
                @test length(ellipse_loop.curves) == 4
                @test length(sector_loop.curves) == 4
                @test all(
                    gmsh.model.get_type(1, curve) == "Ellipse"
                for curve in ellipse_loop.curves
                )
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Circle",
                    sector_loop.curves
                ) == 2
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Circle",
                    rounded_loop.curves
                ) == 7
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Line",
                    rounded_loop.curves
                ) == 2

                first_run = extension_module._create_run(runtime_root)
                first_geometry = extension_module._build_geometry!(model, "fem-mesh-first")
                first_mesh = extension_module._select_mesh!(
                    first_run, model, first_geometry, formulation, runtime_root
                )
                @test isfile(first_mesh)
                @test first_run.mesh_source === :generated
                @test !isempty(first_run.mesh_fingerprint)
                @test length(first_geometry.infinite_surfaces) == 2
                @test length(first_geometry.inner_shell_curves) == 4
                extension_module._validate_mesh_file(model, first_mesh)

                second_run = extension_module._create_run(runtime_root)
                second_geometry = extension_module._build_geometry!(model, "fem-mesh-second")
                second_mesh = extension_module._select_mesh!(
                    second_run, model, second_geometry, formulation, runtime_root
                )
                @test isfile(second_mesh)
                @test second_run.mesh_source === :cache
                @test second_run.mesh_fingerprint == first_run.mesh_fingerprint

                remesh = Formulation(
                    :LineCableModelsFEM;
                    options = (ideal_transposition = false,),
                    fem_options = (mesh_policy = :remesh, gmsh_verbosity = 0)
                )
                third_run = extension_module._create_run(runtime_root)
                third_geometry = extension_module._build_geometry!(model, "fem-mesh-third")
                third_mesh = extension_module._select_mesh!(
                    third_run, model, third_geometry, remesh, runtime_root
                )
                @test isfile(third_mesh)
                @test third_run.mesh_source === :generated

                if Sys.isunix()
                    failing_getdp = joinpath(runtime_root, "getdp-failure")
                    open(failing_getdp, "w") do io
                        println(io, "#!/bin/sh")
                        println(io, "if [ \"\$1\" = \"-info\" ]; then")
                        println(io, "  echo 'GetDP Version 3.5.0'")
                        println(io, "  exit 0")
                        println(io, "fi")
                        println(io, "echo 'intentional GetDP failure marker' >&2")
                        println(io, "exit 17")
                    end
                    chmod(failing_getdp, 0o700)
                    failing_formulation = Formulation(
                        :LineCableModelsFEM;
                        options = (ideal_transposition = false,),
                        fem_options = (
                            getdp_executable = failing_getdp,
                            gmsh_verbosity = 0,
                            getdp_verbosity = 0
                        )
                    )
                    failure_run = extension_module._create_run(runtime_root)
                    extension_module._prepare_run_inputs!(failure_run, model)
                    material_functions = read(
                        joinpath(
                            failure_run.path,
                            "input",
                            "material_functions.pro"
                        ), String)
                    @test occursin("sigma_dc[MaterialRegion1]", material_functions)
                    @test occursin("tan_delta[MaterialRegion1]", material_functions)
                    failure = try
                        extension_module._run_getdp!(
                            failure_run,
                            model,
                            failing_formulation,
                            first_mesh
                        )
                        nothing
                    catch exception
                        exception
                    end
                    @test failure isa LineCableModelsFEMError
                    @test failure.category === :getdp
                    @test failure.field === :client
                    @test occursin("Gmsh log tail", failure.message)
                    @test failure.run_directory == failure_run.path
                    @test failure_run.getdp_invocations == 1
                    @test isdir(failure_run.path)
                    @test isfile(joinpath(failure_run.path, "logs", "getdp.log"))
                end
            finally
                extension_module._finish_gmsh(session)
            end
        end
        @test Bool(gmsh.is_initialized())
        @test Set(String.(gmsh.model.list())) == caller_models
        @test gmsh.model.get_current() == "caller-owned"
        @test gmsh.view.get_tags() == [caller_view]
        @test gmsh.option.get_number("General.Terminal") == 1
        @test gmsh.option.get_number("General.Verbosity") == 4
        @test gmsh.onelab.get_string(caller_parameter) == ["preserve-me"]
    finally
        Gmsh.finalize()
    end
end

@testitem "Gmsh FEM / optional real GetDP multi-frequency scan" tags=[
    :extension,
    :integration
] begin
    using LineCableModels
    using Gmsh

    executable = get(ENV, "LINECABLEMODELS_GETDP", something(Sys.which("getdp"), ""))
    if isempty(executable)
        @test true
    else
        copper = Material(kind = :conductor, rho = 1 / 5.8e7)
        dielectric = Material(kind = :insulator, rho = 1.0e15, eps_r = 2.3)
        design = build(CableDesign,
            "fem-integration",
            Stack(
                Group(:core, Region(:core_metal, Disk(0.005), copper)),
                Region(:dielectric, Annulus(0.005, 0.01), dielectric),
                Group(:sheath, Region(:sheath_metal, Annulus(0.01, 0.011), copper)),
                Region(:jacket, Annulus(0.011, 0.013), dielectric)
            ))
        system = build(
            LineCableSystem,
            design,
            (0.0, -0.1);
            connections = Dict(:core => 1, :sheath => 0),
            system_id = "fem-integration",
            line_length = 1.0
        )
        problem = LineParametersProblem(
            system;
            earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
            frequencies = [50.0, 1000.0]
        )
        formulation = Formulation(
            :LineCableModelsFEM;
            options = (ideal_transposition = false,),
            fem_options = (
                getdp_executable = executable,
                getdp_verbosity = 0,
                gmsh_verbosity = 0
            )
        )
        result = compute(problem, formulation; options = (trace = true,))
        @test size(result.Z) == (1, 1, 2)
        @test size(result.Y) == (1, 1, 2)
        @test all(isfinite, result.Z)
        @test all(isfinite, result.Y)
        @test result.f == [50.0, 1000.0]
        @test result.details.fem.run.state === Base.get_extension(
            LineCableModels, :LineCableModelsGmshExt
        ).completed
        @test result.details.fem.run.getdp_invocations == 1
        @test result.details.fem.run.run_directory === nothing

        if !isempty(get(ENV, "DISPLAY", ""))
            ui_formulation = Formulation(
                :LineCableModelsFEM;
                options = (ideal_transposition = false,),
                fem_options = (
                    ui = true,
                    getdp_executable = executable,
                    getdp_verbosity = 0,
                    gmsh_verbosity = 0
                )
            )
            action_name = "LineCableModels/FEM/ui/action"
            mesh_state_name = "LineCableModels/FEM/ui/mesh_state"
            solve_state_name = "LineCableModels/FEM/ui/solve_state"
            driver = @async begin
                deadline = time() + 60
                while !Bool(gmsh.is_initialized()) || isempty(try
                    gmsh.onelab.get_names("^$action_name\$")
                catch
                    String[]
                end)
                    time() < deadline || error("Gmsh UI action did not become ready")
                    sleep(0.02)
                end
                gmsh.onelab.set_string(action_name, ["generate_mesh"])
                while gmsh.onelab.get_string(mesh_state_name) != ["ready"]
                    time() < deadline || error("Gmsh UI mesh action timed out")
                    sleep(0.02)
                end
                gmsh.onelab.set_string(action_name, ["run_model"])
                while gmsh.onelab.get_string(solve_state_name) != ["completed"]
                    time() < deadline || error("Gmsh UI solve action timed out")
                    sleep(0.02)
                end
                gmsh.fltk.finalize()
            end
            ui_result = compute(problem, ui_formulation; options = (trace = true,))
            wait(driver)
            @test ui_result.f == [50.0, 1000.0]
            @test ui_result.details.fem.run.getdp_invocations == 1
            @test ui_result.details.fem.run.run_directory === nothing
        end
    end
end

@testitem "Gmsh FEM / frozen Python reference fixture contract" tags=[
    :extension
] begin
    using JSON3
    using LineCableModels

    path = joinpath(
        pkgdir(LineCableModels),
        "test",
        "fixtures",
        "reference",
        "fem_python_quasi_tem.json"
    )
    fixture = JSON3.read(read(path, String))
    @test fixture.schema == "LineCableModels.FEMPythonReference"
    @test fixture.schema_version == 1
    @test fixture.cases.single_coaxial_kron.frequencies_hz == [
        10.0, 1000.0, 100000.0
    ]
    @test length(fixture.cases.single_coaxial_kron.terminal_order) == 2
    @test length(fixture.cases.two_coaxial_bundle_kron.terminal_order) == 4
    @test fixture.cases.two_coaxial_bundle_kron.phase_map == [1, 0, 1, 0]
    @test fixture.provenance.development_only_repairs !== nothing
end

@testitem "Gmsh FEM / optional frozen Python numerical comparisons" tags=[
    :extension,
    :integration
] begin
    using Gmsh
    using JSON3
    using LineCableModels
    using LinearAlgebra

    executable = get(ENV, "LINECABLEMODELS_GETDP", something(Sys.which("getdp"), ""))
    if isempty(executable)
        @test true
    else
        fixture_path = joinpath(
            pkgdir(LineCableModels),
            "test",
            "fixtures",
            "reference",
            "fem_python_quasi_tem.json"
        )
        fixture = JSON3.read(read(fixture_path, String))

        function complex_scan(values)
            frequency_count = length(values)
            row_count = length(values[1])
            column_count = length(values[1][1])
            scan = Array{ComplexF64, 3}(
                undef, row_count, column_count, frequency_count
            )
            for frequency in 1:frequency_count
                for row in 1:row_count, column in 1:column_count

                    pair = values[frequency][row][column]
                    scan[row, column, frequency] = complex(pair[1], pair[2])
                end
            end
            return scan
        end

        relative_error(actual, reference) = norm(actual - reference) /
                                            max(norm(reference), eps(Float64))

        copper = Material(
            kind = :conductor,
            rho = 1 / 5.8e7,
            alpha = 0.00393,
            T0 = 20.0
        )
        core_insulation = Material(kind = :insulator, rho = Inf, eps_r = 2.3)
        jacket = Material(kind = :insulator, rho = Inf, eps_r = 3.0)
        design = build(CableDesign,
            "python-reference-coax",
            Stack(
                Group(:core, Region(:core_metal, Disk(0.01), copper)),
                Region(
                    :core_insulation,
                    Annulus(0.01, 0.015),
                    core_insulation
                ),
                Group(
                    :sheath,
                    Region(:sheath_metal, Annulus(0.015, 0.017), copper)
                ),
                Region(:jacket, Annulus(0.017, 0.02), jacket)
            ))
        earth = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0)
        limits = fixture.provenance.comparison_metrics

        cases = (
            (
                name = :single_coaxial_kron,
                designs = design,
                positions = (0.0, -0.1),
                connections = Dict(:core => 1, :sheath => 0),
                maps = true
            ),
            (
                name = :two_coaxial_bundle_kron,
                designs = [design, design],
                positions = [(-0.03, -0.1), (0.03, -0.1)],
                connections = [
                    Dict(:core => 1, :sheath => 0),
                    Dict(:core => 1, :sheath => 0)
                ],
                maps = false
            )
        )

        for reference_case in cases
            expected = getproperty(fixture.cases, reference_case.name)
            system = build(
                LineCableSystem,
                reference_case.designs,
                reference_case.positions;
                connections = reference_case.connections,
                system_id = String(reference_case.name),
                line_length = 1.0
            )
            problem = LineParametersProblem(
                system;
                earth_props = earth,
                frequencies = Float64[expected.frequencies_hz...]
            )
            formulation = Formulation(
                :LineCableModelsFEM;
                options = (ideal_transposition = false,),
                fem_options = (
                    getdp_executable = executable,
                    getdp_verbosity = 0,
                    gmsh_verbosity = 0,
                    keep_run_directory = true,
                    plot_field_maps = reference_case.maps
                )
            )
            result = compute(problem, formulation; options = (trace = true,))
            run_directory = result.details.fem.run.run_directory
            try
                actual = result.details.fem.primitive
                expected_z = complex_scan(expected.Z_primitive_ohm_per_m)
                expected_p = complex_scan(
                    expected.P_primitive_inverse_siemens_per_m
                )
                expected_reduced_z = complex_scan(expected.Z_reduced_ohm_per_m)
                expected_reduced_y = complex_scan(expected.Y_reduced_siemens_per_m)

                errors = (
                    primitive_z = relative_error(actual.Z_primitive, expected_z),
                    primitive_p = relative_error(actual.P_primitive, expected_p),
                    reduced_z = relative_error(result.Z.values, expected_reduced_z),
                    reduced_y = relative_error(result.Y.values, expected_reduced_y)
                )
                @info "Frozen Python FEM comparison" case=reference_case.name errors
                @test errors.primitive_z <=
                      limits.primitive_relative_frobenius_limit
                @test errors.primitive_p <=
                      limits.primitive_relative_frobenius_limit
                @test errors.reduced_z <=
                      limits.reduced_relative_frobenius_limit
                @test errors.reduced_y <=
                      limits.reduced_relative_frobenius_limit
                @test result.f == Float64[expected.frequencies_hz...]
                @test result.details.fem.run.getdp_invocations == 1
                @test actual.phase_map == Int[expected.phase_map...]

                terminal_count = length(expected.terminal_order)
                expected_rows = length(expected.frequencies_hz) * terminal_count^2
                for filename in ("Z.tsv", "P.tsv")
                    lines = filter(
                        !isempty,
                        strip.(readlines(joinpath(run_directory, "raw", filename)))
                    )
                    @test length(lines) - 1 == expected_rows
                end
                completion = filter(
                    !isempty,
                    strip.(readlines(joinpath(
                        run_directory, "raw", "scan_complete.tsv"
                    )))
                )
                @test length(completion) == 2

                if reference_case.maps
                    @test length(result.details.fem.run.map_paths) ==
                          9 * terminal_count * length(expected.frequencies_hz)
                    @test all(isfile, result.details.fem.run.map_paths)
                    first_map = read(first(result.details.fem.run.map_paths), String)
                    @test occursin("f=10 Hz", first_map)
                    @test occursin("basis=cable_0001", first_map)
                else
                    @test isempty(result.details.fem.run.map_paths)
                end
            finally
                run_directory === nothing || rm(
                    run_directory; recursive = true, force = true
                )
            end
        end
    end
end
