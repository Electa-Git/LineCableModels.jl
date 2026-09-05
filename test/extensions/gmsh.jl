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
    mktempdir() do directory
        run = extension_module.FEMRun(
            directory, extension_module.created, "fixture", :none, ""
        )
        socket_name = extension_module._getdp_socket_name(run)
        if Sys.iswindows()
            @test socket_name == "127.0.0.1:0"
        else
            @test occursin("linecablemodels-gmsh-", basename(socket_name))
            @test !contains(socket_name, Base.Filesystem.path_separator)
        end
        mesh_source = joinpath(directory, "source.msh")
        mesh_snapshot = joinpath(directory, "snapshot.msh")
        write(mesh_source, "mesh snapshot fixture")
        extension_module._copy_or_link_mesh(mesh_source, mesh_snapshot)
        @test read(mesh_snapshot, String) == "mesh snapshot fixture"
        Sys.iswindows() || @test Base.Filesystem.samefile(
            mesh_source, mesh_snapshot
        )
    end

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
    resume_options = extension_module._fem_computation_options(
        (trace = true, resume_run_directory = :latest)
    )
    @test resume_options.trace === Val(true)
    @test resume_options.resume_run_directory === :latest
    @test_throws ArgumentError extension_module._fem_computation_options(
        (resume_run_directory = :invalid,)
    )
    @test extension_module._resume_value_matches(
        Dict("first" => 1, "second" => [2.0, 3.0]),
        Dict("second" => [2.0, 3.0], "first" => 1)
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

        job_path = joinpath(directory, "getdp-f0001-b0001-Z.tsv")
        open(job_path, "w") do io
            println(io, rows[1])
            println(io, rows[2])
        end
        @test extension_module._valid_job_raw(job_path, 2, 1, 50.0, 1)
        open(job_path, "w") do io
            println(io, rows[1])
            println(io, rows[1])
        end
        @test !extension_module._valid_job_raw(job_path, 2, 1, 50.0, 1)

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

@testitem "Gmsh FEM / nominal Float64 preflight" tags=[:extension] begin
    using Gmsh
    using LineCableModels
    using Measurements

    extension_module = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    @test extension_module !== nothing

    radius = measurement(0.005, 1.0e-6)
    copper = Material(kind = :conductor, rho = 1 / 5.8e7)
    design = build(
        CableDesign,
        "fem-preflight",
        Group(:core, Region(:core_metal, Disk(radius), copper))
    )
    system = build(
        LineCableSystem,
        design,
        (0.0, -0.1);
        connections = Dict(:core => 1),
        system_id = "fem-preflight",
        line_length = 1.0
    )
    problem = LineParametersProblem(
        system;
        temperature = measurement(35.0, 0.5),
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [measurement(50.0, 0.25)]
    )

    runtime_root = extension_module._runtime_root()
    runs = joinpath(runtime_root, "runs")
    before = isdir(runs) ? sort(readdir(runs)) : nothing
    @test !Bool(gmsh.is_initialized())

    normalized = extension_module._preflight_fem_problem(problem)

    @test eltype(problem) === Measurement{Float64}
    @test eltype(normalized) === Float64
    @test eltype(normalized.system) === Float64
    @test eltype(only(normalized.system.designs)) === Float64
    @test normalized.temperature === 35.0
    @test normalized.frequencies == [50.0]
    @test normalized.earth_props.layers[2].rho === 100.0
    @test only(normalized.system.designs).geometry.regions[1].primitive.r === 0.005
    @test Measurements.uncertainty(problem.temperature) == 0.5
    @test Measurements.uncertainty(only(problem.frequencies)) == 0.25
    @test !Bool(gmsh.is_initialized())
    @test (isdir(runs) ? sort(readdir(runs)) : nothing) == before

    float32_design = build(
        CableDesign,
        "fem-preflight-f32",
        Group(
            :core,
            Region(
                :core_metal,
                Disk(0.005f0),
                Material(kind = :conductor, rho = Float32(1 / 5.8e7))
            )
        )
    )
    float32_system = build(
        LineCableSystem,
        float32_design,
        (0.0f0, -0.1f0);
        connections = Dict(:core => 1),
        system_id = "fem-preflight-f32",
        line_length = 1.0f0
    )
    float32_problem = LineParametersProblem(
        float32_system;
        temperature = 20.0f0,
        earth_props = LineCableModels.EarthProps.EarthModel(
            100.0f0, 10.0f0, 1.0f0
        ),
        frequencies = Float32[50.0]
    )
    promoted = extension_module._preflight_fem_problem(float32_problem)
    @test eltype(promoted) === Float64
    @test promoted.frequencies == [50.0]
    @test only(promoted.system.designs).geometry.regions[1].primitive.r ===
          Float64(0.005f0)

    unsupported_problem = LineParametersProblem(
        system;
        temperature = measurement(35.0, 0.5),
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [measurement(50.0, 0.25)],
        Γ = [complex(measurement(0.0, 0.0), measurement(1.0e-12, 1.0e-14))]
    )
    formulation = Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,),
        fem_options = (gmsh_verbosity = 0,)
    )
    exception = try
        compute(unsupported_problem, formulation)
        nothing
    catch caught
        caught
    end
    @test exception isa LineCableModelsFEMError
    @test exception.category === :unsupported
    @test exception.field === :Γ
    @test exception.run_directory === nothing
    @test !Bool(gmsh.is_initialized())
    @test (isdir(runs) ? sort(readdir(runs)) : nothing) == before
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
    @test model.tags.outer_air_boundary == 2_004
    @test model.tags.outer_earth_boundary == 2_005
    @test model.shell_outer_radius == 1.25model.domain_radius
    @test getproperty.(model.region_plans, :mesh_size) ≈ [
        0.001, 0.0025, 0.0005, 0.001
    ]
    @test model.fine_mesh_size ≈ 0.0005
    @test model.cable_outer_mesh_sizes ≈ [0.001]
    @test model.mesh_growth_factor == 1.2
    @test length(model.mesh_plans) == 1
    @test only(model.mesh_plans).domain_mesh_size == model.domain_radius / 20
    @test only(model.mesh_plans).infinite_mesh_size == model.domain_radius / 10
    @test only(model.mesh_plans).cable_interface_mesh_sizes ≈ [
        min(
        model.domain_radius / 20,
        0.001 + 0.2 * (0.1 - 0.013)
    )
    ]

    multifrequency_problem = LineParametersProblem(
        system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0, 1000.0]
    )
    multifrequency_model = extension_module._resolved_fem_model(
        multifrequency_problem, formulation
    )
    @test getproperty.(multifrequency_model.mesh_plans, :frequency) == [
        50.0, 1000.0
    ]
    @test multifrequency_model.mesh_plans[1].domain_radius >
          multifrequency_model.mesh_plans[2].domain_radius
    @test multifrequency_model.domain_radius ==
          multifrequency_model.mesh_plans[2].domain_radius

    matrix = Material(kind = :insulator, rho = Inf, eps_r = 1.0)
    sector_core = Stack(
        Group(:core, Region(:core_centre, Disk(0.001), copper)),
        Group(
            :core,
            Region(
                :core_sectors,
                Sector(span = pi / 3, r_base = 0.001, r_back = 0.002),
                copper
            );
            pattern = Ring(6; r = 0.0)
        )
    )
    filled_design = build(
        CableDesign,
        "fem-filled-sector",
        Enclosure(
            :matrix,
            sector_core;
            primitive = Disk(0.002),
            fill = matrix
        ),
        Region(:outer_insulation, Annulus(0.002, 0.003), dielectric)
    )
    filled_system = build(
        LineCableSystem,
        filled_design,
        (0.0, -0.1);
        connections = Dict(:core => 1),
        system_id = "fem-filled-sector",
        line_length = 1.0
    )
    filled_problem = LineParametersProblem(
        filled_system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0]
    )
    filled_model = extension_module._resolved_fem_model(
        filled_problem, formulation
    )
    @test length(filled_model.region_plans) == 9
    @test count(plan -> plan.field === :matrix_fill,
        filled_model.material_plans) == 1

    armor_wire_radius = 0.002
    armor_inner_radius = 0.010
    armor_outer_radius = 0.014
    tangent_fill_design = build(
        CableDesign,
        "fem-tangent-circular-fill",
        Stack(
            Region(:inner_bedding, Disk(armor_inner_radius), dielectric),
            Enclosure(
                :armor_matrix,
                Group(
                    :armor,
                    Region(:armor_wires, Disk(armor_wire_radius), copper);
                    pattern = Ring(8; r = 0.012)
                );
                primitive = Annulus(armor_inner_radius, armor_outer_radius),
                fill = matrix
            ),
            Region(
                :outer_jacket,
                Annulus(armor_outer_radius, 0.016),
                dielectric
            )
        )
    )
    tangent_fill_system = build(
        LineCableSystem,
        tangent_fill_design,
        (0.0, -0.1);
        connections = Dict(:armor => 1),
        system_id = "fem-tangent-circular-fill",
        line_length = 1.0
    )
    tangent_fill_problem = LineParametersProblem(
        tangent_fill_system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0]
    )
    tangent_fill_model = extension_module._resolved_fem_model(
        tangent_fill_problem, formulation
    )
    tangent_fill_index = findfirst(
        plan -> plan.field === :armor_matrix_fill,
        tangent_fill_model.material_plans
    )
    @test tangent_fill_index !== nothing

    tape_inner_radius = armor_outer_radius
    tape_outer_radius = tape_inner_radius + 0.0002
    mixed_fill_design = build(
        CableDesign,
        "fem-tangent-fill-with-tape",
        Stack(
            Region(:inner_bedding, Disk(armor_inner_radius), dielectric),
            Enclosure(
                :screen_matrix,
                Stack(
                    Group(
                        :screen,
                        Region(:screen_wires, Disk(armor_wire_radius), copper);
                        pattern = Ring(8; r = 0.012)
                    ),
                    Group(
                        :screen,
                        Region(
                            :screen_tape,
                            Rectangle(
                                0.4 * (tape_inner_radius + tape_outer_radius) / 2,
                                tape_outer_radius - tape_inner_radius
                            ),
                            copper
                        );
                        pattern = Ring(
                            1;
                            r = (tape_inner_radius + tape_outer_radius) / 2
                        )
                    )
                );
                primitive = Annulus(armor_inner_radius, tape_outer_radius),
                fill = matrix
            ),
            Region(
                :outer_jacket,
                Annulus(tape_outer_radius, 0.016),
                dielectric
            )
        )
    )
    mixed_fill_system = build(
        LineCableSystem,
        mixed_fill_design,
        (0.0, -0.1);
        connections = Dict(:screen => 1),
        system_id = "fem-tangent-fill-with-tape",
        line_length = 1.0
    )
    mixed_fill_problem = LineParametersProblem(
        mixed_fill_system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0]
    )
    mixed_fill_model = extension_module._resolved_fem_model(
        mixed_fill_problem, formulation
    )
    mixed_fill_index = findfirst(
        plan -> plan.field === :screen_matrix_fill,
        mixed_fill_model.material_plans
    )
    @test mixed_fill_index !== nothing

    incomplete_design = build(
        CableDesign,
        "fem-incomplete-partition",
        Group(:core, Region(:core_metal, Disk(0.001), copper)),
        Region(:outer_insulation, Annulus(0.002, 0.003), dielectric)
    )
    incomplete_system = build(
        LineCableSystem,
        incomplete_design,
        (0.0, -0.1);
        connections = Dict(:core => 1),
        system_id = "fem-incomplete-partition",
        line_length = 1.0
    )
    incomplete_problem = LineParametersProblem(
        incomplete_system;
        earth_props = LineCableModels.EarthProps.EarthModel(100.0, 10.0, 1.0),
        frequencies = [50.0]
    )
    partition_error = try
        extension_module._resolved_fem_model(incomplete_problem, formulation)
        nothing
    catch exception
        exception
    end
    @test partition_error isa LineCableModelsFEMError
    @test partition_error.category === :adaptation
    @test partition_error.field === :material_partition

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
    gmsh.option.set_string("Solver.SocketName", "caller-owned-socket")
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
                    LineCableModels.DataModel.SectorShape(
                        Sector(
                            span = π / 3,
                            r_base = 0.005,
                            r_back = 0.02,
                            fillet = 0.001
                        ),
                        Pose2(0.05, 0.0, -π / 11)
                    )
                )
                gmsh.model.geo.synchronize()
                @test length(ellipse_loop.curves) == 4
                @test length(sector_loop.curves) == 10
                @test all(
                    gmsh.model.get_type(1, curve) == "Ellipse"
                for curve in ellipse_loop.curves
                )
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Circle",
                    sector_loop.curves
                ) == 8
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Line",
                    sector_loop.curves
                ) == 2

                gmsh.model.add("fem-zero-inner-strip")
                zero_inner_registry = extension_module.FEMLoopRegistry(1.0e-3)
                zero_inner_shape = LineCableModels.DataModel.BentStrip(
                    0.0, 0.01, π / 3
                )
                extension_module._register_shape_breaks!(
                    zero_inner_registry, zero_inner_shape
                )
                @test all(
                    key -> !iszero(key[3]),
                    keys(zero_inner_registry.circle_breaks)
                )
                zero_inner_strip = extension_module._boundary_loop!(
                    zero_inner_registry,
                    zero_inner_shape
                )
                gmsh.model.geo.synchronize()
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Circle",
                    zero_inner_strip.curves
                ) >= 1
                @test count(
                    curve -> gmsh.model.get_type(1, curve) == "Line",
                    zero_inner_strip.curves
                ) == 2

                gmsh.model.add("fem-shared-circle-arcs")
                shared_registry = extension_module.FEMLoopRegistry(1.0e-3)
                contact = 0.5343318251524625
                first_contact = extension_module._point!(
                    shared_registry, (contact, -1.0)
                )
                second_contact = extension_module._point!(
                    shared_registry, (nextfloat(contact), -1.0)
                )
                @test first_contact == second_contact
                full_circle = Disk(0.02)
                partial_sector = LineCableModels.DataModel.SectorShape(
                    Sector(span = 0.8, r_base = 0.01, r_back = 0.02)
                )
                extension_module._register_shape_breaks!(
                    shared_registry, full_circle
                )
                extension_module._register_shape_breaks!(
                    shared_registry, partial_sector
                )
                extension_module._register_circle_contacts!(shared_registry)
                full_loop = extension_module._boundary_loop!(
                    shared_registry, full_circle
                )
                partial_loop = extension_module._boundary_loop!(
                    shared_registry, partial_sector; mesh_size = 2.0e-4
                )
                gmsh.model.geo.synchronize()
                shared_curves = intersect(full_loop.curves, partial_loop.curves)
                @test length(shared_curves) == 2
                @test all(
                    point -> shared_registry.point_sizes[point] == 2.0e-4,
                    Iterators.flatten(
                        shared_registry.curve_points[curve]
                    for curve in shared_curves
                    )
                )

                gmsh.model.add("fem-tangent-circles")
                tangent_registry = extension_module.FEMLoopRegistry(1.0e-3)
                first_circle = Disk(0.01, Pose2(-0.01, 0.0, 0.0))
                second_circle = Disk(0.01, Pose2(0.01, 0.0, 0.0))
                extension_module._register_shape_breaks!(
                    tangent_registry, first_circle
                )
                extension_module._register_shape_breaks!(
                    tangent_registry, second_circle
                )
                extension_module._register_circle_contacts!(tangent_registry)
                first_loop = extension_module._boundary_loop!(
                    tangent_registry, first_circle
                )
                second_loop = extension_module._boundary_loop!(
                    tangent_registry, second_circle
                )
                gmsh.model.geo.synchronize()
                first_points = reduce(
                    union,
                    (Set(gmsh.model.get_adjacencies(1, curve)[2])
                    for curve in first_loop.curves);
                    init = Set{Int32}()
                )
                second_points = reduce(
                    union,
                    (Set(gmsh.model.get_adjacencies(1, curve)[2])
                    for curve in second_loop.curves);
                    init = Set{Int32}()
                )
                @test length(intersect(first_points, second_points)) == 1

                filled_geometry = extension_module._build_geometry!(
                    filled_model, "fem-filled-sector-geometry"
                )
                fill_index = findfirst(
                    plan -> plan.field === :matrix_fill,
                    filled_model.material_plans
                )
                fill_surface = only(filled_geometry.material_surfaces[fill_index])
                fill_curves = Set(extension_module._entity_boundary([fill_surface]))
                terminal_curves = Set(only(filled_geometry.terminal_curves))
                @test !isempty(intersect(terminal_curves, fill_curves))
                @test all(terminal_curves) do curve
                    adjacent, _ = gmsh.model.get_adjacencies(1, curve)
                    length(unique(adjacent)) >= 2
                end

                tangent_fill_geometry = extension_module._build_geometry!(
                    tangent_fill_model, "fem-tangent-fill-geometry"
                )
                tangent_fill_surfaces = tangent_fill_geometry.material_surfaces[
                    tangent_fill_index
                ]
                @test length(tangent_fill_surfaces) == 8
                tangent_fill_curves = extension_module._entity_boundary(
                    tangent_fill_surfaces
                )
                @test all(tangent_fill_curves) do curve
                    adjacent, _ = gmsh.model.get_adjacencies(1, curve)
                    length(unique(adjacent)) >= 2
                end
                minimum_curve_extent = minimum(tangent_fill_curves) do curve
                    bounds = gmsh.model.get_bounding_box(1, curve)
                    hypot(bounds[4] - bounds[1], bounds[5] - bounds[2])
                end
                @test minimum_curve_extent > 1.0e-8
                gmsh.model.set_current(tangent_fill_geometry.model_name)
                extension_module._configure_mesh!(
                    tangent_fill_model,
                    tangent_fill_geometry,
                    only(tangent_fill_model.mesh_plans)
                )
                gmsh.model.mesh.generate(2)
                @test !isempty(first(gmsh.model.mesh.get_nodes()))

                mixed_fill_geometry = extension_module._build_geometry!(
                    mixed_fill_model, "fem-tangent-fill-with-tape-geometry"
                )
                mixed_fill_surfaces = mixed_fill_geometry.material_surfaces[
                    mixed_fill_index
                ]
                @test length(mixed_fill_surfaces) == 9
                gmsh.model.set_current(mixed_fill_geometry.model_name)
                extension_module._configure_mesh!(
                    mixed_fill_model,
                    mixed_fill_geometry,
                    only(mixed_fill_model.mesh_plans)
                )
                gmsh.model.mesh.generate(2)
                @test !isempty(first(gmsh.model.mesh.get_nodes()))

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
                @test !isempty(first_geometry.outer_air_curves)
                @test !isempty(first_geometry.outer_earth_curves)
                @test isempty(intersect(
                    first_geometry.outer_air_curves,
                    first_geometry.outer_earth_curves
                ))
                @test sort!(unique([first_geometry.outer_air_curves;
                                    first_geometry.outer_earth_curves])) ==
                      first_geometry.outer_curves
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
        @test gmsh.option.get_string("Solver.SocketName") ==
              "caller-owned-socket"
        @test gmsh.onelab.get_string(caller_parameter) == ["preserve-me"]
    finally
        Gmsh.finalize()
    end
end

@testitem "Gmsh FEM / bounded formations and complete material ownership" tags=[
    :extension
] begin
    using Gmsh
    using LineCableModels

    const DM = LineCableModels.DataModel
    extension_module = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    formulation = Formulation(
        :LineCableModelsFEM; options = (ideal_transposition = false,)
    )
    copper = Material(
        kind = :conductor, rho = 1.72e-8, eps_r = 1.0, mu_r = 1.0
    )
    dielectric = Material(
        kind = :insulator, rho = Inf, eps_r = 2.3, mu_r = 1.0
    )

    function problem(design, connections, id)
        system = build(
            LineCableSystem,
            design,
            Pose2(0.0, -0.1);
            connections,
            system_id = id,
            line_length = 1.0
        )
        return LineParametersProblem(
            system;
            earth_props = homogeneous(rho = 100.0, eps_r = 10.0),
            frequencies = [50.0]
        )
    end

    strand_radius = 0.5e-3
    strand_count = 19
    circular_boundary = Disk(sqrt(strand_count) * strand_radius)
    circular = stranded(
        copper;
        shape = Disk(strand_radius),
        compact = true,
        boundary = circular_boundary
    )
    circular_design = build(
        CableDesign,
        "fem-compacted-circular-core",
        terminal(:core, circular),
        Region(:insulation, Shell(1e-3), dielectric)
    )
    circular_model = extension_module._resolved_fem_model(
        problem(circular_design, Dict(:core => 1), circular_design.cable_id),
        formulation
    )
    @test length(circular_model.region_plans) == 2
    @test circular_model.region_plans[1].shape isa DM.Disk
    @test circular_model.region_plans[1].shape.r == circular_boundary.r
    @test circular_model.region_plans[1].shape.at == Pose2(0.0, -0.1)
    @test circular_model.region_plans[1].terminal_index == 1

    natural = terminal(
        :core,
        stranded(
            copper;
            shape = Disk(strand_radius),
            boundary = Disk(3strand_radius),
            fill = dielectric
        )
    )
    natural_design = build(CableDesign, "fem-natural-circular-core", natural)
    natural_model = extension_module._resolved_fem_model(
        problem(natural_design, Dict(:core => 1), natural_design.cable_id),
        formulation
    )
    @test length(natural_model.region_plans) == 8
    @test count(
        plan -> plan.shape isa DM.Disk && plan.terminal_index == 1,
        natural_model.region_plans
    ) == 7
    @test count(
        plan -> plan.shape isa DM.DifferenceShape && plan.terminal_index == 0,
        natural_model.region_plans
    ) == 1

    partial_boundary = Disk(sqrt(7 / 0.9) * strand_radius)
    partial = terminal(
        :core,
        stranded(
            copper;
            shape = Disk(strand_radius),
            compact = true,
            boundary = partial_boundary,
            fill = dielectric
        )
    )
    filled_design = build(
        CableDesign,
        "fem-filled-compaction",
        partial
    )
    filled_model = extension_module._resolved_fem_model(
        problem(filled_design, Dict(:core => 1), filled_design.cable_id),
        formulation
    )
    @test length(filled_model.region_plans) == 8
    @test count(
        plan -> plan.shape isa DM.Polygon,
        filled_model.region_plans
    ) == 7
    @test count(
        plan -> plan.shape isa DM.DifferenceShape,
        filled_model.region_plans
    ) == 1

    sector = Sector(
        span = 2pi / 3,
        r_base = 0.6e-3,
        r_back = 4.0e-3,
        fillet = 0.2e-3
    )
    sector_shape = DM.resolve(DM.EmptyBoundary(), sector)
    sector_strand_radius = sqrt(0.9DM.area(sector_shape) / (7pi))
    sector_part = stranded(
        copper;
        shape = Disk(sector_strand_radius),
        boundary = sector
    )
    sector_names = (:a, :b, :c)
    sectors = assembly((
        at(
            terminal(name, sector_part),
            0.15e-3cos(2pi * (index - 1) / 3),
            0.15e-3sin(2pi * (index - 1) / 3);
            φ = 2pi * (index - 1) / 3
        ) for (index, name) in enumerate(sector_names)
    )...)
    sector_design = build(
        CableDesign,
        "fem-sector-formations",
        Enclosure(
            :sector_matrix,
            sectors;
            primitive = Disk(4.5e-3),
            fill = dielectric
        )
    )
    sector_model = extension_module._resolved_fem_model(
        problem(
            sector_design,
            Dict(:a => 1, :b => 2, :c => 3),
            sector_design.cable_id
        ),
        formulation
    )
    @test length(sector_model.region_plans) == 25
    @test count(
        plan -> plan.shape isa DM.Polygon,
        sector_model.region_plans
    ) == 21
    @test [count(==(terminal), getproperty.(sector_model.region_plans, :terminal_index))
           for terminal in 1:3] == [7, 7, 7]
    matrix_plan = only(filter(
        plan -> plan.terminal_index == 0 &&
                plan.shape isa DM.DifferenceShape &&
                plan.shape.outer isa DM.Disk,
        sector_model.region_plans
    ))
    @test matrix_plan.shape isa DM.DifferenceShape
    @test count(
        hole -> hole isa DM.SectorShape,
        matrix_plan.shape.holes
    ) == 3
    @test all(hole -> !(hole isa DM.Polygon), matrix_plan.shape.holes)
    @test count(
        plan -> plan.terminal_index == 0 &&
                plan.shape isa DM.DifferenceShape &&
                plan.shape.outer isa DM.SectorShape,
        sector_model.region_plans
    ) == 3

    four_sector = Sector(
        span = pi / 2,
        r_base = 0.6e-3,
        r_back = 4.0e-3,
        fillet = 0.2e-3
    )
    four_sector_shape = DM.resolve(DM.EmptyBoundary(), four_sector)
    four_strand_radius = sqrt(0.9DM.area(four_sector_shape) / (7pi))
    four_part = stranded(
        copper;
        shape = Disk(four_strand_radius),
        boundary = four_sector
    )
    four_names = (:d, :e, :f, :g)
    four_assembly = assembly((
        at(
            terminal(name, four_part),
            0.15e-3cos(2pi * (index - 1) / 4),
            0.15e-3sin(2pi * (index - 1) / 4);
            φ = 2pi * (index - 1) / 4
        ) for (index, name) in enumerate(four_names)
    )...)
    four_design = build(
        CableDesign,
        "fem-four-sector-formations",
        Enclosure(
            :four_sector_matrix,
            four_assembly;
            primitive = Disk(4.5e-3),
            fill = dielectric
        )
    )
    four_model = extension_module._resolved_fem_model(
        problem(
            four_design,
            Dict(:d => 1, :e => 2, :f => 3, :g => 4),
            four_design.cable_id
        ),
        formulation
    )
    @test length(four_model.region_plans) == 33
    @test count(plan -> plan.shape isa DM.Polygon, four_model.region_plans) == 28
    @test [count(==(terminal), getproperty.(four_model.region_plans, :terminal_index))
           for terminal in 1:4] == [7, 7, 7, 7]
    @test count(
        plan -> plan.terminal_index == 0 &&
                plan.shape isa DM.DifferenceShape &&
                plan.shape.outer isa DM.SectorShape,
        four_model.region_plans
    ) == 4

    milliken_core = milliken(
        copper;
        shape = Disk(0.33e-3),
        segment = Sector(
            span = pi / 3,
            r_base = 0.85e-3,
            r_back = 3.0e-3,
            fillet = 0.1e-3
        )
    )
    milliken_design = build(
        CableDesign,
        "fem-milliken-core",
        terminal(:core, milliken_core)
    )
    milliken_model = extension_module._resolved_fem_model(
        problem(
            milliken_design,
            Dict(:core => 1),
            milliken_design.cable_id
        ),
        formulation
    )
    @test count(
        plan -> plan.terminal_index == 1,
        milliken_model.region_plans
    ) > 7
    milliken_fill = only(filter(
        plan -> plan.terminal_index == 0,
        milliken_model.region_plans
    ))
    @test milliken_fill.shape isa DM.DifferenceShape
    @test length(milliken_fill.shape.holes) ==
          count(plan -> plan.terminal_index == 1, milliken_model.region_plans)

    Gmsh.initialize(String[]; finalize_atexit = false)
    try
        gmsh.option.set_number("General.Terminal", 1)
        gmsh.option.set_number("General.Verbosity", 0)
        for (name, model) in (
            ("fem-compacted-circular-core", circular_model),
            ("fem-natural-circular-core", natural_model),
            ("fem-filled-compaction", filled_model),
            ("fem-sector-formations", sector_model),
            ("fem-four-sector-formations", four_model),
            ("fem-milliken-core", milliken_model)
        )
            geometry = extension_module._build_geometry!(model, name)
            gmsh.model.set_current(geometry.model_name)
            @test all(!isempty, geometry.material_surfaces)
            @test all(!isempty, geometry.terminal_surfaces)
            name in (
                "fem-compacted-circular-core",
                "fem-sector-formations",
                "fem-milliken-core"
            ) || continue
            extension_module._configure_mesh!(
                model, geometry, only(model.mesh_plans)
            )
            gmsh.model.mesh.generate(2)
            @test !isempty(first(gmsh.model.mesh.get_nodes()))
        end
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
                gmsh_verbosity = 0,
                keep_run_directory = true
            )
        )
        result = compute(problem, formulation; options = (trace = true,))
        run_directory = result.details.fem.run.run_directory
        @test size(result.Z) == (1, 1, 2)
        @test size(result.Y) == (1, 1, 2)
        @test all(isfinite, result.Z)
        @test all(isfinite, result.Y)
        @test result.f == [50.0, 1000.0]
        expected_capacitance = 2π * 8.8541878128e-12 * dielectric.eps_r / log(0.01 / 0.005)
        fem_capacitance = [imag(result.Y[1, 1, index]) / (2π * frequency)
                           for (index, frequency) in pairs(result.f)]
        @test all(
            isapprox(value, expected_capacitance; rtol = 0.02)
        for value in fem_capacitance
        )
        @test result.details.fem.run.state === Base.get_extension(
            LineCableModels, :LineCableModelsGmshExt
        ).completed
        @test result.details.fem.run.getdp_invocations ==
              length(problem.frequencies) * length(result.details.fem.terminal_ids)
        @test run_directory !== nothing
        expected_rows = 1 +
                        length(problem.frequencies) *
                        length(result.details.fem.terminal_ids)^2
        for filename in ("Z.tsv", "P.tsv")
            content = read(joinpath(run_directory, "raw", filename), String)
            @test !occursin("\n\n", content)
            @test length(readlines(IOBuffer(content))) == expected_rows
        end

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
            @test ui_result.details.fem.run.getdp_invocations ==
                  length(problem.frequencies) *
                  length(ui_result.details.fem.terminal_ids)
            @test ui_result.details.fem.run.run_directory === nothing
        end
        rm(run_directory; recursive = true, force = true)
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

    getdp_model = read(
        joinpath(
            pkgdir(LineCableModels),
            "ext",
            "LineCableModelsGmshExt",
            "getdp",
            "model.pro"
        ),
        String)
    @test occursin("AirInfJacobian = Region[{AIR_INF}]", getdp_model)
    @test occursin("EarthInfJacobian = Region[{EARTH_INF}]", getdp_model)
    @test occursin(
        "DomainInf = Region[{AirInfJacobian, EarthInfJacobian}]",
        getdp_model
    )
    @test !occursin("DomainInf = Region[{DOMAIN_INF}]", getdp_model)
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
                terminal_count = length(expected.terminal_order)
                @test result.details.fem.run.getdp_invocations ==
                      length(expected.frequencies_hz) * terminal_count
                @test actual.phase_map == Int[expected.phase_map...]

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
