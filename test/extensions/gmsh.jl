@testitem "Gmsh FEM / public API and strict parsing" tags=[:extension] begin
    import LineCableModels
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    using LinearAlgebra

    extension_module = Base.get_extension(
        LineCableModels, :LineCableModelsGmshExt
    )
    @test extension_module !== nothing
    @test Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM <:
          LineCableModels.AbstractFormulation
    @test supertype(Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM) === LineCableModels.AbstractFormulation
    mktempdir() do directory
        run = extension_module.FEMRun(
            directory, extension_module.created, "fixture", :none, ""
        )
        mesh_source = joinpath(directory, "source.msh")
        mesh_snapshot = joinpath(directory, "snapshot.msh")
        write(mesh_source, "mesh snapshot fixture")
        extension_module._copy_or_link_mesh(mesh_source, mesh_snapshot)
        @test read(mesh_snapshot, String) == "mesh snapshot fixture"
        Sys.iswindows() || @test Base.Filesystem.samefile(
            mesh_source, mesh_snapshot
        )
    end

    options = computation_options(LineCableModelsFEM, ComputationOptions((;mesh_policy = :remesh,
        keep_run_directory = true,
        gmsh_verbosity = 0,
        getdp_verbosity = 5)))
    @test options.data.mesh_policy === :remesh
    @test options.data.keep_run_directory
    @test options.data.getdp_verbosity == 5
    @test options.data.frequency_workers == clamp(Sys.CPU_THREADS ÷ 4, 1, 8)
    @test options.data.mumps_ordering == 0
    @test options.data.solver_threads == 1
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((;frequency_workers = 0)))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((;solver_threads = -1)))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((;mesh_policy = :invalid)))
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((;getdp_verbosity = 6)))
    resume_options = computation_options(LineCableModelsFEM, ComputationOptions((trace = true, resume_run_directory = :latest)))
    @test resume_options.data.trace === Val(true)
    @test resume_options.data.resume_run_directory === :latest
    @test_throws ArgumentError computation_options(LineCableModelsFEM, ComputationOptions((resume_run_directory = :invalid,)))
    @test extension_module._resume_value_matches(
        Dict("first" => 1, "second" => [2.0, 3.0]),
        Dict("second" => [2.0, 3.0], "first" => 1)
    )

    formulation = LineCableModels.Formulation(
        :LineCableModelsFEM;
        options = (ideal_transposition = false,))
    @test formulation isa Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    @test !formulation.options.data.ideal_transposition
    @test !hasproperty(formulation, :execution)

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
        @test_throws Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError extension_module._parse_raw_matrix(
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
        @test missing_completion isa Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
        @test missing_completion.run_directory == directory

        open(path, "w") do io
            println(io, header)
            println.(Ref(io), ["1\t50\t1\t1\tNaN\t0", rows[2], rows[3], rows[4]])
        end
        @test_throws Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError extension_module._parse_raw_matrix(
            Float64, path, [50.0], 2, run
        )
    end

    potential = Array{ComplexF64, 3}(undef, 2, 2, 2)
    potential[:, :, 1] = [2.0 0.25; 0.25 1.5]
    potential[:, :, 2] = [3.0 + 0.1im 0.5; 0.5 2.0 + 0.2im]
    inversion = extension_module.potential_to_admittance(
        potential; diagnostics = true
    )
    for frequency in axes(potential, 3)
        @test potential[:, :, frequency] * inversion.Y[:, :, frequency] ≈ I
        @test inversion.residuals[frequency] ≤ 100eps(Float64)
    end
    # An overflowing condition estimate is not a failed solve. The physical
    # inverse is finite and remains the result of the original LU operation.
    ill_conditioned = reshape(ComplexF64[1e-200 0; 0 1e200], 2, 2, 1)
    warned = @test_logs (:warn, r"condition estimate is not finite") begin
        extension_module.potential_to_admittance(ill_conditioned; diagnostics = true)
    end
    @test all(isfinite, warned.Y)
    @test warned.Y[:, :, 1] == inv(ill_conditioned[:, :, 1])
    @test isinf(only(warned.condition_numbers))
    @test only(warned.residuals) ≤ 100eps(Float64)
    @test_throws ArgumentError extension_module.potential_to_admittance(
        fill(ComplexF64(NaN), 2, 2, 1))
    @test_throws SingularException extension_module.potential_to_admittance(
        zeros(ComplexF64, 2, 2, 1))
end

@testitem "Gmsh FEM / nominal Float64 preflight" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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
        earth_props = LineCableModels.Earth.EarthModel(100.0, 10.0, 1.0),
        frequencies = [measurement(50.0, 0.25)]
    )

    runtime_root = extension_module._runtime_root()
    runs = extension_module._system_run_root(runtime_root,problem.system.system_id)
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
        earth_props = LineCableModels.Earth.EarthModel(
            100.0f0, 10.0f0, 1.0f0
        ),
        frequencies = Float32[50.0]
    )
    promoted = extension_module._preflight_fem_problem(float32_problem)
    @test eltype(promoted) === Float64
    @test promoted.frequencies == [50.0]
    @test only(promoted.system.designs).geometry.regions[1].primitive.r ===
          Float64(0.005f0)

    @test_throws MethodError LineParametersProblem(
        system;
        temperature = measurement(35.0, 0.5),
        earth_props = LineCableModels.Earth.EarthModel(100.0, 10.0, 1.0),
        frequencies = [measurement(50.0, 0.25)],
        Γ = [complex(measurement(0.0, 0.0), measurement(1.0e-12, 1.0e-14))]
    )
    @test !Bool(gmsh.is_initialized())
    @test (isdir(runs) ? sort(readdir(runs)) : nothing) == before
end

@testitem "Gmsh FEM / bounded formations and complete material ownership" tags=[
    :extension
] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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

    annular_design = build(CableDesign, "fem-rectangular-last-strip",
        terminal(:core,
            stranded(copper; shape = Rectangle(0.3e-3, 0.1e-3),
                center = Disk(0.2e-3), boundary = Disk(0.6e-3), lay = LayRatio(12)),
            insulation(dielectric; t = 0.2e-3)))
    annular_model = extension_module._resolved_fem_model(
        problem(annular_design, Dict(:core => 1), annular_design.cable_id),
        formulation)
    @test any(
        region -> region.source.primitive isa Rectangle &&
                  region.primitive isa Annulus,
        annular_design.geometry.regions)
    @test all(region -> region.source.tag !== :stranded_fill,
        annular_design.geometry.regions)
    @test length(annular_model.region_plans) == 2
    @test first(annular_model.region_plans).shape isa Disk
    @test first(annular_model.region_plans).shape.r ==
          last(annular_design.geometry.regions).primitive.ri

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
            ("fem-milliken-core", milliken_model),
            ("fem-rectangular-last-strip", annular_model)
        )
            @testset "$name" begin
                geometry = extension_module._build_physical_geometry!(model,name)
                @test all(!isempty,geometry.material_surfaces)
                @test all(!isempty,geometry.terminal_surfaces)
                for surfaces in geometry.material_surfaces,surface in surfaces
                    @test !isempty(gmsh.model.get_boundary([(2,surface)],false,false,false))
                end

            end
        end
        @testset "contact edges share subdivisions in either direction" begin
            gmsh.model.add("contact-edge-subdivision")
            registry = extension_module.FEMLoopRegistry(1e-12)
            for point in ((0.0, 0.0), (0.5e-3, 0.0), (1e-3, 0.0))
                extension_module._point!(registry, point; mesh_size = nothing)
            end
            @test isempty(registry.point_sizes)
            forward = extension_module._line_path!(
                registry, (0.0, 0.0), (1e-3, 0.0); mesh_size = 1e-4
            )
            backward = extension_module._line_path!(
                registry, (1e-3, 0.0), (0.0, 0.0); mesh_size = 1e-4
            )
            @test length(forward) == 2
            @test backward == -reverse(forward)
            @test all(==(1e-4), values(registry.point_sizes))
        end
    finally
        Gmsh.finalize()
    end
end


@testitem "Gmsh FEM / native multi-frequency scan and completed reuse" tags=[
    :extension,
    :integration,
    :fem_numerical
] setup=[TemporaryFEMRuntime] begin
    using LineCableModels
    using Gmsh
    using SHA
    using Serialization, JSON3, Logging

    cd(fem_test_runtime_directory)
    try
        copper = Material(kind = :conductor, rho = 1 / 5.8e7)
        dielectric = Material(kind = :insulator, rho = 1.0e8, eps_r = 2.3,
            tan_delta = 0.025)
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
            earth_props = LineCableModels.Earth.EarthModel(100.0, 10.0, 1.0),
            frequencies = [50.0, 1000.0]
        )
        formulation = Formulation(
            :LineCableModelsFEM;
            options = (ideal_transposition = false,))
        formulation_controls = (
                getdp_verbosity = 0,
                gmsh_verbosity = 0,
                keep_run_directory = true, timing = true
            )
        formulation_space = Formulation(:LineCableModelsFEM;
            earth_properties = Grid((formula(:default), nothing)),
            insulation_admittance = Grid((:default, :lossy)),
            options = formulation.options)
        selected = collect(formulation_space)
        completions = Tuple[]
        on_result = (resolved, index, result) -> push!(completions, (index, result))
        batch = compute(ParametricProblem(problem, ComputationOptions((;formulation_controls..., trace = true, on_result = on_result))),
            Combinatorial(formulation_space; options = (retain_details = true,)))
        @test length(batch) == 4
        @test first.(completions) == collect(eachindex(selected))
        @test all(index -> completions[index][2].Z.values == batch[index].Z.values,
            eachindex(selected))
        run_directories = [value.details.data.fem.run.run_directory for value in batch]
        @test length(unique(run_directories)) == 2
        for index in eachindex(selected)
            files = details(batch[index]).data.files
            @test any(file -> file.path == "run.json", files)
            @test all(file -> bytes2hex(open(sha256, file.source)) == file.sha256, files)
            soil = batch[index].details.data.formulations.methods.earth_properties
            @test soil === nothing ? selected[index].methods.earth_properties === nothing :
                  soil.identifier === formula_id(selected[index].methods.earth_properties)
            @test keys(batch[index].details.data.formulations.methods) == keys(selected[index].methods)
            @test details(batch).data.points[index] == details(batch[index])
        end
        default_indices = findall(
            value -> formula_id(value.methods.insulation_admittance) === :lossless, selected)
        lossy_indices = findall(
            value -> formula_id(value.methods.insulation_admittance) === :lossy, selected)
        for indices in (default_indices, lossy_indices)
            first_result, second_result = batch[indices[1]], batch[indices[2]]
            @test first_result.Z.values == second_result.Z.values
            @test first_result.Y.values == second_result.Y.values
            @test first_result.Z.values !== second_result.Z.values
            @test first_result.Y.values !== second_result.Y.values
            @test first_result.details.data.fem.primitive.Z_primitive !==
                  second_result.details.data.fem.primitive.Z_primitive
            @test run_directories[indices[1]] == run_directories[indices[2]]
            @test !isempty(details(first_result).data.timing)
            @test isempty(details(second_result).data.timing)
            @test details(second_result).data.fem.run.reused
            @test typeof(first_result) === typeof(second_result)
        end
        result = batch[first(default_indices)]
        run_directory = result.details.data.fem.run.run_directory
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
        @test result.details.data.fem.run.state === Base.get_extension(
            LineCableModels, :LineCableModelsGmshExt
        ).completed
        @test result.details.data.fem.run.getdp_invocations ==
              length(problem.frequencies)
        @test result.details.data.fem.run.completed_columns ==
              length(problem.frequencies) * length(result.details.data.fem.terminal_ids)
        @test run_directory !== nothing
        timing=result.details.data.timing
        @test result.details.data.fem.run.columns == result.details.data.fem.run.completed_columns
        @test result.details.data.fem.run.recovered_columns == 0
        @test !result.details.data.fem.run.reused
        @test keys(timing) == (:wall_seconds, :constraint_seconds, :assembly_seconds,
            :solve_seconds, :output_seconds, :worker_wall_seconds)
        @test timing.wall_seconds >= 0
        @test timing.worker_wall_seconds >= 0
        @test !haskey(result.details.data.fem, :timing)
        extension=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
        native=[extension._column_timing(extension._column_paths(
                run_directory,frequency,basis,false).timing,frequency,basis)
            for frequency in eachindex(problem.frequencies)
            for basis in eachindex(result.details.data.fem.terminal_ids)]
        for phase in (:constraint_seconds,:assembly_seconds,:solve_seconds,:output_seconds)
            @test getproperty(timing,phase) ≈ sum(getproperty(row,phase) for row in native)
        end
        expected_rows = 1 +
                        length(problem.frequencies) *
                        length(result.details.data.fem.terminal_ids)^2
        for filename in ("Z.tsv", "P.tsv")
            content = read(joinpath(run_directory, "results", filename), String)
            @test !occursin("\n\n", content)
            @test length(readlines(IOBuffer(content))) == expected_rows
        end

        # Same geometry and terminal conditions, now with explicitly requested
        # conduction plus polarization loss. The shared .pro is unchanged.
        lossy = batch[first(lossy_indices)]
        for (index, frequency) in pairs(problem.frequencies)
            expected_g = 2π / log(0.01 / 0.005) * (
                inv(dielectric.rho) +
                2π * frequency * 8.8541878128e-12 *
                dielectric.eps_r * dielectric.tan_delta)
            @test real(lossy.Y[1, 1, index]) ≈ expected_g rtol = 0.03
            @test abs(real(result.Y[1, 1, index])) < 0.03 * expected_g
            @test imag(lossy.Y[1, 1, index]) / (2π * frequency) ≈
                  expected_capacitance rtol = 0.02
        end
        # A fresh call can reconstruct the selected result from a completed run,
        # without starting Gmsh, touching historical files, or another solve.
        function snapshot_files(root)
            Dict(relpath(joinpath(directory, file),
                     root) => (bytes2hex(open(sha256, joinpath(directory, file))),
                     stat(joinpath(directory, file)).mtime)
            for (directory, _, files) in walkdir(root) for file in files)
        end
        before_reuse = snapshot_files(run_directory)
        @test !Bool(Gmsh.gmsh.is_initialized())
        repeated = compute(problem, selected[last(default_indices)];
            options = (;formulation_controls..., trace = true, resume_run_directory = run_directory))
        @test repeated.Z.values == result.Z.values
        @test repeated.Y.values == result.Y.values
        @test repeated.details.data.formulations.methods.earth_properties === nothing
        @test repeated.details.data.fem.run.run_directory == run_directory
        @test repeated.details.data.fem.run.reused
        @test isempty(repeated.details.data.timing)
        @test repeated.details.data.fem.inputs.getdp_identity.sha256 != ""
        @test !Bool(Gmsh.gmsh.is_initialized())
        @test snapshot_files(run_directory) == before_reuse
        batch_log = Test.TestLogger()
        reused_batch = with_logger(batch_log) do
            compute(problem, fill(selected[last(default_indices)], 3); options=(;
                formulation_controls..., trace=true, resume_run_directory=run_directory,
                verbosity=(default=0, progress=1)))
        end
        @test all(value->isempty(details(value).data.timing), reused_batch)
        @test all(value->Z(value) == Z(result) && Y(value) == Y(result), reused_batch)
        @test last(batch_log.logs).message == "FEM computation completed successfully"
        @test last(batch_log.logs).kwargs[:completed] == 3
        @test snapshot_files(run_directory) == before_reuse
        mktempdir() do temporary
            expected = joinpath(temporary, "expected.bin")
            serialize(expected, (result.Z.values, result.Y.values))
            code = raw"""
                using LineCableModels, Gmsh, JSON3, Serialization
                path, expected = ARGS
                problem = LineCableModels.ImportExport.deserialize_value(
                    JSON3.read(read(joinpath(path, "input", "problem.json"), String)))
                formulation = Formulation(:LineCableModelsFEM; earth_properties=nothing,
                    options=(ideal_transposition=false,))
                formulation_controls = (getdp_verbosity=0, gmsh_verbosity=0,
                        keep_run_directory=true)
                result = compute(problem, formulation; options=merge(formulation_controls, (trace=true, resume_run_directory=path)))
                Z, Y = deserialize(expected)
                @assert result.Z.values == Z && result.Y.values == Y
                @assert result.details.data.fem.run.run_directory == path
                @assert !Bool(Gmsh.gmsh.is_initialized())
                println("completed-run reuse verified in fresh Julia")
                """
            command = `$(Base.julia_cmd()) --startup-file=no --project=$(dirname(Base.active_project())) -e $code $run_directory $expected`
            @test occursin("completed-run reuse verified in fresh Julia", read(command, String))
            @test snapshot_files(run_directory) == before_reuse
        end
        # Timing does not alter native compatibility, including completed reuse.
        untimed = compute(problem, selected[last(default_indices)]; options=(;
            formulation_controls..., timing=false, trace=true, resume_run_directory=run_directory))
        @test !haskey(details(untimed).data, :timing)
        @test JSON3.write(details(untimed).data.fem.inputs) == JSON3.write(details(repeated).data.fem.inputs)
        @test Z(untimed) == Z(result) && Y(untimed) == Y(result)
        @test snapshot_files(run_directory) == before_reuse
        # Reopen this test-owned run with one missing column. The remaining
        # completed work must be recovered without becoming a fresh timing sample.
        state = JSON3.read(read(joinpath(run_directory, "run.json"), String), Dict{String,Any})
        state["state"] = "failed"
        write(joinpath(run_directory, "run.json"), JSON3.write(state))
        paths = extension._column_paths(run_directory, 1, 1, false)
        rm(paths.checkpoint)
        rm(paths.marker; force=true)
        for attempt in (isdir(joinpath(run_directory,"work")) ? readdir(joinpath(run_directory,"work");join=true) : String[])
            marker = extension._column_paths(attempt, 1, 1, false).marker
            rm(marker; force=true)
        end
        recovered = compute(problem, selected[last(default_indices)]; options=(;
            formulation_controls..., trace=true, resume_run_directory=run_directory))
        @test isempty(details(recovered).data.timing)
        @test details(recovered).data.fem.run.recovered_columns > 0
        @test details(recovered).data.fem.run.getdp_invocations ==
            details(result).data.fem.run.getdp_invocations + 1
        @test Z(recovered) == Z(result) && Y(recovered) == Y(result)
        rm(lossy.details.data.fem.run.run_directory; recursive = true, force = true)

        rm(run_directory; recursive = true, force = true)
    finally
        cd(fem_test_working_directory)
        rm(fem_test_runtime_directory;recursive=true)
    end

end
