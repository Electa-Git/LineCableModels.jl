@testitem "Gmsh FEM / polygon contacts preserve complete material partitions" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const DM = LineCableModels.DataModel
    const gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    polygons = (
        Polygon(((-0.003,-0.001), (0.0,-0.001), (0.0,0.001), (-0.003,0.001))),
        Polygon(((0.0,-0.001), (0.003,-0.001), (0.003,0.001), (0.0,0.001))),
        Polygon(((0.003,0.001), (0.005,0.001), (0.005,0.003), (0.003,0.003)))
    )
    body = terminal(:core, assembly((solid(copper, shape) for shape in polygons)...))
    design = build(CableDesign, "polygon-contacts",
        Enclosure(:matrix, body; primitive=Disk(0.01), fill=dielectric))
    session = FEM._start_gmsh(0)
    try
        for (index, pose) in enumerate((Pose2(0.0,-0.1), Pose2(0.06,-0.2,0.43)))
            system = build(LineCableSystem, design, pose; connections=Dict(:core=>1))
            problem = LineParametersProblem(system; frequencies=[50.0], earth_props=homogeneous(rho=100.0))
            model = FEM._resolved_fem_model(problem, LineCableModelsFEM())
            @test length(model.material_plans) == 2
            @test getproperty.(model.region_plans, :material_index) == [1,1,1,2]
            @test getproperty.(model.region_plans, :terminal_index) == [1,1,1,0]
            @test all(model.region_plans[i].shape == system.geometry[i].primitive for i in 1:3)
            geometry = FEM._build_geometry!(model, "polygon-contacts-$index")
            @test length(only(geometry.terminal_surfaces)) == 3
            FEM._configure_mesh!(model, geometry, only(model.mesh_plans))
            gmsh.model.mesh.generate(2)
            @test isnothing(FEM._inspect_loaded_mesh(model, "polygon-contacts"))
            mktempdir() do directory
                path = joinpath(directory, "polygon-contacts.msh")
                gmsh.write(path)
                @test isnothing(FEM._validate_mesh_file(model, path))
            end
            # A nonempty but incomplete filler must be rejected, including on
            # import. The former validator accepted exactly this situation.
            surface = first(geometry.material_surfaces[2])
            elements, _ = gmsh.model.mesh.get_elements_by_type(2, surface)
            @test length(elements) > 1
            gmsh.model.mesh.remove_elements(2, surface, elements[1:1])
            @test_throws LineCableModelsFEMError FEM._inspect_loaded_mesh(model, "partial-filler")
            mktempdir() do directory
                path = joinpath(directory, "partial.msh")
                gmsh.write(path)
                @test_throws LineCableModelsFEMError FEM._validate_mesh_file(model, path)
            end
        end
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / Milliken filler coverage and invariant scan topology" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const DM = LineCableModels.DataModel
    const gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor, rho=1.72e-8)
    core = milliken(copper; shape=Disk(0.33e-3),
        segment=Sector(span=pi/3, r_base=0.85e-3, r_back=3e-3, fillet=0.1e-3))
    design = build(CableDesign, "partition-milliken", terminal(:core, core))
    system = build(LineCableSystem, design, Pose2(0.0,-0.1); connections=Dict(:core=>1))
    problem = LineParametersProblem(system; frequencies=[0.1,50.0,1e7], earth_props=homogeneous(rho=100.0))
    model = FEM._resolved_fem_model(problem, LineCableModelsFEM())
    @test length(model.material_plans) == 2
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model, "partition-milliken")
        points = unique(gmsh.model.get_boundary(
            [(1,c) for c in only(geometry.terminal_curves)], false, false, true))
        coordinates = [gmsh.model.get_value(dim, tag, Float64[]) for (dim,tag) in points]
        entities = gmsh.model.get_entities()
        model_names = gmsh.model.list()
        for plan in model.mesh_plans
            FEM._update_exterior_mesh!(model, geometry, plan)
            @test length(gmsh.model.get_entities()) == length(entities)
            @test [gmsh.model.get_value(dim, tag, Float64[]) for (dim,tag) in points] == coordinates
            for (curves, radius) in ((geometry.inner_shell_curves,plan.domain_radius),
                                    (geometry.outer_curves,plan.shell_outer_radius))
                vertices = unique(gmsh.model.get_boundary([(1,c) for c in curves],false,false,true))
                @test all(vertices) do (dim,tag)
                    p = gmsh.model.get_value(dim,tag,Float64[])
                    hypot(p[1]-model.centre[1],p[2]) ≈ radius
                end
                @test all(curves) do curve
                    lower, upper = gmsh.model.get_parametrization_bounds(1, curve)
                    p = gmsh.model.get_value(1, curve, (lower + upper) / 2)
                    hypot(p[1]-model.centre[1],p[2]) ≈ radius
                end
            end
            FEM._configure_mesh!(model, geometry, plan)
            gmsh.model.mesh.generate(2)
            @test isnothing(FEM._inspect_loaded_mesh(model, "milliken-frequency-$(plan.frequency)"))
            @test gmsh.model.list() == model_names
        end
        @test length(geometry.material_surfaces[2]) > 1
    finally
        FEM._finish_gmsh(session)
    end
end
