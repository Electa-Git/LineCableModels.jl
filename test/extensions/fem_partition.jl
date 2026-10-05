@testitem "Gmsh FEM / polygon contacts preserve complete material partitions" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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
            geometry = FEM._build_physical_geometry!(model, "polygon-contacts-$index")
            @test length(only(geometry.terminal_surfaces)) == 3

        end
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / Milliken filler coverage and invariant scan topology" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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
        geometry = FEM._build_physical_geometry!(model, "partition-milliken")
        @test length(geometry.material_surfaces[2]) > 1
    finally
        FEM._finish_gmsh(session)
    end
end
