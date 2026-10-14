@testitem "Gmsh FEM / polygon contacts preserve complete material partitions" tags=[:extension, :fem] begin
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
            model = FEM._resolved_fem_model(problem, Formulation(:fem))
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

@testitem "Gmsh FEM / Milliken filler coverage and invariant scan topology" tags=[:extension, :fem] begin
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
    model = FEM._resolved_fem_model(problem, Formulation(:fem))
    @test length(model.material_plans) == 2
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_physical_geometry!(model, "partition-milliken")
        @test length(geometry.material_surfaces[2]) > 1
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / native polylines retain polygon CAD and resolve its area" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N = NativeFEMFixtures
    radius = 0.003
    points = [(radius*cos(2π*i/128), radius*sin(2π*i/128)) for i in 0:127]
    source_area = 128radius^2*sin(2π/128)/2
    copper = Material(kind=:conductor, rho=1.72e-8)
    design = build(CableDesign, "polygon-wire", terminal(:core, solid(copper, Polygon(points))))
    system = build(LineCableSystem, design, Pose2(0.,-0.1); connections=Dict(:core=>1))
    problem = LineParametersProblem(system; frequencies=[50.], earth_props=homogeneous(rho=100.))
    for tolerance in (1e-3, 1e-5)
        N.geometry(problem; options=(overrides=(ConductorGeometryTolerance=tolerance,),), mesh=true) do g, log
            curves = Int.(g.parser.get_number("FEMPolygonCurves"))
            @test length(curves) == 4
            # Every retained source vertex is on the native CAD, even though
            # it is no longer an independent curve endpoint.
            controls = Tuple{Float64,Float64}[]
            for curve in curves
                lo, hi = g.model.get_parametrization_bounds(1,curve)
                xyz = g.model.get_value(1,curve,collect(range(only(lo),only(hi);length=33)))
                append!(controls, [(xyz[i],xyz[i+1]) for i in 1:3:length(xyz)])
            end
            for point in points
                distance = minimum(q -> hypot(q[1]-point[1], q[2]-(point[2]-0.1)), controls)
                @test distance < 1e-12
            end
            tags, xyz, _ = g.model.mesh.get_nodes()
            positions = Dict(tag => (xyz[3i-2],xyz[3i-1]) for (i,tag) in enumerate(tags))
            area = 0.0
            for surface in g.model.get_entities_for_physical_group(2,3001)
                types, _, elements = g.model.mesh.get_elements(2,surface)
                @test types == [2]
                for nodes in elements, i in 1:3:length(nodes)
                    a,b,c = (positions[nodes[i+j]] for j in 0:2)
                    area += abs((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))/2
                end
            end
            @test abs(area/source_area-1) <= tolerance
        end
    end
end

@testitem "Gmsh FEM / polygon polylines preserve shared curved seams" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N = NativeFEMFixtures
    copper = Material(kind=:conductor, rho=1.72e-8)
    seam = [(0.,0.), (.001,.0001), (.002,.00015), (.003,0.)]
    lower = Polygon(vcat(reverse(seam), [(0.,-.002), (.003,-.002)]))
    upper = Polygon(vcat(seam, [(.003,.002), (0.,.002)]))
    design = build(CableDesign, "polygon-seams", terminal(:core,
        assembly(solid(copper,lower),solid(copper,upper))))
    system = build(LineCableSystem, design, Pose2(0.,-0.1); connections=Dict(:core=>1))
    problem = LineParametersProblem(system; frequencies=[50.], earth_props=homogeneous(rho=100.))
    N.geometry(problem; mesh=true) do g, log
        surfaces = g.model.get_entities_for_physical_group(2,3001)
        @test length(surfaces) == 2
        boundaries = [Set(abs(tag) for (_,tag) in
            g.model.get_boundary([(2,s)],false,true,false)) for s in surfaces]
        @test length(intersect(boundaries...)) == length(seam)-1
    end
end

@testitem "Gmsh FEM / managed polygon-only material area coverage" tags=[:extension, :integration, :fem_numerical, :fem] setup=[TemporaryFEMRuntime] begin
    using Gmsh, LinearAlgebra
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    g = Gmsh.gmsh
    cd(fem_test_runtime_directory)
    try
        copper = Material(kind=:conductor, rho=1.72e-8)
        outline = Polygon(Tuple((0.003cospi(i/6),0.003sinpi(i/6)) for i in 0:11))
        wire = build(CableDesign, "polygon-only", terminal(:core, solid(copper,outline)))
        system = build(LineCableSystem, [wire,wire], [(0.,0.1),(0.2,-0.1)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
        problem = LineParametersProblem(system; frequencies=[50.],
            earth_props=homogeneous(rho=100.,eps_r=10.))
        form = Formulation(:fem; options=(
            reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        result = compute(problem,form; options=(
            overrides=(PmlSideLayers=8,PmlTopLayers=6,PmlBottomLayers=4,
                PmlSideGrading=3.,PmlTopGrading=2.,PmlBottomGrading=1.,
                MumpsOrdering=0,PetscPrealloc=256),
            frequency_workers=1,solver_threads=1,keep_run_directory=true,
            gmsh_verbosity=0,getdp_verbosity=0))
        @test all(isfinite,Z(result))
        @test all(isfinite,Y(result))
        model = FEM._resolved_fem_model(problem,form)
        material = only(filter(m -> m.kind === :conductor,model.material_plans))
        regions = filter(r -> model.material_plans[r.material_index] === material,model.region_plans)
        @test length(regions) == 2
        @test all(r -> r.shape isa LineCableModels.DataModel.Polygon,regions)
        mesh = joinpath(details(result).data.fem.run.run_directory,"mesh","frequency_0001.msh")
        session = FEM._start_gmsh(0)
        try
            g.open(mesh)
            tags,xyz,_ = g.model.mesh.get_nodes()
            points = Dict(tag => (xyz[3i-2],xyz[3i-1]) for (i,tag) in enumerate(tags))
            actual = 0.0
            h2 = 0.0
            for surface in g.model.get_entities_for_physical_group(2,material.physical_tag)
                kinds,_,blocks = g.model.mesh.get_elements(2,surface)
                @test kinds == [2]
                for nodes in blocks, i in 1:3:length(nodes)
                    a,b,c = (points[nodes[i+j]] for j in 0:2)
                    actual += abs((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))/2
                end
                for (_,curve) in g.model.get_boundary([(2,surface)],false,false,false)
                    kinds,_,blocks = g.model.mesh.get_elements(1,curve)
                    @test kinds == [1]
                    for nodes in blocks, i in 1:2:length(nodes)
                        a,b = points[nodes[i]],points[nodes[i+1]]
                        h2 = max(h2,(a[1]-b[1])^2+(a[2]-b[2])^2)
                    end
                end
            end
            expected = 5.4e-5 # Two regular dodecagons of radius 3 mm.
            @test sum(r -> area(r.shape),regions) ≈ expected
            # This fixture must exceed the former zero-allowance assertion.
            @test abs(actual-expected) > 1e-8expected
            @test abs(actual-expected) <= 1e-3expected
            allowance = sum(r -> FEM._mesh_area_allowance(r.shape,h2),regions)
            @test abs(actual-expected) <= 1e-8expected + allowance
            @test FEM._mesh_area_allowance(Rectangle(0.003,0.002),h2) == 0.0
            println("Polygon coverage: actual=",actual," expected=",expected," h2=",h2," allowance=",allowance)
        finally
            FEM._finish_gmsh(session)
        end
    finally
        cd(fem_test_working_directory)
        rm(fem_test_runtime_directory;recursive=true)
    end
end
