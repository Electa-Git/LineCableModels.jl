@testitem "Gmsh FEM / rectangular boundaries and filled ring ownership" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    const DM = LineCableModels.DataModel
    const gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor,rho=1.72e-8)
    matrix = Material(kind=:semicon,rho=100.0,eps_r=10.0)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    limit, wire_radius = 0.6e-3, 0.05e-3
    body = stranded(copper;center=Disk(0.2e-3),shape=Rectangle(0.3e-3,0.1e-3),
        boundary=Disk(limit))
    core = terminal(:core,body)
    occupied = outer_radius(build(CableDesign,"core",core))
    designs = [build(CableDesign,"rectangular-bare",core),
        build(CableDesign,"rectangular-insulated",core,insulation(dielectric;t=0.1e-3))]
    for radius in (nothing,limit+wire_radius)
        ring = Group(:core,Region(:round_wire,Disk(wire_radius),copper);
            pattern=Ring(12;r=radius))
        filled = Enclosure(:ring_matrix,ring;
            primitive=Annulus(limit,limit+2wire_radius),fill=matrix)
        push!(designs,build(CableDesign,"rectangular-ring-$(radius)",core,filled,
            insulation(dielectric;t=0.1e-3)))
    end
    session = FEM._start_gmsh(0)
    try
        for design in designs
            system = build(LineCableSystem,design,Pose2(0.,-0.1);connections=Dict(:core=>1))
            problem = LineParametersProblem(system;frequencies=[50.0],earth_props=homogeneous(rho=100.0))
            model = FEM._resolved_fem_model(problem,LineCableModelsFEM())
            @test first(model.region_plans).shape isa Disk
            @test first(model.region_plans).shape.r == occupied
            @test sum(r -> area(r.shape),model.region_plans) ≈ area(design.geometry.outer)
            @test sum(r -> area(r.shape),filter(r -> r.terminal_index>0,model.region_plans)) ≈
                sum(area,filter(r -> r.source.material.kind===:conductor,design.geometry.regions))
            geometry = FEM._build_geometry!(model,design.cable_id)
            FEM._configure_mesh!(model,geometry,only(model.mesh_plans))
            gmsh.model.mesh.generate(2)
            @test isnothing(FEM._inspect_loaded_mesh(model,design.cable_id))
            mktempdir() do directory
                mesh = joinpath(directory,"model.msh")
                gmsh.write(mesh)
                @test isnothing(FEM._validate_mesh_file(model,mesh))
            end
        end
        # An explicitly declared thin layer is not absorbed into the conductor,
        # even in the interval between the former 2 ppm / 5 ppm cutoffs.
        for fraction in (1e-6,3e-6,6e-6)
            design = build(CableDesign,"explicit-thin-wrap",Enclosure(:paper,core;
                primitive=Disk(occupied/sqrt(1-fraction)),fill=dielectric))
            system = build(LineCableSystem,design,Pose2(0.,-0.1);connections=Dict(:core=>1))
            problem = LineParametersProblem(system;frequencies=[50.],earth_props=homogeneous(rho=100.))
            model = FEM._resolved_fem_model(problem,LineCableModelsFEM())
            metal, fill_region = model.region_plans
            @test metal.shape.r == fill_region.shape.ri == occupied
            @test area(fill_region.shape)>0
            @test metal.shape.r < fill_region.shape.ro
            @test FEM._build_geometry!(model,"explicit-thin-$fraction") isa FEM.FEMGeometry
        end
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / rectangular core capacitance uses occupied radius" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh
    copper = Material(kind=:conductor,rho=1.72e-8)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    body = stranded(copper;center=Disk(0.2e-3),shape=Rectangle(0.3e-3,0.1e-3),
        boundary=Disk(0.6e-3))
    design = build(CableDesign,"rectangular-coaxial",
        terminal(:core,body),insulation(dielectric;t=0.2e-3),
        terminal(:sheath,Region(:screen,Shell(0.05e-3),copper)),
        insulation(dielectric;t=0.1e-3))
    occupied = design.geometry.regions[end-2].primitive.ri
    system = build(LineCableSystem,design,Pose2(0.,-0.1);connections=Dict(:core=>1,:sheath=>0))
    problem = LineParametersProblem(system;frequencies=[50.,1000.],earth_props=homogeneous(rho=100.))
    result = compute(problem,LineCableModelsFEM();options=(gmsh_verbosity=0,getdp_verbosity=0))
    expected = 2pi*8.8541878128e-12*2.3/log((occupied+0.2e-3)/occupied)
    @test all(isfinite,result.Z)
    @test all(isfinite,result.Y)
    @test all(isapprox(imag(result.Y[1,1,i])/(2pi*f),expected;rtol=0.02)
        for (i,f) in enumerate(result.f))
end
