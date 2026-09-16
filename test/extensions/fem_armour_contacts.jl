@testitem "Gmsh FEM / touching armour wires retain every filler face" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor, rho=1.72e-8)
    bedding = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    matrix = Material(kind=:insulator, rho=Inf, eps_r=1.0)
    session = FEM._start_gmsh(0)
    try
        for count in (6, 8, 68), pose in (Pose2(0.0,-0.1), Pose2(-0.5,-1.0,0.37))
            wire_radius = 0.00291
            radius = wire_radius / sinpi(1/count)
            inner, outer = radius-wire_radius, radius+wire_radius
            ring = Group(:armour, Region(:wire, Disk(wire_radius), copper);
                pattern=Ring(count; r=radius))
            design = build(CableDesign, "touching-armour",
                Region(:bedding, Disk(inner), bedding),
                Enclosure(:matrix, ring; primitive=Annulus(inner,outer), fill=matrix),
                Region(:jacket, Shell(0.002), bedding))
            system = build(LineCableSystem, design, pose; connections=Dict(:armour=>1))
            problem = LineParametersProblem(system; frequencies=[0.1,50.0],
                earth_props=homogeneous(rho=100.0))
            model = FEM._resolved_fem_model(problem, LineCableModelsFEM())
            fill_index = only(findall(m -> m.field === :matrix_fill, model.material_plans))
            geometry = FEM._build_geometry!(model, "touching-armour-$count")
            surfaces = geometry.material_surfaces[fill_index]
            # Every pair of neighbouring wires encloses an inner and an outer
            # filler lobe. A point contact cannot join them into one CAD face.
            @test length(surfaces) == 2count
            for surface in surfaces
                incidence = Dict{Int,Int}()
                for (_,curve) in gmsh.model.get_boundary([(2,surface)], false, false, false)
                    for (_,point) in gmsh.model.get_boundary([(1,curve)], false, false, false)
                        incidence[point] = get(incidence,point,0)+1
                    end
                end
                @test all(==(2), values(incidence))
            end
            for plan in model.mesh_plans
                FEM._update_exterior_mesh!(model, geometry, plan)
                FEM._configure_mesh!(model, geometry, plan)
                gmsh.model.mesh.generate(2)
                @test isnothing(FEM._inspect_loaded_mesh(model, "touching-armour"))
                mktempdir() do directory
                    mesh = joinpath(directory,"armour.msh")
                    gmsh.write(mesh)
                    @test isnothing(FEM._validate_mesh_file(model,mesh))
                end
            end
        end
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / partial and separated armour rings preserve physical gaps" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor, rho=1.72e-8)
    bedding = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    matrix = Material(kind=:insulator, rho=Inf, eps_r=1.0)
    wire_radius = 0.00291
    session = FEM._start_gmsh(0)
    try
        for (count,span,gap) in ((7,6pi/4,0.0), (2,pi/4,0.0),
                                (8,2pi,1e-6), (1,2pi,0.0)), padding in (0.0,0.001)
            radius = (wire_radius+gap/2)/sinpi(1/8)
            inner,outer = radius-wire_radius,radius+wire_radius+padding
            ring = Group(:armour,Region(:wire,Disk(wire_radius),copper);
                pattern=Ring(count;r=radius,span,φ0=0.23))
            design = build(CableDesign,"partial-armour",
                Region(:bedding,Disk(inner),bedding),
                Enclosure(:matrix,ring;primitive=Annulus(inner,outer),fill=matrix),
                Region(:jacket,Shell(0.002),bedding))
            system = build(LineCableSystem,design,Pose2(0.0,-0.1);connections=Dict(:armour=>1))
            problem = LineParametersProblem(system;frequencies=[50.0],
                earth_props=homogeneous(rho=100.0))
            model = FEM._resolved_fem_model(problem,LineCableModelsFEM())
            fill_index = only(findall(m -> m.field === :matrix_fill,model.material_plans))
            geometry = FEM._build_geometry!(model,"partial-armour")
            expected = (count == 1 || gap > 0 ? count : 2count-1) + (padding > 0)
            @test length(geometry.material_surfaces[fill_index]) == expected
            @test sum(r -> area(r.shape),filter(r -> r.material_index==fill_index,
                model.region_plans)) ≈ pi*(outer^2-inner^2)-count*pi*wire_radius^2
            FEM._configure_mesh!(model,geometry,only(model.mesh_plans))
            gmsh.model.mesh.generate(2)
            @test isnothing(FEM._inspect_loaded_mesh(model,"partial-armour"))
        end
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / touching armour filler reaches GetDP" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh
    copper = Material(kind=:conductor,rho=1.72e-8)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    matrix = Material(kind=:insulator,rho=Inf,eps_r=1.0)
    count,wire_radius = 8,0.00291
    radius = wire_radius/sinpi(1/count)
    ring = Group(:armour,Region(:wire,Disk(wire_radius),copper);
        pattern=Ring(count;r=radius))
    design = build(CableDesign,"armour-contact-solve",
        Region(:bedding,Disk(radius-wire_radius),dielectric),
        Enclosure(:matrix,ring;primitive=Annulus(radius-wire_radius,radius+wire_radius),fill=matrix),
        Region(:jacket,Shell(0.002),dielectric))
    system = build(LineCableSystem,design,Pose2(0.0,-0.1);connections=Dict(:armour=>1))
    problem = LineParametersProblem(system;frequencies=[0.1,50.0],earth_props=homogeneous(rho=100.0))
    result = compute(problem,LineCableModelsFEM();options=(gmsh_verbosity=0,getdp_verbosity=0))
    @test result.f == problem.frequencies
    @test all(isfinite,result.Z)
    @test all(isfinite,result.Y)
    @test all(real(z)>0 for z in result.Z)
end
