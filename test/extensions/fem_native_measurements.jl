@testitem "Gmsh FEM / native line integration and complex matrix inversion" tags=[:extension, :fem_numerical] begin
    using Gmsh, LinearAlgebra
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    execution = computation_options(LineCableModelsFEM, ComputationOptions())
    getdp = FEM._getdp_selection(execution).path
    fixtures = joinpath(@__DIR__, "..", "fixtures", "data", "fem", "native_getdp")
    mktempdir() do root
        for name in ("model.geo", "model.pro", "inverse.pro")
            cp(joinpath(fixtures, name), joinpath(root, name))
        end
        session = FEM._start_gmsh(0)
        try
            gmsh.open(joinpath(root, "model.geo"))
            gmsh.model.mesh.generate(2)
            gmsh.write(joinpath(root, "model.msh"))
        finally
            FEM._finish_gmsh(session)
        end
        # Analytic field (1+2j)*(1-y,x,0), represented exactly by edge elements.
        # The measurement curve is embedded in the volume mesh.
        run(Cmd(`$getdp model.pro -msh model.msh -solve Solve -v 2`; dir=root))
        values(name) = parse.(Float64, split(read(joinpath(root,name),String)))
        @test values("open.txt")[2:3] ≈ [0.24,0.48] atol=1e-13 rtol=0
        @test values("loop.txt")[2:3] ≈ [2.,4.] atol=1e-13 rtol=0
        @test values("offmesh.txt")[9:14] ≈ [.56,.31,0.,1.12,.62,0.] atol=1e-13 rtol=0
        # Arbitrary-size native algebra: four unknowns, no Tensor/Inv size limit.
        run(Cmd(`$getdp inverse.pro -msh model.msh -solve Invert -v 2`; dir=root))
        P = [complex(5*(i==j)+.1i+.2j, .15i-.07j) for i in 1:4, j in 1:4]
        Y = zeros(ComplexF64,4,4)
        for line in eachline(joinpath(root,"arbitrary.txt"))
            i,j,re,im = parse.(Float64,split(line))
            Y[Int(i),Int(j)] = complex(re,im)
        end
        @test Y ≈ inv(P) atol=1e-13 rtol=0
        @test maximum(abs, P*Y-I) < 1e-13
    end
end

@testitem "Gmsh FEM / conforming path through a coated rotated ellipse" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    copper = Material(kind=:conductor,rho=2e-8)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    wire = build(CableDesign,"ellipse-path",Enclosure(:insulation,
        terminal(:core,solid(copper,Ellipse(.02,.01,Pose2(0.,0.,π/7))));
        primitive=Disk(.04),fill=dielectric))
    system = build(LineCableSystem,wire,Pose2(.12,-.2,.21);connections=Dict(:core=>1))
    problem = LineParametersProblem(system;frequencies=[1e4],earth_props=homogeneous(rho=100.))
    model = FEM._resolved_fem_model(problem,LineCableModelsFEM(),
        computation_options(LineCableModelsFEM,ComputationOptions(pml_layers=8)))
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"ellipse-path")
        FEM._configure_mesh!(model,geometry)
        gmsh.model.mesh.generate(2)
        @test FEM._inspect_loaded_mesh(model,"ellipse-path") === nothing
        passive = [s for (i,m) in enumerate(model.material_plans) if m.kind !== :conductor
            for s in geometry.material_surfaces[i]]
        @test any(!isempty(gmsh.model.mesh.get_embedded(2,s)) for s in passive)
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / direct edge trace across a field discontinuity" tags=[:extension,:fem_numerical] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    fixtures = joinpath(@__DIR__,"..","fixtures","data","fem","native_getdp")
    mktempdir() do root
        for name in ("support.geo","support.pro")
            cp(joinpath(fixtures,name),joinpath(root,name))
        end
        session = FEM._start_gmsh(0)
        try
            gmsh.open(joinpath(root,"support.geo"))
            gmsh.model.mesh.generate(2)
            gmsh.write(joinpath(root,"support.msh"))
            values(name) = parse.(Float64,split(read(joinpath(root,name),String)))
            # Integral of (1+2j)*(2 below y=.37, 6 above it). The field is
            # discontinuous; native trace integration resolves its pieces.
            for points in (1,2,4,10,20)
                run(Cmd(`$getdp support.pro -msh support.msh -solve Solve -setnumber LinePoints $points -v 2`;dir=root))
                forward = complex(values("forward.txt")[2:3]...)
                reverse = complex(values("reverse.txt")[2:3]...)
                upwards = complex(values("upwards.txt")[2:3]...)
                @test forward ≈ 4.52+9.04im atol=1e-13 rtol=0
                @test reverse ≈ -forward atol=1e-13 rtol=0
                @test upwards ≈ forward atol=1e-13 rtol=0
            end
        finally
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / explicit voltage curves and reference points" tags=[:extension] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    wire = build(CableDesign, "native-path", terminal(:core,
        core(Material(kind=:conductor,rho=2e-8);r=.0425)))
    system = build(LineCableSystem,[wire,wire],[(0.,1.),(1.,-1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[1e4],
        earth_props=homogeneous(rho=100.,eps_r=1.))
    model = FEM._resolved_fem_model(problem,LineCableModelsFEM(),
        computation_options(LineCableModelsFEM,ComputationOptions(pml_layers=8)))
    session = FEM._start_gmsh(0)
    try
        geometry = FEM._build_geometry!(model,"native-voltage-curves")
        for i in 1:2
            curves = gmsh.model.get_entities_for_physical_group(1,model.tags.voltage_path_base+i)
            ref = only(gmsh.model.get_entities_for_physical_group(0,model.tags.voltage_reference_base+i))
            from = gmsh.model.get_value(0,ref,Float64[])
            endpoints = gmsh.model.get_boundary([(1,c) for c in curves],true,false,false)
            coords = [gmsh.model.get_value(0,last(p),Float64[]) for p in endpoints]
            to = only(filter(!=(from),coords))
            @test from[1] ≈ to[1] atol=1e-14
            @test to[1] ≈ i-1 atol=1e-14
            @test to[2] ≈ (i==1 ? 1. : -1.)-.0425
            plan = only(model.mesh_plans)
            @test from[2] ≈ (i==1 ? 0. : -plan.domain_halfwidth-plan.pml_thickness[3])
            @test to[2] > from[2]
        end
        FEM._configure_mesh!(model,geometry)
        gmsh.model.mesh.generate(2)
        @test FEM._inspect_loaded_mesh(model,"conforming") === nothing
        gmsh.model.mesh.set_order(2)
        @test FEM._inspect_loaded_mesh(model,"conforming quadratic mesh") === nothing
        gmsh.model.mesh.set_order(1)
        # Same coordinates with unrelated node identities reproduce the old
        # detached-line error. Such a mesh must fail before solver assembly.
        curve = first(gmsh.model.get_entities_for_physical_group(1,model.tags.voltage_path_base+1))
        elements,nodes = gmsh.model.mesh.get_elements_by_type(1,curve)
        tags = maximum(first(gmsh.model.mesh.get_nodes())) .+ UInt64[1,2]
        coordinates = [gmsh.model.mesh.get_node(n)[1] for n in nodes[1:2]]
        gmsh.model.mesh.remove_elements(1,curve,elements[1:1])
        gmsh.model.mesh.add_nodes(1,curve,tags,reduce(vcat,coordinates))
        gmsh.model.mesh.add_elements_by_type(curve,1,elements[1:1],tags)
        @test_throws LineCableModelsFEMError FEM._inspect_loaded_mesh(model,"detached")
    finally
        FEM._finish_gmsh(session)
    end
end
