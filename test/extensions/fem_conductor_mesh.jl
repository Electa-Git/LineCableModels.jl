@testitem "Gmsh FEM / merged arc endpoints do not create a full circle" tags=[:extension, :fem] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    session = FEM._start_gmsh(0)
    try
        Gmsh.gmsh.model.add("merged-arc")
        registry = FEM.FEMLoopRegistry(1e-4)
        centre,radius = (0.00212,1.0),0.00212
        # Independent contact calculations differ in angle, but the physical
        # endpoints are the same within the registry's Float64 tolerance.
        FEM._register_circle_break!(registry,centre,radius,π+1.36e-12)
        curves = FEM._circle_arc_path!(registry,centre,radius,2.6,1.0)
        @test length(curves) == 2
        @test all(c -> begin a,b=registry.curve_points[abs(c)]; a!=b end,curves)
    finally
        FEM._finish_gmsh(session)
    end
end


@testitem "FEM / conductor contours resolve adjacent annular layers" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    copper=Material(kind=:conductor,rho=1.72e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3,mu_r=1.)
    for (radius,thickness,expected) in ((.0425,.001,384),(.01784,.0083,96))
        design=build(CableDesign,"insulated-wire",Stack(
            terminal(:core,Region(:metal,Disk(radius),copper)),
            Region(:insulation,Shell(thickness),dielectric)))
        system=build(LineCableSystem,[design],[Pose2(0.,1.)];
            connections=[Dict(:core=>1)],system_id="annular-layer-mesh",line_length=1.)
        problem=LineParametersProblem(system;temperature=20.,frequencies=[1e3],
            earth_props=homogeneous(rho=1000.))
        N.geometry(problem;mesh=true) do g,log
            curves=g.model.get_entities_for_physical_group(1,4001)
            nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in curves)))
            @test length(nodes)==expected
            @test g.parser.get_number("FEMConductorLayerThickness_0")≈fill(thickness,length(curves))
        end
    end
end


@testitem "FEM / circular non-conductor layer boundaries resolve layer gaps" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    aluminium=Material(kind=:conductor,rho=2.82e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3,mu_r=1.)
    semicon=Material(kind=:semicon,rho=1.,eps_r=10.,mu_r=1.)
    for (radius,shell,insulation) in ((.01784,0.,.0083),(.02,.001,.01),(.0425,0.,.001))
        layers=Any[terminal(:core,Region(:metal,Disk(radius),aluminium))]
        shell>0 && push!(layers,Region(:semicon,Shell(shell),semicon))
        push!(layers,Region(:insulation,Shell(insulation),dielectric))
        design=build(CableDesign,"layer-boundary-wire",Stack(layers...))
        system=build(LineCableSystem,[design],[Pose2(0.,1.)];
            connections=[Dict(:core=>1)],system_id="layer-boundary-mesh",line_length=1.)
        problem=LineParametersProblem(system;temperature=20.,frequencies=[1e3],
            earth_props=homogeneous(rho=1000.))
        # NativeFEMFixtures.geometry parses the freshly exported ONELAB bundle.
        N.geometry(problem;mesh=true) do g,log
            curves=Int.(g.parser.get_number("FEMLayerBoundaryCurves"))
            radii=g.parser.get_number("FEMLayerBoundaryRadii")
            thicknesses=g.parser.get_number("FEMLayerBoundaryThicknesses")
            targets=Int.(g.parser.get_number("FEMLayerBoundaryTargetCurves"))
            core_curves=g.model.get_entities_for_physical_group(1,4001)
            @test !isempty(curves)
            @test length(curves)==length(radii)==length(thicknesses)
            @test isempty(intersect(curves,core_curves))
            for i in eachindex(curves)
                @test !isempty(g.parser.get_number("FEMLayerBoundaryPoints_$(i-1)"))
                cables=g.parser.get_number("FEMLayerBoundaryCables_$(i-1)")
                @test isapprox(radii[i],radius+shell+insulation) ? cables==[0.] : isempty(cables)
            end
            if radius==.01784
                outer_curves=g.model.get_entities_for_physical_group(1,5001)
                outer_nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in outer_curves)))
                core_nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in core_curves)))
                @test all(isapprox.(radii,radius+insulation))
                @test all(isapprox.(thicknesses,insulation))
                @test Set(targets)==Set(outer_curves)
                @test length(outer_nodes)>=192
                @test length(core_nodes)==96
                println("Layer boundary: outer nodes=$(length(outer_nodes)), core nodes=$(length(core_nodes))")
            elseif shell>0
                shared=findall(r -> isapprox(r,radius+shell),radii)
                shared_curves=curves[shared]
                shared_nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in shared_curves)))
                @test !isempty(shared)
                @test all(isapprox.(thicknesses[shared],shell))
                @test all(in(targets),shared_curves)
                @test length(shared_nodes)>=384
                println("Layer boundary: semicon/insulation nodes=$(length(shared_nodes))")
            else
                @test all(isapprox.(radii,radius+insulation))
                @test all(isapprox.(thicknesses,insulation))
                @test isempty(targets)
                println("Layer boundary: already-fine outer contour, target curves=$(length(targets))")
            end
        end
    end
    N.geometry(N.problem()) do g,log
        @test isempty(g.parser.get_number("FEMLayerBoundaryCurves"))
        @test isempty(g.parser.get_number("FEMLayerBoundaryTargetCurves"))
    end
end


@testitem "FEM / native conductor contour and skin controls" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    values=[]
    for factor in (2.5,1.25)
        N.geometry(N.problem(;frequency=1e6);options=(mesh_size_factor=factor,)) do g,log
            a=Dict("first"=>only(g.parser.get_number("FEMFirst")),
                "bulk"=>only(g.parser.get_number("FEMBulk")))
            @test a["first"]>0 && a["bulk"]>0
            push!(values,a)
        end
    end
    @test values[2]["first"]<=values[1]["first"]
    @test values[2]["bulk"]<=values[1]["bulk"]
end
