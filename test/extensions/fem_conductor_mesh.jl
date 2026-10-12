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


@testitem "FEM / equal-size targets share fields and preserve contour counts" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    copper=Material(kind=:conductor,rho=1.72e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3,mu_r=1.)
    count,wire_radius,enclosure_radius,thickness=6,.005,.03,.005
    wires=Group(:wires,Region(:wire,Disk(wire_radius),copper);
        pattern=Ring(count;r=.02))
    design=build(CableDesign,"grouped-targets",
        Enclosure(:matrix,wires;primitive=Disk(enclosure_radius),fill=dielectric),
        Region(:insulation,Shell(thickness),dielectric))
    system=build(LineCableSystem,design,Pose2(0.,1.);connections=Dict(:wires=>1))
    problem=LineParametersProblem(system;temperature=20.,frequencies=[1e3],
        earth_props=homogeneous(rho=1000.))
    N.geometry(problem;mesh=true) do g,log
        regions=Int(only(g.parser.get_number("NumPhysicalRegions")))
        conductors=filter(r -> only(g.parser.get_number("FEMConductor_$r"))==1,0:regions-1)
        curves=Int.(g.parser.get_number("FEMLayerBoundaryCurves"))
        targets=Int.(g.parser.get_number("FEMLayerBoundaryTargetCurves"))
        contour_groups=Int(only(g.parser.get_number("FEMContourGroups")))
        layer_groups=Int(only(g.parser.get_number("FEMLayerGroups")))
        @test length(conductors)==count
        @test 0<contour_groups<length(conductors)
        @test layer_groups==2<length(targets)
        @test Set(targets)==Set(curves)
        @test Set(reduce(vcat,(Int.(g.parser.get_number("FEMContourGroupCurves_$i")) for i in 0:contour_groups-1)))==
            Set(reduce(vcat,(Int.(g.parser.get_number("FEMRegionCurves_$r")) for r in conductors)))
        @test Set(reduce(vcat,(Int.(g.parser.get_number("FEMLayerGroupCurves_$i")) for i in 0:layer_groups-1)))==Set(targets)

        # Expected counts come from the unchanged per-curve resolution criteria.
        tolerance=only(g.parser.get_number("ConductorGeometryTolerance"))
        circle_segments=12*2^max(0,ceil(Int,log2(sqrt(2π^2/(3tolerance))/12)))
        for region in conductors
            contour=Int.(g.parser.get_number("FEMRegionCurves_$region"))
            fractions=g.parser.get_number("FEMConductorArcFractions_$region")
            size=only(g.parser.get_number("RegionSize_$region"))
            expected=sum(max(1,ceil(Int,circle_segments*f-64eps(circle_segments*f)),
                ceil(Int,2π*wire_radius*f/size)) for f in fractions)
            nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in contour)))
            @test length(nodes)==expected
            println("Grouped targets: conductor region=$region, nodes=$(length(nodes)), expected=$expected")
        end
        radii=g.parser.get_number("FEMLayerBoundaryRadii")
        thicknesses=g.parser.get_number("FEMLayerBoundaryThicknesses")
        @test all(isapprox.(thicknesses,thickness))
        for radius in (enclosure_radius,enclosure_radius+thickness)
            boundary=curves[findall(r -> isapprox(r,radius),radii)]
            @test !isempty(boundary)
            expected=12*2^max(0,ceil(Int,log2(π*sqrt(radius/(3tolerance*thickness))/12)))
            nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in boundary)))
            @test length(nodes)==expected
            println("Grouped targets: radius=$radius, layer nodes=$(length(nodes)), expected=$expected")
        end
        println("Grouped targets: conductor groups=$contour_groups/$(length(conductors)), layer groups=$layer_groups/$(length(targets))")
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


@testitem "FEM / round conductor point contacts mesh with active boundary layers" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    copper=Material(kind=:conductor,rho=1.72e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3,mu_r=1.)
    inner,wire_radius=.02922,.0004
    outer=inner+2wire_radius
    # This contact angle reproduces an unrecovered edge with the former 0.8 cap.
    wire=Group(:wire,Region(:metal,Disk(wire_radius),copper);
        pattern=Ring(1;r=inner+wire_radius,φ0=π/2-.025))
    design=build(CableDesign,"round-point-contacts",
        Region(:centre,Disk(inner-.0003),dielectric),
        Region(:bedding,Annulus(inner-.0003,inner),dielectric),
        Enclosure(:matrix,wire;primitive=Annulus(inner,outer),fill=dielectric),
        terminal(:wire,Region(:screen,Shell(.0001),copper)),
        Region(:jacket,Shell(.0003),dielectric))
    system=build(LineCableSystem,design,Pose2(0.,-.9);connections=Dict(:wire=>1))
    problem=LineParametersProblem(system;temperature=20.,
        frequencies=[251188.6431509582],earth_props=homogeneous(rho=1000.))
    N.geometry(problem) do g,log
        g.option.set_number("General.Verbosity",3)
        number(name)=only(g.parser.get_number(name))
        @test number("FrequencyHz")==251188.6431509582
        regions=0:Int(number("NumPhysicalRegions"))-1
        round=only(filter(r -> number("FEMConductor_$r")==1 &&
            number("FEMConductorKind_$r")==1,regions))
        fields=g.model.mesh.field
        curves=g.parser.get_number("FEMRegionCurves_$round")
        layer=only(filter(f -> fields.get_type(f)=="BoundaryLayer" &&
            Set(fields.get_numbers(f,"CurvesList"))==Set(curves),fields.list()))
        surfaces=g.parser.get_number("FEMRegionSurfaces_$round")
        @test isempty(intersect(surfaces,fields.get_numbers(layer,"ExcludedSurfacesList")))
        g.model.mesh.generate(1)
        g.model.mesh.generate(2)
        @test !isempty(first(g.model.mesh.get_elements(2)))
        @test !any(line -> occursin("intersections in the 1D mesh",line) ||
            occursin("Edge not recovered",line),g.logger.get())
    end
end
