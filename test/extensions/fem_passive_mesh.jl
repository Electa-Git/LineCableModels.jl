@testitem "FEM / passive interiors use their own layer targets" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    copper=Material(kind=:conductor,rho=1.72e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3,mu_r=1.)
    semicon=Material(kind=:semicon,rho=1.,eps_r=10.,mu_r=1.)
    radius,thickness,outer_thickness=.01,.01,.0003
    for outer in (false,true)
        layers=Any[terminal(:core,Region(:metal,Disk(radius),copper)),
            Region(:insulation,Shell(thickness),dielectric)]
        outer && push!(layers,Region(:semicon,Shell(outer_thickness),semicon))
        design=build(CableDesign,"passive-layer-wire",Stack(layers...))
        system=build(LineCableSystem,[design],[Pose2(0.,1.)];
            connections=[Dict(:core=>1)],system_id="passive-layer-mesh",line_length=1.)
        problem=LineParametersProblem(system;temperature=20.,frequencies=[1e3],
            earth_props=homogeneous(rho=1000.))
        N.geometry(problem;mesh=true) do g,log
            number(name)=only(g.parser.get_number(name))
            numbers(name)=Int.(g.parser.get_number(name))
            regions=0:Int(number("NumPhysicalRegions"))-1
            passive=filter(r -> number("FEMPassive_$r")==1,regions)
            fields=g.model.mesh.field
            extensions=filter(f -> fields.get_type(f)=="Extend",fields.list())
            @test length(extensions)==length(passive)==(outer ? 2 : 1)
            for region in passive
                curves=numbers("FEMRegionCurves_$region")
                surfaces=numbers("FEMRegionSurfaces_$region")
                extension=only(filter(f -> Int.(fields.get_numbers(f,"CurvesList"))==curves,extensions))
                target=number("RegionSize_$region")
                @test fields.get_number(extension,"SizeMax")==target
                @test fields.get_number(extension,"DistMax")==target/number("MeshGrowthSlope")
                restriction=only(filter(fields.list()) do f
                    fields.get_type(f)=="Restrict" && fields.get_number(f,"InField")==extension
                end)
                @test Int.(fields.get_numbers(restriction,"SurfacesList"))==surfaces
                @test restriction in numbers("FEMBackground")
                for medium in 0:1
                    @test extension ∉ numbers("FEMMeasurementSources_$medium")
                    @test restriction ∉ numbers("FEMMeasurementSources_$medium")
                end
            end
            if outer
                # In the middle half of the insulation, distance from either
                # boundary is at least thickness/4. Its own native extension
                # permits growth over that distance rather than the outer cap.
                thick_target=number("RegionSize_1")
                thin_target=number("RegionSize_2")
                factor=number("MeshSizeFactor")
                @test thick_target≈factor*thickness/2
                @test thin_target≈factor*2π*(radius+thickness+outer_thickness)/number("FEMCircleSegments")
                threshold=min(thick_target,number("MeshGrowthSlope")*thickness/4)
                @test thin_target>factor*outer_thickness/2
                tags,xyz,_=g.model.mesh.get_nodes()
                coordinates=Dict(tag=>xyz[3i-2:3i] for (i,tag) in enumerate(tags))
                edges=Float64[]
                for surface in numbers("FEMRegionSurfaces_1")
                    kinds,_,blocks=g.model.mesh.get_elements(2,surface)
                    @test all(==(2),kinds)
                    for nodes in blocks, i in 1:3:length(nodes)
                        points=[coordinates[nodes[i+j]] for j in 0:2]
                        centre=sum(points)/3
                        r=hypot(centre[1],centre[2]-1)
                        radius+thickness/4<r<radius+3thickness/4 || continue
                        append!(edges,[sqrt(sum(abs2,points[j]-points[mod1(j+1,3)])) for j in 1:3])
                    end
                end
                @test !isempty(edges)
                @test maximum(edges)>threshold
                println("Passive interior: max middle-layer edge=$(maximum(edges)) m, derived threshold=$threshold m, outer target=$thin_target m")
            else
                # Structural equivalence to the former rule: with one passive
                # region the cable-wide boundary/surfaces and cap are identical.
                region=only(passive)
                surfaces=numbers("FEMRegionSurfaces_$region")
                @test numbers("FEMRegionCurves_$region")==N.FEM._entity_boundary(surfaces)
                @test number("CableSize_0")==number("RegionSize_$region")
                @test numbers("CableSizingRegions_0")==[region]
            end
        end
    end
end

@testitem "FEM / thin annuli separate normal and tangential resolution" tags=[:extension, :fem] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    copper=Material(kind=:conductor,rho=1.72e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=Inf,eps_r=2.3,mu_r=1.)
    counts=Int[]
    for thickness in (25e-6,50e-6,100e-6)
        radius=.01
        design=build(CableDesign,"thin-annular-wire",Stack(
            terminal(:core,Region(:metal,Disk(radius),copper)),
            Region(:face,Shell(thickness),dielectric),
            Region(:jacket,Shell(.003),dielectric)))
        system=build(LineCableSystem,design,Pose2(0.,1.);connections=Dict(:core=>1))
        problem=LineParametersProblem(system;temperature=20.,frequencies=[1e3],
            earth_props=homogeneous(rho=1000.))
        N.geometry(problem;mesh=true) do g,log
            numbers(name)=Int.(g.parser.get_number(name))
            surfaces=numbers("FEMRegionSurfaces_1")
            nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(2,s,true)) for s in surfaces)))
            push!(counts,length(nodes))
            tags,xyz,_=g.model.mesh.get_nodes()
            coordinates=Dict(tag=>xyz[3i-2:3i] for (i,tag) in enumerate(tags))
            # Two normal elements, even when the tangential edges are much longer.
            radii=sort!([hypot(coordinates[n][1],coordinates[n][2]-1) for n in nodes])
            @test all(r -> any(isapprox(r,expected;atol=1e-12,rtol=0)
                for expected in (radius,radius+thickness/2,radius+thickness)),radii)
            @test all(expected -> any(isapprox(r,expected;atol=1e-12,rtol=0)
                for r in radii),(radius,radius+thickness/2,radius+thickness))
            area=0.
            for surface in surfaces
                kinds,_,blocks=g.model.mesh.get_elements(2,surface)
                @test all(==(2),kinds)
                for block in blocks, i in 1:3:length(block)
                    a,b,c=[coordinates[block[i+j]] for j in 0:2]
                    area+=abs((b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1]))/2
                end
            end
            tolerance=only(g.parser.get_number("ConductorGeometryTolerance"))
            @test abs(area/(π*((radius+thickness)^2-radius^2))-1)<tolerance
            # At corresponding edge midpoints, the polygonal gap is t*cos(π/N).
            curves=numbers("FEMRegionCurves_1")
            ring_nodes=unique(reduce(vcat,(first(g.model.mesh.get_nodes(1,c,true)) for c in curves)))
            angular=length(ring_nodes)÷2
            @test 1-cos(π/angular)<tolerance
            @test length(nodes)==3angular
        end
    end
    @test all(==(first(counts)),counts)
end
