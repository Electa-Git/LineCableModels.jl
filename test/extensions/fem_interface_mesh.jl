@testitem "FEM / editable interface footprint preserves terminal contours" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    p=N.problem(;frequency=1e6,rho=.1,positions=[(0.,1.),(2.,1.)])
    records=[]
    for factor in (1.,12.)
        N.geometry(p;options=(overrides=(InterfaceRefinementFactor=factor,),),mesh=true) do g,log
            @test only(g.parser.get_number("InterfaceRefinementFactor"))==factor
            contours=[sum(length(block) for c in g.model.get_entities_for_physical_group(1,4000+i) for block in g.model.mesh.get_elements(1,c)[2]) for i in 1:2]
            push!(records,contours)
            @test !any(l->occursin("closer than the geometrical tolerance",l),log)
        end
    end
    @test records[1]==records[2]
end


@testitem "FEM / active interface footprint materially refines the mesh" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    p=N.problem(;frequency=17_782_794.100389227,rho=1000.,eps_r=1.,
        radius=.0425,positions=[(0.,1.),(1.,-1.)])
    records=[]
    for factor in (1.,12.)
        N.geometry(p;options=(mesh_size_factor=1.,overrides=(InterfaceRefinementFactor=factor,),),mesh=true) do g,log
            @test only(g.parser.get_number("InterfaceRefinementFactor"))==factor
            @test only(g.parser.get_number("FEMEarthLayerActive"))==0
            if factor==1.
                @test only(g.parser.get_number("MeshWaveEarth"))<only(g.parser.get_number("MeshRemoteEarth"))
            end
            contours=[sum(length(block) for c in g.model.get_entities_for_physical_group(1,4000+i) for block in g.model.mesh.get_elements(1,c)[2]) for i in 1:2]
            push!(records,(nodes=length(first(g.model.mesh.get_nodes())),contours=contours))
            @test !any(l->occursin("closer than the geometrical tolerance",l),log)
        end
    end
    # A 1% increase is material, well above the observed <0.1% node-count noise.
    @test records[2].nodes>=ceil(Int,1.01 * records[1].nodes)
    @test records[1].contours==records[2].contours
end
