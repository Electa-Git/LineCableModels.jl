@testitem "FEM / native PML cell families and quadrature transport" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    inventories=[]
    for (family,physical,quadrature) in ((:triangle,3,9),(:quadrangle,3,9),(:quadrangle,12,16))
        N.geometry(N.problem(;frequency=1e6,rho=100.);
                options=(overrides=(PmlQuadrangles=Int(family===:quadrangle),PhysicalVolumeQuadrature=physical,PmlQuadrature=quadrature),),mesh=true) do g,log
            @test only(g.parser.get_number("PhysicalVolumeQuadrature"))==physical
            @test only(g.parser.get_number("PmlQuadrature"))==quadrature
            pml=g.model.get_entities_for_physical_group(2,1005)
            kinds=[g.model.mesh.get_elements(2,s)[1] for s in pml]
            @test all(==([family===:triangle ? 2 : 3]),kinds)
            @test all(g.model.mesh.get_elements(2,s)[1]==[2] for (_,s) in g.model.get_entities(2) if s ∉ pml)
            push!(inventories,sum(length(block) for s in pml for block in g.model.mesh.get_elements(2,s)[2]))
        end
    end
    @test inventories[1]==2inventories[2]
    @test inventories[2]==inventories[3]
end
