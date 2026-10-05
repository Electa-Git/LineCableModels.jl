@testitem "FEM / editable interface footprint preserves terminal contours" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    p=N.problem(;frequency=1e6,rho=.1,positions=[(0.,1.),(2.,1.)])
    records=[]
    for factor in (1.,12.)
        N.geometry(p;options=(interface_refinement_factor=factor,),mesh=true) do g,log
            @test only(g.parser.get_number("InterfaceRefinementFactor"))==factor
            contours=[sum(length(block) for c in g.model.get_entities_for_physical_group(1,4000+i) for block in g.model.mesh.get_elements(1,c)[2]) for i in 1:2]
            push!(records,(nodes=length(first(g.model.mesh.get_nodes())),contours=contours))
            @test !any(l->occursin("closer than the geometrical tolerance",l),log)
        end
    end
    @test records[2].nodes>=records[1].nodes
    @test records[1].contours==records[2].contours
end
