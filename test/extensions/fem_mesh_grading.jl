@testitem "FEM / directional PML progression and interval nodes" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    N.geometry(N.problem(;frequency=1e6,rho=100.);
            options=(overrides=(PmlSideLayers=16,PmlTopLayers=24,PmlBottomLayers=48,PmlSideGrading=1.,PmlTopGrading=2.,PmlBottomGrading=3.,PmlSideThicknessFactor=1.,PmlTopThicknessFactor=2.,PmlBottomThicknessFactor=3.),),mesh=true) do g,log
        for (name,grading) in (("Side",1.),("Top",2.),("Bottom",3.))
            @test only(g.parser.get_number("Pml$(name)Grading"))==grading
            @test only(g.parser.get_number("FEMPml$(name)Layers"))>=Dict("Side"=>16,"Top"=>24,"Bottom"=>48)[name]
        end
        D=only(g.parser.get_number("DomainHalfwidth"));L=only(g.parser.get_number("PmlBottomThickness"))
        count=Int(only(g.parser.get_number("FEMPmlBottomLayers")))
        edges=filter(g.model.get_entities(1)) do e
            b=g.model.get_bounding_box(e...)
            abs(b[2]+D+L)<1e-6 && abs(b[5]+D)<1e-6 && abs(b[1]-b[4])<1e-6
        end
        @test !isempty(edges)
        for edge in edges
            coords=reshape(g.model.mesh.get_nodes(1,edge[2],true)[2],3,:)
            actual=sort((-D .- coords[2,:])./L)
            expected=expm1.(3 .* (0:count)./count)./expm1(3.)
            @test length(actual)==count+1
            @test maximum(abs,actual-expected)<2e-7
        end
    end
end
