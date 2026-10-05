@testitem "FEM / touching trefoil clearances survive native mesh export" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh, Measurements
    N=NativeFEMFixtures
    copper=Material(kind=:conductor,rho=1.7e-8)
    for r in (.01,measurement(.01,1e-4))
        design=build(CableDesign,"touching",Group(:core,Region(:metal,Disk(r),copper)))
        system=@test_logs (:warn,r"Cable placements adjusted") build(LineCableSystem,
            trefoil(design;center=at(0,-.1),spacing=2r,connections=(core=(1,2,3),)))
        p=LineParametersProblem(system;earth_props=homogeneous(rho=100.),frequencies=[1e6])
        normalized=N.FEM._preflight_fem_problem(p)
        @test normalized.system.clearances≈nominal.(system.clearances)
        N.geometry(normalized;mesh=true) do g,log
            model=N.FEM._resolved_fem_model(normalized,N.FEM.LineCableModelsFEM())
            ratios=N.FEM._inspect_loaded_mesh(model,"touching trefoil")
            @test length(ratios)==3 && all(isfinite,ratios)
            @test !isempty(first(g.model.mesh.get_nodes()))
        end
    end
end
