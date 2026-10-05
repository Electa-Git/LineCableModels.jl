@testitem "Gmsh FEM / merged arc endpoints do not create a full circle" tags=[:extension] begin
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


@testitem "FEM / native conductor contour and skin controls" tags=[:extension] setup=[NativeFEMFixtures] begin
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
