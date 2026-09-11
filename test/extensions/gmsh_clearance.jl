@testitem "Gmsh FEM / exterior clearance / touching trefoil mesh" tags=[:extension] begin
    using Gmsh
    using Measurements: measurement
    using LineCableModels
    ext = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    copper = Material(kind = :conductor, rho = 1.7e-8)
    formulation = Formulation(:LineCableModelsFEM;
        options = (ideal_transposition = false,))
    execution = computation_options(LineCableModelsFEM, (mesh_policy=:remesh, gmsh_verbosity=0,))
    for radius in (0.01, measurement(0.01, 1e-4))
        design = build(CableDesign, "touching-mesh",
            Group(:core, Region(:metal, Disk(radius), copper)))
        system = @test_logs (:warn, r"Cable placements adjusted") build(LineCableSystem,
            trefoil(design; center = at(0, -0.1), spacing = 2radius,
                connections = (core = (1, 2, 3),)))
        problem = LineParametersProblem(system;
            earth_props = homogeneous(rho = 100.0), frequencies = [50.0])
        normalized = ext._preflight_fem_problem(problem)
        @test normalized.system.clearances ≈ nominal.(system.clearances)
        @test getproperty.(normalized.system.positions, :x) ≈ nominal.(getproperty.(system.positions, :x))
        model = ext._resolved_fem_model(normalized, formulation)
        session = ext._start_gmsh(0)
        try
            mktempdir() do directory
                run = ext._create_run(directory)
                geometry = ext._build_geometry!(model, "touching-trefoil")
                mesh = ext._select_mesh!(run, model, geometry, execution, directory)
                @test isfile(mesh)
                @test run.mesh_source === :generated
                @test !isempty(first(Gmsh.gmsh.model.mesh.get_nodes()))
                ext._validate_mesh_file(model, mesh)
            end
        finally
            ext._finish_gmsh(session)
        end
    end
end
