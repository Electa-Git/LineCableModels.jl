@testitem "Gmsh FEM / production edge paths / affine field and Stokes controls" tags=[:extension] begin
    using Gmsh
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    triangles=[((0.,0.),(2.,0.),(2.,1.)),((0.,0.),(2.,1.),(0.,1.))]
    a=(3.0,-2.0); b=.7
    integrate(vertices;kwargs...)=sum((a[1]-b*y)*dx+(a[2]+b*x)*dy
        for (x,y,dx,dy) in FEM._voltage_path_points(triangles,vertices;kwargs...))
    p=(.3,0.0);q=(1.7,1.0)
    exact=a[1]*(q[1]-p[1])+a[2]*(q[2]-p[2])+b*(p[1]*q[2]-p[2]*q[1])
    @test integrate([p,q]) ≈ exact rtol=1e-10
    @test integrate([q,p]) ≈ -exact rtol=1e-10
    @test integrate([(0.,0.),(2.,1.)]) ≈ 2a[1]+a[2]
    @test integrate([(0.,0.),(2.,0.),(2.,1.),(0.,1.),(0.,0.)]) ≈ 2b*2
    @test integrate([(.5,0.),(.5,1.)];active=[false,true]) ≈ (a[2]+.5b)*.75
    @test_throws ErrorException integrate([(0.,0.),(0.,2.)])
    @test_throws DimensionMismatch integrate([p,q];active=[true])
end

@testitem "Gmsh FEM / physics selection, native maps, and resume isolation" tags=[:extension, :integration, :fem_numerical] begin
    using LineCableModels, Gmsh, LinearAlgebra, JSON3
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    @test formulation_options(LineCableModelsFEM, FormulationOptions()).data.physics === Symbol("quasi-tem")
    for value in (:quasi_tem, "quasi-tem", Symbol("quasi-tem"))
        @test formulation_options(LineCableModelsFEM, FormulationOptions(physics=value)).data.physics === Symbol("quasi-tem")
    end
    for value in (:quasi_fw, "quasi-fw", Symbol("quasi-fw"))
        @test formulation_options(LineCableModelsFEM, FormulationOptions(physics=value)).data.physics === Symbol("quasi-fw")
    end
    @test_throws ArgumentError formulation_options(LineCableModelsFEM, FormulationOptions(physics=:fullwave))
    wire = build(CableDesign, "physics-selection", terminal(:core,
        core(Material(kind=:conductor, rho=2e-8); r=0.005)))
    system = build(LineCableSystem, [wire,wire], [(0.,1.), (.4,1.5)];
        connections=[Dict(:core=>1),Dict(:core=>2)], line_length=1.)
    problem = LineParametersProblem(system; frequencies=[1e4],
        earth_props=homogeneous(rho=100.0,eps_r=10.,mu_r=1.))
    reductions = (reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    formulations = [Formulation(:LineCableModelsFEM; options=(;reductions...,physics))
        for physics in (:quasi_tem,:quasi_fw)]
    controls = (gmsh_verbosity=0,getdp_verbosity=0,keep_run_directory=true,plot_field_maps=true)
    results = compute(problem, formulations; options=(;controls...,trace=true))
    tem, fw = results
    @test all(isfinite,Z(fw)) # No cross-formulation numerical authority.
    @test tem.details.data.fem.inputs.options.physics !== fw.details.data.fem.inputs.options.physics
    @test all(isfinite, Y(fw))
    @test fw.details.data.fem.inputs.options.physics === Symbol("quasi-fw")
    @test fw.details.data.fem.run.run_directory != tem.details.data.fem.run.run_directory
    @test details(fw).data.fem.run isa NamedTuple
    @test fw.details.data.fem.run.getdp_invocations == 1
    @test fw.details.data.fem.run.completed_columns == 2
    @test fw.details.data.fem.timing.factorized_columns == 1
    @test length(fw.details.data.fem.run.map_paths) == 24
    @test length(tem.details.data.fem.run.map_paths) == 18
    @test occursin("Coupled", fw.details.data.formulations.assumptions.admittance)
    path = fw.details.data.fem.run.run_directory
    @test isfile(joinpath(path,"raw/jobs/getdp-f0001-b0001-Pscalar.tsv"))
    @test isfile(joinpath(path,"input/paths-f0001.pro"))
    checkpoint = JSON3.read(read(joinpath(path,"raw/jobs/getdp-f0001-b0001.json"),String))
    @test checkpoint.physics == "quasi-fw"
    @test length(checkpoint.checksums) == 16 # Z/P, timing, twelve maps, scalar diagnostic
    repeated = compute(problem, last(formulations); options=(;controls...,resume_run_directory=path))
    @test Y(repeated) == Y(fw)
    @test repeated.details.data.fem.timing.reused
    @test_throws ArgumentError compute(problem, first(formulations);
        options=(;controls...,resume_run_directory=path))
    # Both physics use the same named ONELAB selector and resolution.
    session = FEM._start_gmsh(0)
    try
        model = FEM._resolved_fem_model(FEM._preflight_fem_problem(problem), last(formulations))
        run = FEM.FEMRun(path,FEM.completed,"test",:none,"")
        FEM._publish_transport!(run,model,joinpath(path,"input/model_data.pro"),
            joinpath(path,"mesh/model.msh"),last(formulations), computation_options(LineCableModelsFEM, ComputationOptions(controls)))
        @test gmsh.onelab.get_number("LineCableModels/FEM/physics") == [1.0]
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / voltage path preparation preserves the caller session" tags=[:extension] begin
    using LineCableModels, Gmsh
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    wire = build(CableDesign, "bare", terminal(:core,
        core(Material(kind=:conductor, rho=2e-8); r=0.005)))
    positions = [(0.,1.), (.4,1.5)]
    system = build(LineCableSystem, [wire,wire], positions;
        connections=[Dict(:core=>1),Dict(:core=>2)], line_length=1.)
    problem = LineParametersProblem(system; frequencies=[1e4],
        earth_props=homogeneous(rho=100.0,eps_r=10.,mu_r=1.))
    fem = Formulation(:LineCableModelsFEM;
        options=(physics=:quasi_fw,reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    fem_controls = (gmsh_verbosity=0,getdp_verbosity=3,mesh_policy=:remesh)
    model = FEM._resolved_fem_model(FEM._preflight_fem_problem(problem), fem)
    mktempdir() do root
        info = FEM._create_run(root)
        path = joinpath(info.path, "paths.pro")
        lock(FEM.FEM_SESSION_LOCK) do
            session = FEM._start_gmsh(0)
            try
                geometry = FEM._build_geometry!(model, "quasi-full-test")
                mesh = only(FEM._select_meshes!(info, model, geometry, computation_options(LineCableModelsFEM, ComputationOptions(fem_controls)), root))
                # Exercise a caller model whose name collides with the mesh
                # basename; path preparation must preserve both its identity
                # and its contents rather than opening another "model".
                gmsh.model.add(splitext(basename(mesh))[1])
                caller_model = gmsh.model.get_current()
                gmsh.model.add_discrete_entity(0, 991)
                initial_models = sort(gmsh.model.list())
                FEM._write_voltage_paths(path,mesh,only(model.mesh_plans),model)
                @test gmsh.model.get_current() == caller_model
                @test sort(gmsh.model.list()) == initial_models
                @test gmsh.model.get_entities(0) == [(0,991)]
                gmsh.model.remove()
            finally
                FEM._finish_gmsh(session)
            end
        end
    end
end
