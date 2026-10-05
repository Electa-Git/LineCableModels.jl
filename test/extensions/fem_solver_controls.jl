@testitem "Gmsh FEM / native solver controls in managed and detached execution" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"solver-controls",
        terminal(:core,core(Material(kind=:conductor,rho=1.72e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[50.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    form = Formulation(:LineCableModelsFEM)
    solver = (mumps_ordering=0,petsc_prealloc=256)
    configured = computation_options(LineCableModelsFEM,ComputationOptions(;solver_threads=2,solver...))
    defaults = computation_options(LineCableModelsFEM,ComputationOptions())
    model = FEM._resolved_fem_model(problem,form)
    mktempdir() do root
        run = FEM._create_run(root)
        for options in (defaults,configured)
            dir = mktempdir(root)
            command = FEM._getdp_command("getdp","model.pro","mesh.msh",run,
                form,options,1,[1,2],dir)
            args = command.exec
            @test "-setnumber" in args
            @test "LinearSolver" in args
            @test ("MumpsOrdering" in args) 
            @test ("PetscPrealloc" in args) 
            if options===configured
                @test args[findfirst(==("MumpsOrdering"),args)+1] == "0"
                @test args[findfirst(==("PetscPrealloc"),args)+1] == "256"
                @test "OPENBLAS_NUM_THREADS=2" in command.env
            end
        end
        entry = export_data(:onelab,problem,form;file_name=joinpath(root,"bundle","study.pro"),
            mesh_options=(pml_layers=8,),solver_options=solver)
        data = read(joinpath(dirname(entry),"study_data.pro"),String)
        @test occursin("GetDPThreads = {1,",data)
        @test occursin("MumpsOrdering = {0,",data)
        @test occursin("PetscPrealloc = {256,",data)
        native = read(joinpath(dirname(entry),"formulations/helmholtz.pro"),String)*
            read(joinpath(dirname(entry),"formulations/solver.pro"),String)
        @test occursin("SetGlobalSolverOptions[FEMSolverOptions]",native)
        @test occursin("-mat_mumps_icntl_7",native)
        @test occursin("-petsc_prealloc",native)
        @test_throws ArgumentError export_data(:onelab,problem,form;
            file_name=joinpath(root,"bad","study.pro"),solver_options=(frequency_workers=2,))
        @test !isdir(joinpath(root,"bad"))
    end
end
