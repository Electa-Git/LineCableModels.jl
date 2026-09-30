@testitem "Gmsh FEM / physics selection, native maps, and resume isolation" tags=[
    :extension, :integration, :fem_numerical] begin
    using LineCableModels, Gmsh, LinearAlgebra, JSON3
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    wire = build(CableDesign, "physics-selection",
        terminal(:core,
            core(Material(kind = :conductor, rho = 2e-8); r = 0.005)))
    system = build(LineCableSystem, [wire, wire], [(0.0, 1.0), (0.4, 1.5)];
        connections = [Dict(:core=>1), Dict(:core=>2)], line_length = 1.0)
    problem = LineParametersProblem(system; frequencies = [1e4],
        earth_props = homogeneous(rho = 100.0, eps_r = 10.0, mu_r = 1.0))
    reductions = (
        reduce_bundle = false, kron_reduction = false, ideal_transposition = false)
    formulation = Formulation(:LineCableModelsFEM; options=reductions)
    controls = (gmsh_verbosity = 0, getdp_verbosity = 0,
        keep_run_directory = true, plot_field_maps = true)
    fw = compute(problem, formulation; options = (; controls..., trace = true))
    @test all(isfinite, Z(fw))
    @test all(isfinite, Y(fw))
    @test fw.details.data.fem.inputs.options.physics === Symbol("quasi-fw")
    @test details(fw).data.fem.run isa NamedTuple
    @test fw.details.data.fem.run.getdp_invocations == 1
    @test fw.details.data.fem.run.completed_columns == 2
    @test fw.details.data.fem.run.factorized_columns == 1
    @test length(fw.details.data.fem.run.map_paths) == 32
    # The axial map is total current, including displacement in air. Compare
    # native values with the independently exported electric field there.
    let maps=fw.details.data.fem.run.map_paths
        current_path=only(filter(p->endswith(p,"/jz_f0001_b0001.pos"),maps))
        electric_path=only(filter(p->endswith(p,"/ez_f0001_b0001.pos"),maps))
        found_air=false
        for (current_line,electric_line) in zip(eachline(current_path),eachline(electric_path))
            current=match(r"^ST\(([^)]*)\)\{([^}]*)\};",current_line)
            current===nothing && continue
            coordinates=parse.(Float64,split(current[1],','))
            minimum(coordinates[2:3:end])>2. || continue # Above both conductors.
            electric=match(r"^ST\(([^)]*)\)\{([^}]*)\};",electric_line)
            @test electric!==nothing
            @test current[1]==electric[1]
            j=parse.(Float64,split(current[2],','))
            e=parse.(Float64,split(electric[2],','))
            @test complex.(j[1:3],j[4:6]) ≈
                  (2π*1e4*im*8.8541878128e-12).*complex.(e[1:3],e[4:6]) rtol=1e-10
            found_air=true
            break
        end
        @test found_air
    end
    @test occursin("Coupled", fw.details.data.formulations.assumptions.admittance)
    path = fw.details.data.fem.run.run_directory
    @test isfile(joinpath(path, "raw/jobs/getdp-f0001-b0001-Pscalar.tsv"))
    @test !isfile(joinpath(path, "input/paths-f0001.pro"))
    checkpoint = JSON3.read(read(joinpath(path, "raw/jobs/getdp-f0001-b0001.json"), String))
    @test checkpoint.physics == "quasi-fw"
    @test length(checkpoint.checksums) == 20 # Z/P, timing, sixteen maps, scalar diagnostic
    repeated = compute(problem, formulation; options = (;
        controls..., resume_run_directory = path))
    @test Y(repeated) == Y(fw)
    @test repeated.details.data.fem.run.reused
    changed = Formulation(:LineCableModelsFEM; options=(; reductions..., Γ=.01im))
    @test_throws ArgumentError compute(problem, changed;
        options = (; controls..., resume_run_directory = path))

end

@testitem "Gmsh FEM / native complex reference extraction without maps" tags=[:extension,:integration,:fem_numerical] begin
    using LineCableModels, Gmsh, Printf
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire=build(CableDesign,"measurement",terminal(:core,
        core(Material(kind=:conductor,rho=2e-8);r=0.005)))
    system=build(LineCableSystem,[wire,wire],[(0.,1.),(0.4,-1.)];
        connections=[Dict(:core=>1),Dict(:core=>2)],line_length=1.)
    problem=LineParametersProblem(system;frequencies=[1e4],
        earth_props=homogeneous(rho=100.,eps_r=1.,mu_r=1.))
    execution=computation_options(LineCableModelsFEM,ComputationOptions(
        plot_field_maps=false,gmsh_verbosity=0,getdp_verbosity=3,solver_threads=1))
    lock(FEM.FEM_SESSION_LOCK) do
        session=FEM._start_gmsh(0)
        try
            for Γ in (0., .001+.002im)
                formulation=Formulation(:LineCableModelsFEM;options=(;physics=:quasi_fw,Γ,
                    reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
                model=FEM._resolved_fem_model(FEM._preflight_fem_problem(problem),formulation)
                mktempdir() do root
                    run=FEM._create_run(root)
                    geometry=FEM._build_geometry!(model,"native-measurement")
                    meshes=FEM._select_meshes!(run,model,geometry,execution,root)
                    FEM._prepare_run_inputs!(run,model)
                    FEM._write_json_atomic(joinpath(run.path,"input/computation.json"),
                        FEM._fem_input_record(model,formulation,execution))
                    assets=FEM._getdp_assets(joinpath(run.path,"input/getdp"))
                    # Manufacture complex postquantities in the actual owned
                    # extraction blocks. The assembled PDEs remain native.
                    source=assets.quasi_full
                    code=read(source,String)
                    code=replace(code,
                        "Re[{V} / UnitSource]"=>"Re[(0*{V}+Complex[10+\$FEMBasisTerminal,-8+\$FEMBasisTerminal]) / UnitSource]",
                        "Im[{V} / UnitSource]"=>"Im[(0*{V}+Complex[10+\$FEMBasisTerminal,-8+\$FEMBasisTerminal]) / UnitSource]",
                        "Re[{v}]"=>"Re[0*{v}+Complex[1+0.2*X[],2-0.1*X[]]]",
                        "Im[{v}]"=>"Im[0*{v}+Complex[1+0.2*X[],2-0.1*X[]]]")
                    code=replace(code,
                        "-Re[({U}+gamma2[]*{V}) / UnitSource]"=>"Re[(0*{U}+Complex[4,7])/UnitSource]",
                        "-Im[({U}+gamma2[]*{V}) / UnitSource]"=>"Im[(0*{U}+Complex[4,7])/UnitSource]",
                        "Re[CompY[{bt}]]"=>"Re[0*CompY[{bt}]+\$FEMBasisTerminal*Complex[-4,0.5]]",
                        "Im[CompY[{bt}]]"=>"Im[0*CompY[{bt}]+\$FEMBasisTerminal*Complex[-4,0.5]]")
                    write(source,code)
                    write(assets.model,replace(read(assets.model,String),
                        "UnitSource = 1.0;"=>"UnitSource = 2.5;"))
                    executable=FEM._getdp_selection(execution).path
                    checkdir=mkpath(joinpath(root,"parse"))
                    command=FEM._getdp_command(executable,assets.model,only(meshes),run,
                        formulation,execution,only(model.mesh_plans),[1,2],checkdir)
                    @test "PathDataPath" ∉ collect(command)
                    arguments=collect(command)
                    index=findfirst(==("-solve"),arguments)
                    splice!(arguments,index:index+1,["-check"])
                    physics_index=findfirst(==("Physics"),arguments)
                    splice!(arguments,physics_index-1:physics_index+1)
                    @test success(pipeline(Cmd(arguments),stdin=devnull,stdout=devnull,stderr=devnull))
                    FEM._run_getdp!(run,model,formulation,execution,meshes)
                    scan=FEM._parse_scan(run,model,formulation,execution)
                    @test isempty(scan.map_paths)
                    @test run.getdp_invocations==1
                    @test run.completed_columns==2
                    @test !any(startswith("paths"),readdir(joinpath(run.path,"input")))
                    for j in 1:2, i in 1:2
                        surface = i==1 ? 1+2im : 0im
                        plan = only(model.mesh_plans)
                        path_length = i==1 ? 1.0-0.005 : plan.domain_halfwidth+plan.pml_thickness[3]-1.0-0.005
                        line = j*(-4+0.5im)*path_length
                        expected=(complex(10+j,-8+j)-surface+2π*im*1e4*line)/2.5
                        @test scan.P[i,j,1] ≈ expected rtol=1e-11
                        @test scan.Z[i,j,1] ≈ (4+7im)/2.5+Γ^2*expected rtol=1e-11
                    end
                    for j in 1:2
                        timing=parse.(Float64,split(strip(read(joinpath(run.path,"raw/jobs",
                            @sprintf("getdp-f0001-b%04d-timing.tsv",j)),String)),'\t'))
                        @test last(timing)==(j==1)
                    end
                end
            end
        finally
            FEM._finish_gmsh(session)
        end
    end
end
