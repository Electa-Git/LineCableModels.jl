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
    @test formulation_options(LineCableModelsFEM, (;)).physics === Symbol("quasi-tem")
    for value in (:quasi_tem, "quasi-tem", Symbol("quasi-tem"))
        @test formulation_options(LineCableModelsFEM, (;physics=value)).physics === Symbol("quasi-tem")
    end
    for value in (:quasi_fw, "quasi-fw", Symbol("quasi-fw"))
        @test formulation_options(LineCableModelsFEM, (;physics=value)).physics === Symbol("quasi-fw")
    end
    @test_throws ArgumentError formulation_options(LineCableModelsFEM, (;physics=:fullwave))
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
    @test tem.details.fem.inputs.options.physics !== fw.details.fem.inputs.options.physics
    @test all(isfinite, Y(fw))
    @test fw.details.fem.inputs.options.physics === Symbol("quasi-fw")
    @test fw.details.fem.run.run_directory != tem.details.fem.run.run_directory
    @test details(fw).fem.run isa NamedTuple
    @test fw.details.fem.run.getdp_invocations == 1
    @test fw.details.fem.run.completed_columns == 2
    @test fw.details.fem.timing.factorized_columns == 1
    @test length(fw.details.fem.run.map_paths) == 24
    @test length(tem.details.fem.run.map_paths) == 18
    @test occursin("Coupled", fw.details.formulations.assumptions.admittance)
    path = fw.details.fem.run.run_directory
    @test isfile(joinpath(path,"raw/jobs/getdp-f0001-b0001-Pscalar.tsv"))
    @test isfile(joinpath(path,"input/paths-f0001.pro"))
    checkpoint = JSON3.read(read(joinpath(path,"raw/jobs/getdp-f0001-b0001.json"),String))
    @test checkpoint.physics == "quasi-fw"
    @test length(checkpoint.checksums) == 16 # Z/P, timing, twelve maps, scalar diagnostic
    repeated = compute(problem, last(formulations); options=(;controls...,resume_run_directory=path))
    @test Y(repeated) == Y(fw)
    @test repeated.details.fem.timing.reused
    @test_throws ArgumentError compute(problem, first(formulations);
        options=(;controls...,resume_run_directory=path))
    # Both physics use the same named ONELAB selector and resolution.
    session = FEM._start_gmsh(0)
    try
        model = FEM._resolved_fem_model(FEM._preflight_fem_problem(problem), last(formulations))
        run = FEM.FEMRun(path,FEM.completed,"test",:none,"")
        FEM._publish_transport!(run,model,joinpath(path,"input/model_data.pro"),
            joinpath(path,"mesh/model.msh"),last(formulations), computation_options(LineCableModelsFEM, controls))
        @test gmsh.onelab.get_number("LineCableModels/FEM/physics") == [1.0]
    finally
        FEM._finish_gmsh(session)
    end
end

@testitem "Gmsh FEM / coupled voltage paths / retained edge field and path convergence" tags=[:extension,:integration,:fem_numerical] setup=[TestFixtures] begin
    using Gmsh,LinearAlgebra,Printf,TOML
    include(joinpath(pkgdir(LineCableModels),"test/support/edge_field_control.jl"))
    include(joinpath(pkgdir(LineCableModels),"test/support/fem_mesh_control.jl"))
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    BLAS.set_num_threads(1)
    metal=Material(kind=:conductor,rho=2e-8)
    wires=[build(CableDesign,"path-electrode-$i",terminal(:core,core(metal;r)))
        for (i,r) in enumerate((.005,.007))]
    systems=LineCableSystem[build(LineCableSystem,wires,positions;
        connections=[Dict(:core=>1),Dict(:core=>2)])
        for positions in ([(0.,1.),(.4,1.5)],[(0.,1.),(0.,-1.)])]
    push!(systems,build(LineCableSystem,TestFixtures.coaxial_design(),Pose2(0.,-1.);
        connections=Dict(:core=>1,:sheath=>2)))
    formulation=Formulation(:LineCableModelsFEM;options=(physics=:quasi_fw,
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    controls=(gmsh_verbosity=0,getdp_verbosity=0,keep_run_directory=true,trace=true,
        plot_field_maps=true,frequency_workers=1,solver_threads=1)
    directory=mktempdir(;prefix="lcm-fem-edge-provisional-",cleanup=false)
    println("FEM edge-field evidence: ",directory);flush(stdout)
    notes=Dict{String,Any}[]
    @testset "layout $layout" for (layout,system) in enumerate(systems)
        problem=LineParametersProblem(system;frequencies=[10000.],
            earth_props=homogeneous(rho=100.,eps_r=10.))
        current=compute(problem,formulation;options=controls)
        run_directory=current.details.fem.run.run_directory
        model=FEM._resolved_fem_model(FEM._preflight_fem_problem(problem),formulation)
        plan=only(model.mesh_plans)
        info=FEM.FEMRun(run_directory,FEM.completed,"path control",:none,"")
        mesh=joinpath(run_directory,"mesh/model.msh")
        source=joinpath(run_directory,"input/getdp/model.pro")
        executable=current.details.fem.inputs.getdp_provenance.path
        fields=Vector{Any}();endpoints=Tuple{Float64,Float64}[]
        session=FEM._start_gmsh(0)
        try
            gmsh.open(mesh)
            tags,coordinates,_=gmsh.model.mesh.get_nodes()
            nodes=Dict(tag=>(coordinates[3i-2],coordinates[3i-1]) for (i,tag) in enumerate(tags))
            for terminal in 1:2
                curves=gmsh.model.get_entities_for_physical_group(1,model.tags.terminal_contour_base+terminal)
                boundary=unique(reduce(vcat,[first(gmsh.model.mesh.get_nodes(1,c,true)) for c in curves]))
                # Independently establish the current physical endpoint rule.
                ordered=sort([nodes[tag] for tag in boundary];by=p->(p[2],p[1]))
                push!(endpoints,first(ordered))
            end
            for basis in 1:2
                map_path=only(filter(path->endswith(path,@sprintf("bt_mesh_f0001_b%04d.pos",basis)),
                    current.details.fem.run.map_paths))
                push!(fields,EdgeFieldControl.read_field(map_path))
            end
        finally
            FEM._finish_gmsh(session)
        end
        values=Matrix{ComplexF64}[];independent=Matrix{ComplexF64}[]
        roundoff_bounds=Matrix{Float64}[]
        for segments in (64,128,256,512)
            work=joinpath(directory,"layout$layout-segments$segments");mkpath(work)
            path=joinpath(work,"paths.pro")
            session=FEM._start_gmsh(0)
            try
                FEM._write_voltage_paths(path,mesh,plan,model;shell_segments=segments)
            finally
                FEM._finish_gmsh(session)
            end
            command=FEM._getdp_command(executable,source,mesh,info,formulation,
                computation_options(LineCableModelsFEM,controls),plan,[1,2],work;reuse_factorization=false)
            open(joinpath(work,"getdp.log"),"w") do io
                run(pipeline(`$command -setstring PathDataPath $path`;stdout=io,stderr=io))
            end
            raw=zeros(ComplexF64,2,2);scalar=similar(raw)
            for (quantity,matrix) in (("P",raw),("Pscalar",scalar)),basis in 1:2
                file=joinpath(work,"raw/jobs",@sprintf("getdp-f0001-b%04d-%s.tsv",basis,quantity))
                @test FEM._valid_job_raw(file,2,1,10000.,basis)
                for line in eachline(file)
                    row=split(line);matrix[parse(Int,row[3]),basis]=complex(parse(Float64,row[5]),parse(Float64,row[6]))
                end
            end
            expected=similar(raw);bounds=zeros(2,2)
            for terminal in 1:2,basis in 1:2
                vertices=EdgeFieldControl.ray(endpoints[terminal],model.centre,
                    plan.domain_radius,plan.shell_outer_radius,segments)
                line=EdgeFieldControl.integrate(fields[basis],vertices)
                expected[terminal,basis]=scalar[terminal,basis]+2pi*10000im*line.value
                bounds[terminal,basis]=2pi*10000line.bound
                for component in (real,imag)
                    budget=.01abs(component(expected[terminal,basis]))
                    uncertainty=2pi*10000line.bound
                    @test uncertainty<=budget/4
                    @test abs(component(raw[terminal,basis]-expected[terminal,basis]))+uncertainty<=budget
                end
            end
            push!(values,raw);push!(independent,expected);push!(roundoff_bounds,bounds)
        end
        resolved=true
        for terminal in 1:2,basis in 1:2,component in (real,imag)
            control=component(last(independent)[terminal,basis])
            convergence=if endpoints[terminal][1]==model.centre[1]
                # On the centreline, the pullback ray is exactly straight at
                # every segmentation. There is no geometric path truncation;
                # the independently accumulated field arithmetic bound applies.
                u=maximum(v[terminal,basis] for v in roundoff_bounds)
                (resolved=u<=.01abs(control)/4,uncertainty=u)
            else
                EdgeFieldControl.richardson([component(v[terminal,basis]) for v in values],.01abs(control))
            end
            resolved &= convergence.resolved
            @test abs(component(current.details.fem.primitive.P_primitive[terminal,basis,1])-control)+convergence.uncertainty<=.01abs(control)
        end
        push!(notes,Dict("layout"=>layout,"run"=>run_directory,"segments"=>[64,128,256,512],
            "path_convergence_resolved"=>resolved,"status"=>"provisional extraction record; see path/mesh/domain checks"))
        open(joinpath(directory,"outcomes.toml"),"w") do io
            TOML.print(io,Dict("layouts"=>notes))
        end
        @test resolved
        study_directory=joinpath(directory,"mesh-domain-layout$layout");mkpath(study_directory)
        study=FEMMeshControl.study(FEM,problem,formulation,controls;directory=study_directory,
            observable=(z,p)->p)
        for entry in eachindex(last(study.mesh_values)),component in (real,imag)
            scale=abs(component(last(independent)[entry]))
            mesh_convergence=EdgeFieldControl.richardson([component(v[entry]) for v in study.mesh_values],.01scale)
            @test mesh_convergence.resolved
            differences=abs.(diff([component(v[entry]) for v in study.domain_values]))
            @test length(differences)==2 && last(differences)<=first(differences)/2
            @test 2last(differences)<=.01scale/4
        end
    end
end

@testitem "Manual quasi-full / native gauge and current-basis invariance" tags=[:extension, :integration, :fem_numerical] begin
    using LineCableModels, Gmsh, LinearAlgebra, Printf
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
    root=mktempdir(;prefix="lcm-fem-path-provisional-",cleanup=false)
    println("FEM path evidence: ",root);flush(stdout)
    begin
        info = FEM._create_run(root)
        path = joinpath(info.path, "paths.pro")
        mesh, renumbered = lock(FEM.FEM_SESSION_LOCK) do
            session = FEM._start_gmsh(0)
            try
                geometry = FEM._build_geometry!(model, "quasi-full-test")
                mesh = only(FEM._select_meshes!(info, model, geometry, computation_options(LineCableModelsFEM, fem_controls), root))
                FEM._prepare_run_inputs!(info, model)
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
                # Same elements and coordinates; a different node ordering
                # gives GetDP a different admissible spanning tree.
                gmsh.open(mesh)
                tags, _, _ = gmsh.model.mesh.get_nodes()
                gmsh.model.mesh.renumber_nodes(tags, reverse(tags))
                renumbered = joinpath(info.path, "mesh", "renumbered.msh")
                gmsh.write(renumbered)
                # GetDP enumerates edges from the elements. Reverse the MSH 4.1
                # element blocks and their numbering as well as the nodes.
                lines = readlines(renumbered)
                first = findfirst(==("\$Elements"), lines)
                last = findfirst(==("\$EndElements"), lines)
                blocks = Tuple{String,Vector{String}}[]
                cursor = first + 2
                while cursor < last
                    count = parse(Int, split(lines[cursor])[4])
                    push!(blocks, (lines[cursor], lines[cursor+1:cursor+count]))
                    cursor += count + 1
                end
                open(renumbered, "w") do io
                    foreach(line -> println(io,line), lines[1:first+1])
                    element = 0
                    for (header, rows) in reverse(blocks)
                        println(io, header)
                        for row in reverse(rows)
                            element += 1
                            println(io, element, " ", split(row;limit=2)[2])
                        end
                    end
                    foreach(line -> println(io,line), lines[last:end])
                end
                (mesh, renumbered)
            finally
                FEM._finish_gmsh(session)
            end
        end
        executable = FEM._getdp_selection(computation_options(LineCableModelsFEM, fem_controls)).path
        source = joinpath(pkgdir(LineCableModels), "ext", "LineCableModelsGmshExt", "getdp", "model.pro")
        function solve_case(name, mesh; reuse=true)
            directory = joinpath(root, name); mkpath(directory)
            command = FEM._getdp_command(executable, source, mesh, info, fem, computation_options(LineCableModelsFEM, fem_controls),
                only(model.mesh_plans), [1,2], directory; reuse_factorization=reuse)
            open(joinpath(directory, "getdp.log"), "w") do io
                run(pipeline(`$command -setstring PathDataPath $path`; stdout=io, stderr=io))
            end
            matrices = map(("Z", "P", "Pscalar")) do quantity
                matrix = zeros(ComplexF64, 2, 2)
                for basis in 1:2
                    raw = joinpath(directory, "raw/jobs", @sprintf("getdp-f0001-b%04d-%s.tsv", basis, quantity))
                    @test FEM._valid_job_raw(raw, 2, 1, 1e4, basis)
                    for line in eachline(raw)
                        row = split(line)
                        matrix[parse(Int,row[3]),basis] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
                    end
                end
                matrix
            end
            matrices
        end
        z, p, scalar = solve_case("original", mesh)
        z2, p2, scalar2 = solve_case("renumbered", renumbered)
        # Gauge-dependent scalar diagnostics have no required change.
        @test z2 ≈ z rtol=1e-7
        @test p2 ≈ p rtol=1e-7
        z3, p3, _ = solve_case("fresh-factors", mesh; reuse=false)
        @test z3 ≈ z rtol=1e-8
        @test p3 ≈ p rtol=1e-8

    end
end
