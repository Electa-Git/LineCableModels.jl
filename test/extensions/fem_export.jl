@testitem "Gmsh FEM / detached native export and caller ownership" tags=[:extension] begin
    using Gmsh
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"export-test",
        terminal(:core,core(Material(kind=:conductor,rho=1.72e-8);r=0.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,0.1),(0.2,-0.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    earth = homogeneous(rho=100.,eps_r=10.)
    problem = LineParametersProblem(system;frequencies=[50.,10000.],earth_props=earth)
    formulation = Formulation(:LineCableModelsFEM;options=(
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    mktempdir() do root
        entry = joinpath(root,"bundle with spaces","study.pro")
        gmsh.initialize(String[],false,false)
        try
            gmsh.model.add("caller-owned")
            gmsh.view.add("caller-view")
            gmsh.option.set_number("Mesh.MeshSizeMax",0.123)
            gmsh.option.set_number("Mesh.Binary",0)
            models,views = gmsh.model.list(),gmsh.view.get_tags()
            withenv("LINECABLEMODELS_GETDP"=>"/not/a/solver","DISPLAY"=>"") do
                @test export_data(:onelab,problem,formulation;
                    file_name=entry,mesh_options=(pml_layers=8,)) == entry
            end
            @test gmsh.model.get_current() == "caller-owned"
            @test gmsh.model.list() == models
            @test gmsh.view.get_tags() == views
            @test gmsh.option.get_number("Mesh.MeshSizeMax") == 0.123
            @test gmsh.option.get_number("Mesh.Binary") == 0
            @test !isdir(joinpath(dirname(entry),"work"))
            @test !any(endswith(".msh"),readdir(dirname(entry)))
            data = read(joinpath(dirname(entry),"study_data.pro"),String)
            @test occursin("Material_1_conductor_Sigma()",data)
            @test occursin("Connection_2 = 2;",data)
            @test occursin("UnitSource = 1.;",data)
            @test occursin("Physics = {1, Choices{1=\"quasi-fw\"}",data)
            @test !isfile(joinpath(dirname(entry),"formulations","quasi-tem.pro"))
            @test !occursin(pkgdir(LineCableModels),data)
            geometry = read(joinpath(dirname(entry),"geometry","case-0001.geo"),String)
            @test occursin("Mesh.Binary = 1;",geometry)
            # Native Gmsh serialization previously truncated these PML ratios
            # and omitted the alternating triangle arrangement entirely.
            ratios = [parse(Float64,m[1]) for m in eachmatch(r"Using Progression ([0-9eE.+-]+);",geometry)]
            @test any(==(exp(((192/191)*log(1536))/8)),ratios)
            @test count("AlternateLeft;",geometry)==12
            @test occursin("In Surface",geometry)
            @test_throws ArgumentError export_data(:onelab,problem,formulation;file_name=entry)
            @test_throws ArgumentError export_data(:onelab,problem,formulation;
                file_name=joinpath(root,"invalid","study.pro"),mesh_options=(frequency_workers=2,))
            @test !isdir(joinpath(root,"invalid"))
            second = export_data(:onelab,system,formulation;earth_props=earth,
                frequencies=problem.frequencies,file_name=joinpath(root,"system","study.pro"),
                mesh_options=(pml_layers=8,))
            @test read(joinpath(dirname(second),"study_data.pro"),String)==data
            write(joinpath(dirname(entry),"user-notes.txt"),"keep me")
            write(joinpath(dirname(entry),"formulations/quasi-tem.pro"),"owned old asset")
            open(joinpath(dirname(entry),".onelab-export-files"),"a") do io
                println(io,"formulations/quasi-tem.pro")
            end
            export_data(:onelab,problem,formulation;file_name=entry,
                mesh_options=(pml_layers=8,),overwrite=true)
            @test read(joinpath(dirname(entry),"user-notes.txt"),String)=="keep me"
            @test !isfile(joinpath(dirname(entry),"formulations/quasi-tem.pro"))
            marker = joinpath(dirname(entry),".onelab-export-files")
            recorded = read(marker,String)
            write(marker,replace(recorded,"README.md\n"=>""))
            before = read(entry,String)
            @test_throws ArgumentError export_data(:onelab,problem,formulation;
                file_name=entry,mesh_options=(pml_layers=8,),overwrite=true)
            @test read(entry,String)==before
            @test !any(startswith(".onelab-export-"),readdir(root))
            @test gmsh.model.get_current()=="caller-owned"
            @test gmsh.model.list()==models && gmsh.view.get_tags()==views
            @test gmsh.option.get_number("Mesh.MeshSizeMax")==0.123
            write(marker,recorded)
            write(marker,recorded*"../outside.txt\n")
            @test_throws ArgumentError export_data(:onelab,problem,formulation;
                file_name=entry,mesh_options=(pml_layers=8,),overwrite=true)
            write(marker,recorded)
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            options = computation_options(LineCableModelsFEM,ComputationOptions(pml_layers=8))
            model = FEM._resolved_fem_model(problem,formulation,options)
            available = Set((d,t,gmsh.model.get_physical_name(d,t)) for (d,t) in gmsh.model.get_physical_groups())
            @test all(group in available for group in FEM._expected_physical_groups(model))
            @test length(gmsh.model.get_entities(2)) == 16
            @test !isempty(gmsh.model.mesh.field.list())
            gmsh.model.mesh.generate(2)
            mesh = joinpath(root,"reopened.msh")
            gmsh.write(mesh)
            open(mesh) do io
                @test readline(io)=="\$MeshFormat"
                @test split(readline(io))[2]=="1"
            end
            @test FEM._inspect_loaded_mesh(model,mesh) === nothing
        finally
            gmsh.finalize()
        end
    end
end

@testitem "Gmsh FEM / editable physical and exterior mesh factors" tags=[:extension] begin
    using Gmsh
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"editable-mesh",terminal(:core,
        core(Material(kind=:conductor,rho=1.72e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[.1,1e4],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    formulation = LineCableModelsFEM()
    controls = (domain_skin_depths=8.,pml_layers=4,mesh_size_factor=3.,exterior_mesh_size_factor=8.)
    mktempdir() do root
        entry = export_data(:onelab,problem,formulation;file_name=joinpath(root,"study.pro"),mesh_options=controls)
        data = read(joinpath(root,"study_data.pro"),String)
        @test occursin("1=\"0.1 Hz\"",data)
        @test occursin("Frequencies() = "*FEM._pro_array(problem.frequencies),data)
        @test !occursin("0.10000000000000001 Hz",data)
        session = FEM._start_gmsh(0)
        function inventory(model)
            pml = gmsh.model.get_entities_for_physical_group(2,model.tags.pml)
            contours = [gmsh.model.get_entities_for_physical_group(1,model.tags.terminal_contour_base+i) for i in 1:2]
            triangles = sum(length(block) for s in pml for block in gmsh.model.mesh.get_elements(2,s)[2])
            lines = Dict(c=>sum(length,gmsh.model.mesh.get_elements(1,c)[2]) for (_,c) in gmsh.model.get_entities(1))
            metals = [sum(lines[c] for c in curves) for curves in contours]
            return (;triangles,lines,metals,nodes=length(first(gmsh.model.mesh.get_nodes())))
        end
        try
            observed = []
            # Include both directions of editing, and a change to the physical
            # factor. Match each detached mesh to a fresh Julia construction.
            for (factor,exterior) in ((3.,8.),(3.,1.),(1.5,4.),(3.,8.))
                options = computation_options(LineCableModelsFEM,ComputationOptions(;
                    controls...,mesh_size_factor=factor,exterior_mesh_size_factor=exterior))
                model = FEM._resolved_fem_model(problem,formulation,options)
                plan = model.mesh_plans[2]
                gmsh.clear()
                geometry = FEM._build_geometry!(model,"managed-controls",plan)
                FEM._configure_mesh!(model,geometry,plan)
                gmsh.model.mesh.generate(2)
                expected = inventory(model)
                gmsh.clear(); gmsh.parser.clear()
                gmsh.onelab.set_number("Inputs/01Frequency case",[2.])
                gmsh.onelab.set_number("Mesh/01Physical size factor",[factor])
                gmsh.onelab.set_number("Mesh/02Exterior size factor",[exterior])
                gmsh.open(replace(entry,r"\.pro$"=>".geo"))
                @test only(gmsh.parser.get_number("MeshBulk")) ≈ plan.domain_mesh_size
                @test only(gmsh.parser.get_number("MeshRemoteAir")) ≈ plan.exterior_mesh_sizes[1]
                @test only(gmsh.parser.get_number("MeshRemoteEarth")) ≈ plan.exterior_mesh_sizes[2]
                @test only(gmsh.parser.get_number("MeshInterface")) ≈ plan.interface_mesh_size
                @test only(gmsh.parser.get_number("ConductorRegion1First")) ≈
                    FEM._conductor_mesh_sizes(model,first(model.region_plans),plan).first_size
                gmsh.model.mesh.generate(2)
                @test FEM._inspect_loaded_mesh(model,"editable-native-mesh") === nothing
                actual = inventory(model)
                @test actual.lines == expected.lines
                @test actual.triangles == expected.triangles
                @test actual.metals == expected.metals
                push!(observed,actual)
            end
            @test observed[1].triangles < observed[2].triangles
            @test observed[1].nodes < observed[2].nodes
            @test observed[1].metals == observed[2].metals
            # Gmsh can choose a different unstructured triangulation on reload;
            # the prescribed edges and structured PML must return identically.
            @test observed[1].lines == observed[4].lines
            @test observed[1].triangles == observed[4].triangles
            @test read(joinpath(root,"study_data.pro"),String) == data
        finally
            gmsh.onelab.clear(); gmsh.parser.clear()
            FEM._finish_gmsh(session)
        end
    end
end

@testitem "Gmsh FEM / detached native frequency scan" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh, JSON3
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    wire = build(CableDesign,"scan-test",terminal(:core,
        core(Material(kind=:conductor,rho=1.72e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[50.,10000.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    formulation = Formulation(:LineCableModelsFEM;options=(
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    mesher = Gmsh.gmsh_jll.gmsh()
    mktempdir() do root
        entry = export_data(:onelab,problem,formulation;
            file_name=joinpath(root,"scan with spaces","study.pro"),mesh_options=(pml_layers=8,))
        bundle = dirname(entry)
        geo = joinpath(bundle,"study.geo")
        data = joinpath(bundle,"study_data.pro")
        session = FEM._start_gmsh(0)
        try
            gmsh.open(geo)
            @test gmsh.onelab.get_number("Inputs/00Run frequency scan") == [0.]
            gmsh.onelab.set_number("Inputs/01Frequency case",[2.])
            # Reparse in one server session, as ONELAB does after a UI edit.
            # Hidden counters let the declaration control Loop even after it
            # was enabled; visible parameters preserve user loop attributes.
            for scan in (1.,0.,1.)
                gmsh.onelab.set_number("Inputs/00Run frequency scan",[scan])
                gmsh.parser.clear()
                gmsh.parser.parse(data)
                counter = JSON3.read(gmsh.onelab.get("Inputs/00Scan frequency index"))
                manual = JSON3.read(gmsh.onelab.get("Inputs/01Frequency case"))
                @test counter.attributes.Loop == (scan == 1 ? "1" : "")
                @test counter.choices == [1.,2.]
                @test !counter.visible
                @test manual.visible == (scan == 0)
                @test manual.values == [2.]
            end
        finally
            gmsh.onelab.clear()
            FEM._finish_gmsh(session)
        end
        # Exercise Gmsh's real ONELAB loop, including native solver discovery
        # and mesh-before-solve ordering, with no project runtime interpreter.
        options = joinpath(root,"options.opt")
        open(options,"w") do io
            println(io,"Solver.Name0 = \"GetDP\";")
            println(io,"Solver.Executable0 = ",FEM._pro_string(getdp),";")
            # Gmsh prefixes SocketName with its home directory, even when the
            # supplied name starts with a slash. Keep concurrent runs isolated.
            socket = relpath(joinpath(root,"onelab.sock"),homedir())
            println(io,"Solver.SocketName = ",FEM._pro_string(socket),";")
            println(io,"Solver.AutoLoadDatabase = 0; Solver.AutoSaveDatabase = 0;")
        end
        command = `$mesher $entry -option $options -setnumber PlotFieldMaps 0 -run`
        run_dir(i,physics) = joinpath(bundle,"results","f"*lpad(i,4,'0')*"-"*physics*"-b0000")
        physics, code = "quasi-fw", 1
        log = read(`$command -setnumber RunFrequencyScan 1 -setnumber Physics $code`,String)
        @test !occursin("Error",log)
        # Check the ordered mesh/save/solve events, not just existence of
        # result files that could have been left by an earlier case.
        events = collect(eachmatch(r"Writing '[^\n]*study\.msh'|Print -> '[^\n]*completed\.txt'",log))
        @test length(events) == 4
        if length(events) == 4
            @test startswith(events[1].match,"Writing")
            @test occursin("f0001-",events[2].match)
            @test startswith(events[3].match,"Writing")
            @test occursin("f0002-",events[4].match)
        end
        for (i,f) in enumerate(problem.frequencies)
            result = run_dir(i,physics)
            @test parse.(Float64,split(read(joinpath(result,"completed.txt"),String))) == [i,f,code,0]
            tables = [joinpath(result,"matrices",name*".tsv") for name in ("Z-primitive","P-primitive","Z","P","Y")]
            scan_tables = read.(tables,String)
            # Independent manual mesh/solve of each frequency must produce
            # the same numerical output as the corresponding scan step.
            run(`$mesher $geo -setnumber BuildMesh 1 -setnumber FrequencyIndex $i -0 -v 2`)
            run(`$getdp $entry -msh $(joinpath(bundle,"study.msh")) -solve LineCableModelsFEM -setnumber FrequencyIndex $i -setnumber Physics $code -setnumber PlotFieldMaps 0 -v 2`)
            @test read.(tables,String) == scan_tables
        end
        # An unchecked native Run visits only the selected frequency.
        log = read(`$command -setnumber RunFrequencyScan 0 -setnumber FrequencyIndex 2`,String)
        @test count("Print -> '",log) == 1
        @test occursin("f0002-quasi-fw-b0000/completed.txt",log)
        @test !occursin("f0001-quasi-fw-b0000/completed.txt",log)
        log = read(`$command -setnumber RunFrequencyScan 1 -setnumber RunAction 0`,String)
        @test count(r"Writing '[^\n]*study\.msh'",log) == 2
        @test !occursin("Print -> '",log)
        @test all(!isfile(joinpath(run_dir(i,"quasi-fw"),"completed.txt")) for i in 1:2)
        singleton = LineParametersProblem(system;frequencies=[50.],earth_props=problem.earth_props)
        single = export_data(:onelab,singleton,formulation;
            file_name=joinpath(root,"single","study.pro"),mesh_options=(pml_layers=8,))
        log = read(`$mesher $single -option $options -setnumber PlotFieldMaps 0 -setnumber RunFrequencyScan 1 -run`,String)
        @test count("Print -> '",log) == 1
        @test isfile(joinpath(dirname(single),"results","f0001-quasi-fw-b0000","completed.txt"))
    end
end

@testitem "Gmsh FEM / detached mesh controls and failed publication" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    gmsh = Gmsh.gmsh
    wire = build(CableDesign,"native-controls",terminal(:core,
        core(Material(kind=:conductor,rho=1.72e-8);r=.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,.1),(.2,-.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[50.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    formulation = Formulation(:LineCableModelsFEM;options=(
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    executable = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
    mesher = Gmsh.gmsh_jll.gmsh()
    mktempdir() do root
        entry = export_data(:onelab,problem,formulation;
            file_name=joinpath(root,"study.pro"),mesh_options=(pml_layers=8,))
        mesh = joinpath(root,"study.msh")
        session = FEM._start_gmsh(0)
        function inventory()
            gmsh.clear(); gmsh.open(mesh)
            types,tags,nodes = gmsh.model.mesh.get_elements(2)
            # Keep a coordinate inventory to check the native mesh controls.
            volume = [[gmsh.model.mesh.get_node(n)[1] for n in block] for block in nodes]
            lines = sum(length(first(gmsh.model.mesh.get_elements(1,c)[2]))
                for (d,p) in gmsh.model.get_physical_groups(1)
                if startswith(gmsh.model.get_physical_name(d,p),"LCM/voltage_path/")
                for c in gmsh.model.get_entities_for_physical_group(d,p))
            (;volume,triangles=sum(length,tags),lines)
        end
        try
            run(`$mesher $(joinpath(root,"study.geo")) -setnumber BuildMesh 1 -0 -v 2`)
            coarse = inventory()
            cp(mesh,joinpath(root,"coarse.msh"))
            @test !occursin("VoltageRefinements",read(joinpath(root,"study_data.pro"),String))
            run(`$mesher $(joinpath(root,"study.geo")) -setnumber BuildMesh 1 -setnumber MeshRefinements 1 -0 -v 2`)
            uniform = inventory()
            @test uniform.triangles == 4coarse.triangles
            @test uniform.lines == 2coarse.lines
        finally
            FEM._finish_gmsh(session)
        end
        # An invalid source must invalidate a prior completion before failing.
        result = joinpath(root,"results","f0001-quasi-fw-b0000")
        mkpath(result); marker = joinpath(result,"completed.txt")
        write(marker,"previous success\n")
        datafile = joinpath(root,"study_data.pro")
        original = read(datafile,String)
        write(datafile,replace(original,"UnitSource = 1.;"=>"UnitSource = 0.;"))
        command = `$executable $entry -msh $mesh -solve LineCableModelsFEM -v 2`
        @test !success(pipeline(command;stdout=devnull,stderr=devnull))
        @test !isfile(marker)
        @test read(datafile,String) == replace(original,"UnitSource = 1.;"=>"UnitSource = 0.;")
        write(datafile,original)
        stale = joinpath(result,"maps","az_f0001_b0001.pos")
        mkpath(dirname(stale)); write(stale,"old field")
        solve = `$executable $entry -msh $(joinpath(root,"coarse.msh")) -solve LineCableModelsFEM -setnumber PlotFieldMaps 0 -v 2`
        run(solve)
        @test isfile(marker)
        @test !isfile(stale)
        @test length(readlines(joinpath(result,"matrices","Y.tsv"))) == 6
        run(solve)
        @test length(readlines(joinpath(result,"matrices","Y.tsv"))) == 6
        @test length(readlines(marker)) == 1
        run(`$solve -setnumber BasisTerminal 1`)
        partial = joinpath(root,"results","f0001-quasi-fw-b0001")
        @test isfile(joinpath(partial,"completed.txt"))
        @test length(readlines(joinpath(partial,"matrices","P-primitive.tsv"))) == 4
        @test !isfile(joinpath(partial,"matrices","Y.tsv"))
        @test read(datafile,String) == original
    end
end

@testitem "Gmsh FEM / detached native numerical parity and relocation" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh, LinearAlgebra, SHA
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    copper = Material(kind=:conductor,rho=1.72e-8)
    wire = build(CableDesign,"export-parity",terminal(:core,core(copper;r=0.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,0.1),(0.2,-0.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[50.,10000.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    options = (pml_layers=(8,6,4),pml_grading=(3.,2.,1.),frequency_workers=1,solver_threads=1,
        mumps_ordering=0,petsc_prealloc=256,
        mesh_policy=:remesh,keep_run_directory=true,trace=true,
        gmsh_verbosity=0,getdp_verbosity=0)
    function read_matrix(path)
        rows = [split(line,'\t') for line in readlines(path)[3:end]]
        n = maximum(parse(Int,row[1]) for row in rows)
        m = maximum(parse(Int,row[2]) for row in rows)
        matrix = zeros(ComplexF64,n,m)
        for row in rows
            matrix[parse(Int,row[1]),parse(Int,row[2])] = complex(parse(Float64,row[5]),parse(Float64,row[6]))
        end
        matrix
    end
    function compare_components(actual,expected)
        @test size(actual)==size(expected)
        for part in (real,imag)
            a,b = part.(actual),part.(expected)
            scale = maximum(abs,b)
            @test all(abs.(a.-b) .<= 2e-9.*abs.(b) .+ 100eps(Float64)*scale)
        end
    end
    mktempdir() do root
        formulation = Formulation(:LineCableModelsFEM;options=(
            reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        export_data(:onelab,problem,formulation;file_name=joinpath(root,"original","study.pro"),
            mesh_options=(pml_layers=(8,6,4),pml_grading=(3.,2.,1.)),
            solver_options=(mumps_ordering=0,petsc_prealloc=256))
        bundle = joinpath(root,"relocated bundle with spaces")
        mv(joinpath(root,"original"),bundle)
        entry = joinpath(bundle,"study.pro")
        files = readlines(joinpath(bundle,".onelab-export-files"))
        before = Dict(file=>bytes2hex(open(sha256,joinpath(bundle,file))) for file in files)
        getdp = FEM._getdp_selection(computation_options(LineCableModelsFEM,ComputationOptions())).path
        @test all(!endswith(file,".py") for file in files)
        @test !isfile(joinpath(bundle,"requirements.txt"))
        @test success(Cmd(`$getdp $entry -v 0`;dir=root))
        @test !success(pipeline(Cmd(`$getdp $entry -setnumber Physics 0 -v 0`;dir=root),stdout=devnull,stderr=devnull))

        physics, code = :quasi_fw, 1
        selected = Formulation(:LineCableModelsFEM;options=(
            physics,reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        reference = compute(problem,selected;options)
        record = details(reference).data.fem
        for index in eachindex(problem.frequencies)
            mesh = joinpath(record.run.run_directory,"mesh",index==length(problem.frequencies) ? "model.msh" : "frequency_0001.msh")
            command = `$getdp $entry -solve LineCableModelsFEM -msh $mesh -setnumber Physics $code -setnumber FrequencyIndex $index -v 2`
            @test success(addenv(Cmd(command;dir=root),"PATH"=>mktempdir(root)))
            run = joinpath(bundle,"results","f"*lpad(index,4,'0')*"-"*replace(String(physics),'_'=>'-')*"-b0000")
            @test isfile(joinpath(run,"completed.txt"))

            compare_components(read_matrix(joinpath(run,"matrices","Z-primitive.tsv")),record.primitive.Z_primitive[:,:,index])
            compare_components(read_matrix(joinpath(run,"matrices","P-primitive.tsv")),record.primitive.P_primitive[:,:,index])
            compare_components(read_matrix(joinpath(run,"matrices","Z.tsv")),Z(reference)[:,:,index])
            compare_components(read_matrix(joinpath(run,"matrices","Y.tsv")),Y(reference)[:,:,index])
            @test !isfile(joinpath(run,"paths.pro"))
            @test !isfile(joinpath(bundle,"measurements.py"))

        end
        @test all(bytes2hex(open(sha256,joinpath(bundle,file)))==hash for (file,hash) in before)
    end
end

@testitem "Gmsh FEM / detached polygon and disconnected terminal geometry" tags=[:extension] begin
    using Gmsh
    gmsh = Gmsh.gmsh
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    copper = Material(kind=:conductor,rho=1.72e-8)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    polygons = (
        Polygon(((-.003,-.001),(0.,-.001),(0.,.001),(-.003,.001))),
        Polygon(((0.,-.001),(.003,-.001),(.003,.001),(0.,.001))),
        Polygon(((.004,.002),(.006,.002),(.006,.004),(.004,.004))))
    design = build(CableDesign,"export-polygon",Enclosure(:matrix,
        terminal(:core,assembly((solid(copper,shape) for shape in polygons)...));
        primitive=Disk(.01),fill=dielectric))
    system = build(LineCableSystem,design,Pose2(.06,-.2,.43);connections=Dict(:core=>1))
    problem = LineParametersProblem(system;frequencies=[50.],earth_props=homogeneous(rho=100.))
    form = LineCableModelsFEM()
    model = FEM._resolved_fem_model(problem,form,
        computation_options(LineCableModelsFEM,ComputationOptions(pml_layers=8)))
    mktempdir() do root
        entry = export_data(:onelab,problem,form;file_name=joinpath(root,"bundle","study.pro"),mesh_options=(pml_layers=8,))
        gmsh.initialize(String[],false,false)
        try
            gmsh.option.set_number("General.Terminal",0)
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            @test length(gmsh.model.get_entities_for_physical_group(2,3001))==3
            curves = gmsh.model.get_entities_for_physical_group(1,model.tags.voltage_path_base+1)
            ref = only(gmsh.model.get_entities_for_physical_group(0,model.tags.voltage_reference_base+1))
            from = gmsh.model.get_value(0,ref,Float64[])
            endpoints = [gmsh.model.get_value(0,p,Float64[]) for (_,p) in
                gmsh.model.get_boundary([(1,c) for c in curves],true,false,false)]
            to = only(filter(!=(from),endpoints))
            # Rotated lower-left vertex of the first polygon. The dielectric
            # enclosure and disconnected third polygon do not change it.
            expected = [.06-.003cos(.43)+.001sin(.43),
                -.2-.003sin(.43)-.001cos(.43),0.]
            @test to ≈ expected atol=1e-14 rtol=0
            @test from[1] ≈ expected[1] atol=1e-14 rtol=0
            plan = only(model.mesh_plans)
            @test from[2] ≈ -plan.domain_halfwidth-plan.pml_thickness[3]
            @test any(!isempty(gmsh.model.mesh.get_embedded(2,s)) for (_,s) in gmsh.model.get_entities(2))
            gmsh.model.mesh.generate(2)
            mesh = joinpath(root,"polygon.msh")
            gmsh.write(mesh)
            @test FEM._inspect_loaded_mesh(model,mesh) === nothing
            gmsh.model.remove_physical_groups([(1,model.tags.voltage_path_base+1)])
            @test_throws LineCableModelsFEMError FEM._inspect_loaded_mesh(model,mesh)
        finally
            gmsh.finalize()
        end
    end
end
