@testitem "Gmsh FEM / detached native export and caller ownership" tags=[:extension] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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
                    file_name=entry,options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),)) == entry
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
            @test occursin("Physics = {1, Choices{1=\"Helmholtz\"}",data)
            @test !isfile(joinpath(dirname(entry),"formulations","quasi-tem.pro"))
            @test !occursin(pkgdir(LineCableModels),data)
            geometry = read(joinpath(dirname(entry),"geometry","physical.geo"),String)
            @test occursin("FEMReceiverX",geometry)
            @test !occursin("FEMReceiverColumn",geometry)
            @test_throws ArgumentError export_data(:onelab,problem,formulation;file_name=entry)
            @test_throws ArgumentError export_data(:onelab,problem,formulation;
                file_name=joinpath(root,"invalid","study.pro"),options=(frequency_workers=2,))
            @test !isdir(joinpath(root,"invalid"))
            second = export_data(:onelab,system,formulation;earth_props=earth,
                frequencies=problem.frequencies,file_name=joinpath(root,"system","study.pro"),
                options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),))
            @test read(joinpath(dirname(second),"study_data.pro"),String)==data
            write(joinpath(dirname(entry),"user-notes.txt"),"keep me")
            write(joinpath(dirname(entry),"formulations/quasi-tem.pro"),"owned old asset")
            open(joinpath(dirname(entry),".onelab-export-files"),"a") do io
                println(io,"formulations/quasi-tem.pro")
            end
            export_data(:onelab,problem,formulation;file_name=entry,
                options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),),overwrite=true)
            @test read(joinpath(dirname(entry),"user-notes.txt"),String)=="keep me"
            @test !isfile(joinpath(dirname(entry),"formulations/quasi-tem.pro"))
            marker = joinpath(dirname(entry),".onelab-export-files")
            recorded = read(marker,String)
            write(marker,replace(recorded,"README.md\n"=>""))
            before = read(entry,String)
            @test_throws ArgumentError export_data(:onelab,problem,formulation;
                file_name=entry,options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),),overwrite=true)
            @test read(entry,String)==before
            @test !any(startswith(".onelab-export-"),readdir(root))
            @test gmsh.model.get_current()=="caller-owned"
            @test gmsh.model.list()==models && gmsh.view.get_tags()==views
            @test gmsh.option.get_number("Mesh.MeshSizeMax")==0.123
            write(marker,recorded)
            write(marker,recorded*"../outside.txt\n")
            @test_throws ArgumentError export_data(:onelab,problem,formulation;
                file_name=entry,options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),),overwrite=true)
            write(marker,recorded)
            gmsh.open(replace(entry,r"\.pro$"=>".geo"))
            options = computation_options(LineCableModelsFEM,ComputationOptions(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,)))
            model = FEM._resolved_fem_model(problem,formulation)
            available = Set((d,t,gmsh.model.get_physical_name(d,t)) for (d,t) in gmsh.model.get_physical_groups())
            @test all(group in available for group in FEM._expected_physical_groups(model))
            @test length(gmsh.model.get_entities(2)) >= 12
            @test !isempty(gmsh.model.mesh.field.list())
            gmsh.model.mesh.generate(2)
            mesh = joinpath(root,"reopened.msh")
            gmsh.write(mesh)
            open(mesh) do io
                @test readline(io)=="\$MeshFormat"
                @test split(readline(io))[2]=="1"
            end
            ratios=FEM._inspect_loaded_mesh(model,mesh)
            @test length(ratios)==length(model.terminal_ids)
            # Allow twice the 0.25 target because Delaunay edges can be smaller than their targets.
            @test all(r -> isfinite(r) && 0<r<=.5,ratios)
            @test only(gmsh.parser.get_number("MeasurementLineSizeFactor")) == .25
        finally
            gmsh.finalize()
        end
    end
end

@testitem "Gmsh FEM / detached native frequency scan" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh, JSON3
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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
            file_name=joinpath(root,"scan with spaces","study.pro"),options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),))
        bundle = dirname(entry)
        geo = joinpath(bundle,"study.geo")
        data = joinpath(bundle,"study_data.pro")
        session = FEM._start_gmsh(0)
        try
            gmsh.open(geo)
            @test gmsh.onelab.get_number("Inputs/00Run frequency scan") == [0.]
            for (name,label) in (("Inputs/Air/01Sigma","Sigma [S/m]"),
                    ("Inputs/Air/02Epsilon","Epsilon [F/m]"),("Inputs/Air/03Mu","Mu [H/m]"),
                    ("Inputs/Cases/0001/02Gamma real","Gamma real [1/m]"),
                    ("Inputs/Cases/0001/04Soil sigma","Soil sigma [S/m]"),
                    ("Boundary/Derived/Medium 0/01Root real","Root real [1/m]"))
                @test JSON3.read(gmsh.onelab.get(name)).label==label
            end
            @test !any(name->occursin(r"\[[^]]*/[^]]*\]",name),gmsh.onelab.get_names())
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
        physics, code = "helmholtz", 1
        log = read(`$command -setnumber RunFrequencyScan 1 -setnumber Physics $code`,String)
        @test !occursin("Error",log)
        # Check the ordered mesh/save/solve events, not just existence of
        # result files that could have been left by an earlier case.
        events = collect(eachmatch(r"Writing '[^\n]*study\.msh'|Print -> '[^\n]*completed\.txt'",log))
        # Automatic check writes an initial mesh before the two compute steps.
        @test length(events) == 5
        if length(events) == 5
            @test startswith(events[1].match,"Writing")
            @test startswith(events[2].match,"Writing")
            @test occursin("f0001-",events[3].match)
            @test startswith(events[4].match,"Writing")
            @test occursin("f0002-",events[5].match)
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
            # Gmsh can reorder a mesh when the metadata path length changes; compare numbers at 1e-9 scaled.
            for (actual,expected) in zip(read.(tables,String),scan_tables)
                rows=split.(split(chomp(actual),'\n'),'\t')
                reference=split.(split(chomp(expected),'\n'),'\t')
                @test rows[1:2]==reference[1:2]
                @test [r[1:4] for r in rows[3:end]]==[r[1:4] for r in reference[3:end]]
                values=[complex(parse(Float64,r[5]),parse(Float64,r[6])) for r in rows[3:end]]
                target=[complex(parse(Float64,r[5]),parse(Float64,r[6])) for r in reference[3:end]]
                scale=maximum(abs(target[k]) for (k,r) in enumerate(reference[3:end]) if r[1]==r[2])
                @test maximum(abs.(values.-target))<=1e-9*scale
            end
        end
        # An unchecked native Run visits only the selected frequency.
        log = read(`$command -setnumber RunFrequencyScan 0 -setnumber FrequencyIndex 2`,String)
        @test count(r"Print -> '[^\n]*completed\.txt'",log) == 1
        @test occursin("f0002-helmholtz-b0000/completed.txt",log)
        @test !occursin("f0001-helmholtz-b0000/completed.txt",log)
        log = read(`$command -setnumber RunFrequencyScan 1 -setnumber RunAction 0`,String)
        @test count(r"Writing '[^\n]*study\.msh'",log) == 3 # Check plus two mesh-only compute steps.
        @test !occursin("Print -> '",log)
        @test all(!isfile(joinpath(run_dir(i,"helmholtz"),"completed.txt")) for i in 1:2)
        singleton = LineParametersProblem(system;frequencies=[50.],earth_props=problem.earth_props)
        single = export_data(:onelab,singleton,formulation;
            file_name=joinpath(root,"single","study.pro"),options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),))
        log = read(`$mesher $single -option $options -setnumber PlotFieldMaps 0 -setnumber RunFrequencyScan 1 -run`,String)
        @test count(r"Print -> '[^\n]*completed\.txt'",log) == 1
        @test isfile(joinpath(dirname(single),"results","f0001-helmholtz-b0000","completed.txt"))
    end
end

@testitem "Gmsh FEM / detached mesh controls and failed publication" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
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
            file_name=joinpath(root,"study.pro"),options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),))
        mesh = joinpath(root,"study.msh")
        session = FEM._start_gmsh(0)
        function inventory()
            gmsh.clear(); gmsh.open(mesh)
            types,tags,nodes = gmsh.model.mesh.get_elements(2)
            # Keep a coordinate inventory to check the native mesh controls.
            volume = [[gmsh.model.mesh.get_node(n)[1] for n in block] for block in nodes]
            lines = sum(length(first(gmsh.model.mesh.get_elements(1,c)[2]))
                for (d,p) in gmsh.model.get_physical_groups(1)
                if startswith(gmsh.model.get_physical_name(d,p),"LCM/measurement_line/")
                for c in gmsh.model.get_entities_for_physical_group(d,p))
            (;volume,triangles=sum(length,tags),lines)
        end
        try
            run(`$mesher $(joinpath(root,"study.geo")) -setnumber BuildMesh 1 -0 -v 2`)
            coarse = inventory()
            cp(mesh,joinpath(root,"coarse.msh"))
            @test !occursin("VoltageRefinements",read(joinpath(root,"study_data.pro"),String))
        finally
            FEM._finish_gmsh(session)
        end
        # An invalid source must invalidate a prior completion before failing.
        result = joinpath(root,"results","f0001-helmholtz-b0000")
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
        partial = joinpath(root,"results","f0001-helmholtz-b0001")
        @test isfile(joinpath(partial,"completed.txt"))
        @test length(readlines(joinpath(partial,"matrices","P-primitive.tsv"))) == 4
        @test !isfile(joinpath(partial,"matrices","Y.tsv"))
        @test read(datafile,String) == original
    end
end

@testitem "Gmsh FEM / detached native numerical parity and relocation" tags=[:extension,:integration,:fem_numerical] setup=[TemporaryFEMRuntime] begin
    using Gmsh, LinearAlgebra, SHA
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    cd(fem_test_runtime_directory)
    try
    FEM = Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    copper = Material(kind=:conductor,rho=1.72e-8)
    wire = build(CableDesign,"export-parity",terminal(:core,core(copper;r=0.005)))
    system = build(LineCableSystem,[wire,wire],[(0.,0.1),(0.2,-0.1)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem = LineParametersProblem(system;frequencies=[50.,10000.],
        earth_props=homogeneous(rho=100.,eps_r=10.))
    options = (overrides=(PmlSideLayers=8,PmlTopLayers=6,PmlBottomLayers=4,PmlSideGrading=3.,PmlTopGrading=2.,PmlBottomGrading=1.,MumpsOrdering=0,PetscPrealloc=256),frequency_workers=1,solver_threads=1,
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
            options=(overrides=(PmlSideLayers=8,PmlTopLayers=6,PmlBottomLayers=4,PmlSideGrading=3.,PmlTopGrading=2.,PmlBottomGrading=1.,MumpsOrdering=0,PetscPrealloc=256),))
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

        physics, code = :helmholtz, 1
        selected = Formulation(:LineCableModelsFEM;options=(
            physics,reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        reference = compute(problem,selected;options)
        record = details(reference).data.fem
        for index in eachindex(problem.frequencies)
            mesh = joinpath(record.run.run_directory,"mesh","frequency_$(lpad(index,4,'0')).msh")
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
    finally
        cd(fem_test_working_directory)
        rm(fem_test_runtime_directory;recursive=true)
    end

end

@testitem "Gmsh FEM / checked meshes publish the selected frequency" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh, SHA
    N=NativeFEMFixtures;g=Gmsh.gmsh
    base=N.problem(;frequency=.1,rho=.1,eps_r=1.,radius=.0425,positions=[(0.,1.),(1.,1.)])
    problem=LineParametersProblem(base.system;temperature=20.,frequencies=[.1,1e6],earth_props=base.earth_props)
    N.bundle(problem) do directory,entry
        g.initialize(String[],false,false)
        try
            g.option.set_number("General.Terminal",0)
            g.onelab.set_string("Gmsh/Action",["check"])
            hashes=String[];counts=Int[];widths=Float64[]
            for index in (1,2,1)
                g.onelab.set_number("Inputs/01Frequency case",[Float64(index)])
                g.open(replace(entry,r"\.pro$"=>".geo"))
                @test g.option.get_number("Solver.AutoCheck")==1
                @test g.onelab.get_string("Mesh/Current mesh/00Status")==["Generated"]
                @test g.onelab.get_number("Mesh/Current mesh/01Case index")==[Float64(index)]
                @test g.onelab.get_number("Mesh/Current mesh/02Frequency [Hz]")==[problem.frequencies[index]]
                push!(widths,only(g.parser.get_number("DomainHalfwidth")))
                push!(counts,sum(length,g.model.mesh.get_elements(2)[2]))
                push!(hashes,bytes2hex(open(sha256,replace(entry,r"\.pro$"=>".msh"))))
            end
            @test all(>(0),counts)
            @test widths[1]==widths[3]>widths[2]
            @test hashes[1]!=hashes[2]
            # A failed save must not publish the newly meshed case.
            geometry=replace(entry,r"\.pro$"=>".geo")
            script=read(geometry,String)
            write(geometry,replace(script,"Save StrCat(CurrentDirectory, \"model.msh\")"=>
                "Save StrCat(CurrentDirectory, \"missing-directory/model.msh\")"))
            @test_throws Exception g.open(geometry)
            @test g.onelab.get_string("Mesh/Current mesh/00Status")==["No mesh"]
            @test g.onelab.get_number("Mesh/Current mesh/01Case index")==[0.]
            write(geometry,script)
            # Test parsing failure in a fresh native session after the save error.
            Gmsh.finalize()
            g.initialize(String[],false,false)
            g.option.set_number("General.Terminal",0)
            g.onelab.set_string("Gmsh/Action",["check"])
            g.open(geometry)
            @test g.onelab.get_string("Mesh/Current mesh/00Status")==["Generated"]
            # A failed reparse must invalidate the preceding mesh label.
            write(joinpath(directory,"geometry/physical.geo"),"Include \"missing-geometry.geo\";")
            @test_throws Exception g.open(replace(entry,r"\.pro$"=>".geo"))
            @test g.onelab.get_string("Mesh/Current mesh/00Status")==["No mesh"]
            @test g.onelab.get_number("Mesh/Current mesh/01Case index")==[0.]
            @test g.onelab.get_number("Mesh/Current mesh/02Frequency [Hz]")==[0.]
        finally
            Gmsh.finalize()
        end
    end
end
