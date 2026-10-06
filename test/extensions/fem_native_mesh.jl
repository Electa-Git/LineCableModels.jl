@testitem "FEM / earth sizing length, ceiling and domain option" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh
    N=NativeFEMFixtures
    expressions=["length"=>"FEMEarthBaseLength", "ceiling"=>"FEMEarthSizingCeilingActive",
        "size"=>"FEMEarthSizingLength", "domain"=>"DomainHalfwidth", "bulk"=>"MeshBulk"]
    for (rho,f,eps_r) in ((100.,50.,1.),(1000.,1e8,12.),(Inf,50.,1.),(Inf,1e8,1.))
        p=N.problem(;rho,frequency=f,eps_r)
        N.parameters(p;expressions) do a,b
            @test isequal(a,b)
            w=2pi*f;mu=4pi*1e-7;eps=8.8541878128e-12*eps_r
            q=sqrt(complex(-w*w*mu*eps,w*mu/rho))
            expected=min(real(q)>0 ? inv(real(q)) : Inf,2pi/abs(q))
            cap=sqrt(2e5/(w*mu))
            @test a["length"]≈expected rtol=2e-14
            @test a["ceiling"]==Float64(expected>cap)
            @test a["size"]≈min(expected,cap) rtol=2e-14
            @test a["domain"]≈max(5.,2min(expected,cap)) rtol=2e-14
        end
    end
    E=N.FEM
    @test E.computation_options(E.LineCableModelsFEM,ComputationOptions(domain_size_factor=3.)).data.domain_size_factor==3.
    err=try E.computation_options(E.LineCableModelsFEM,ComputationOptions(domain_skin_depths=2.));catch e;e;end
    @test err isa ArgumentError
    @test sprint(showerror,err)=="ArgumentError: unknown LineCableModelsFEM computation options: (:domain_skin_depths,)"
    N.bundle(N.problem()) do directory,entry
        @test occursin("DomainSizeFactor",read(joinpath(directory,"model_data.pro"),String))
        @test !occursin("DomainSkinDepths",read(joinpath(directory,"model_data.pro"),String))
    end
end

@testitem "FEM / long buried measurement lines retain remote grading" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    N.geometry(N.problem(;frequency=.1,rho=100.,eps_r=1.,positions=[(0.,-.1)]);
            mesh=true) do g,log
        @test only(g.parser.get_number("DomainHalfwidth"))>1e4
        lines=g.model.get_entities_for_physical_group(1,
            Int(only(g.parser.get_number("MEASUREMENT_LINE"))))
        @test !isempty(lines)
        counts=map(lines) do curve
            _,tags,_=g.model.mesh.get_elements(1,curve)
            sum(length,tags)
        end
        # The local conductor target must not seed the whole distant line.
        @test sum(counts)<10_000
        _,_,blocks=g.model.mesh.get_elements(2)
        nodes=Set(vcat(blocks...))
        for curve in lines
            tags,_,_=g.model.mesh.get_nodes(1,curve,true)
            @test isempty(intersect(nodes,Set(tags)))
        end
    end
end

@testitem "FEM / loss-aware PML intervals and native parser parity" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    expressions=["a_air"=>"FEMRootA~{0}","b_air"=>"FEMRootB~{0}",
        "a_earth"=>"FEMRootA~{1}","b_earth"=>"FEMRootB~{1}","target"=>"FEMPmlTarget",
        "ppw_air"=>"FEMPmlPPW~{0}","ppw_earth"=>"FEMPmlPPW~{1}"]
    for (label,d) in (("side","Side"),("top","Top"),("bottom","Bottom"))
        append!(expressions,[label*"_eta"=>"Pml$(d)Eta",label*"_A"=>"Pml$(d)Strength",
            label*"_L"=>"Pml$(d)Thickness",label*"_N"=>"FEMPml$(d)Layers"])
    end
    for f in (50.,1e8,3e8)
        N.parameters(N.problem(;frequency=f);expressions) do a,b
            @test isequal(a,b)
            for m in ("air","earth")
                q=complex(a["a_"*m],a["b_"*m]);rate=real(q)
                factor=rate>0 ? clamp(sqrt(.1abs(q)/rate),1.,3.) : 3.
                @test a["ppw_"*m]≈10factor rtol=2e-14
            end
            for (d,media) in (("side",("air","earth")),("top",("air",)),("bottom",("earth",)))
                X=a[d*"_L"]*(1+a[d*"_A"]/4-im*a[d*"_eta"]*a[d*"_A"]/4)
                phases=map(media) do m
                    q=complex(a["a_"*m],a["b_"*m]);E=real(q*X)
                    a["ppw_"*m]*abs(q)*abs(X)*(E>0 ? min(1,a["target"]/E) : 1)
                end
                @test a[d*"_N"]==max(16,ceil(maximum(phases)/(2pi)))
            end
        end
    end
end

@testitem "FEM / earth interface layer stays clear of buried cable" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    ex=["thickness"=>"FEMEarthLayerThickness","active"=>"FEMEarthLayerActive",
        "clipped"=>"FEMEarthLayerClippedOrOmitted","decay"=>"MeshDecayEarth",
        "wave"=>"MeshWaveEarth","remote"=>"MeshRemoteEarth","cable"=>"CableSize~{0}"]
    for y in (1.,-1.,-.1)
        N.parameters(N.problem(;frequency=1e8,rho=1000.,eps_r=12.,positions=[(0.,y)]);expressions=ex) do a,b
            @test isequal(a,b)
            candidate=y<0 ? min(a["decay"],.5*(abs(y)-.01-a["cable"])) : a["decay"]
            active=a["wave"]<a["remote"] && candidate>=2a["wave"]
            @test a["active"]==active
            @test a["thickness"]≈(active ? candidate : 0.)
            @test a["clipped"]==((a["wave"]<a["remote"]) && a["thickness"]<a["decay"])
            y<0 && @test a["thickness"]<=.5*(abs(y)-.01-a["cable"])
        end
    end
end

@testitem "FEM / box edge grading and independent measurement mesh" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh
    N=NativeFEMFixtures
    problem=N.problem(;frequency=1e6,rho=.1,eps_r=1.,radius=.01)
    N.geometry(problem;options=(mesh_size_factor=2.5,),mesh=true) do g,log
        D=only(g.parser.get_number("DomainHalfwidth"));x=only(g.parser.get_number("Xcenter"))+D
        wave=only(g.parser.get_number("MeshWaveEarth"));bulk=only(g.parser.get_number("MeshBulk"))
        last=only(g.parser.get_number("MeshRemoteEarth"));first=min(wave,bulk,last)
        curves=filter(g.model.get_entities(1)) do entity
            box=g.model.get_bounding_box(entity...)
            abs(box[1]-x)<1e-6 && abs(box[4]-x)<1e-6 && abs(box[2]+D)<1e-6 && abs(box[5])<1e-6
        end
        @test length(curves)==1
        _,xyz,_=g.model.mesh.get_nodes(1,only(curves)[2],true)
        widths=diff(sort(reshape(xyz,3,:)[2,:]))
        @test first>0
        @test widths[end]<=first*(1+1e-8)
        @test widths[end]<widths[1]
        # Floating lines never share a node with a two-dimensional element.
        _,_,blocks=g.model.mesh.get_elements(2);nodes=Set(vcat(blocks...))
        lines=g.model.get_entities_for_physical_group(1,Int(only(g.parser.get_number("MEASUREMENT_LINE"))))
        @test !isempty(lines)
        form=Formulation(:LineCableModelsFEM;options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        model=N.FEM._resolved_fem_model(problem,form)
        ratios=N.FEM._inspect_loaded_mesh(model,"native measurement mesh")
        # Allow twice the 0.25 target because Delaunay edges can be smaller than their targets.
        @test all(r -> isfinite(r) && 0<r<=.5,ratios)
        for curve in lines
            tags,_,_=g.model.mesh.get_nodes(1,curve,true)
            @test isempty(intersect(nodes,Set(tags)))
        end
    end
end

@testitem "FEM / interface seeds merge near-duplicate cable abscissae" tags=[:extension] setup=[NativeFEMFixtures] begin
    N=NativeFEMFixtures
    p=N.problem(;positions=[(0.,1.),(4e-15,1.5)],frequency=1e6,rho=100.)
    N.geometry(p) do g,log
        x=g.parser.get_number("FEMInterfaceX")
        tol=only(g.parser.get_number("FEMInterfaceTolerance"))
        @test count(v->abs(v)<tol,x)==1
        @test minimum(diff(x))>=tol
        @test !any(l->occursin("closer than the geometrical tolerance",l),log)
        for curve in Int.(g.parser.get_number("FEMFiniteInterface"))
            endpoints=g.model.get_boundary([(1,curve)],false,false,false)
            xyz=[g.model.get_value(dim,tag,Float64[]) for (dim,tag) in endpoints]
            @test length(xyz)==2
            @test abs(xyz[2][1]-xyz[1][1])>=tol
        end
    end
end

@testitem "FEM / mesh fields avoid MathEval reentry and mesh with a deadline" tags=[:extension] setup=[NativeFEMFixtures] begin
    using Gmsh
    N=NativeFEMFixtures
    # Inspect dependencies without evaluating fields: MathEval's lock is non-reentrant.
    for (f,rho,positions) in ((17_782_794.100389227,1000.,[(0.,1.),(1.,-1.)]),
            (1e6,.1,[(0.,1.),(1.,-1.)]),(50.,1000.,[(0.,-1.),(1.,-1.)]))
        p=N.problem(;frequency=f,rho,eps_r=1.,radius=.0425,positions)
        N.geometry(p;options=(mesh_size_factor=1.,)) do g,log
            fields=g.model.mesh.field
            types=Dict(Int(tag)=>fields.get_type(tag) for tag in fields.list())
            dependencies=Dict{Int,Vector{Int}}()
            for (tag,kind) in types
                dependencies[tag]=if kind=="MathEval"
                    [parse(Int,m.captures[1]) for m in eachmatch(r"\bF(\d+)\b",fields.get_string(tag,"F"))]
                elseif kind in ("Min","Max")
                    Int.(fields.get_numbers(tag,"FieldsList"))
                elseif kind in ("Threshold","Restrict")
                    [Int(fields.get_number(tag,"InField"))]
                else
                    Int[]
                end
            end
            for (tag,kind) in types
                kind=="MathEval" || continue
                pending=copy(dependencies[tag]);visited=Set{Int}()
                while !isempty(pending)
                    dependency=pop!(pending)
                    dependency in visited && continue
                    push!(visited,dependency)
                    @test types[dependency]!="MathEval"
                    append!(pending,dependencies[dependency])
                end
            end
            if f==17_782_794.100389227
                @test only(g.parser.get_number("FEMEarthLayerActive"))==0
                @test only(g.parser.get_number("MeshWaveEarth"))<only(g.parser.get_number("MeshRemoteEarth"))
            end
        end
    end
    p=N.problem(;frequency=17_782_794.100389227,rho=1000.,eps_r=1.,
        radius=.0425,positions=[(0.,1.),(1.,-1.)])
    N.bundle(p;options=(mesh_size_factor=1.,)) do directory,entry
        geo=replace(entry,r"\.pro$"=>".geo");mesh=joinpath(directory,"deadline.msh")
        open(joinpath(directory,"mesh.log"),"w") do io
            command=`$(Gmsh.gmsh_jll.gmsh()) $geo -2 -o $mesh -nt 1 -v 2`
            process=run(pipeline(command;stdout=io,stderr=io);wait=false)
            completed=timedwait(() -> process_exited(process),120.;pollint=.1)===:ok
            if !completed
                kill(process,Base.SIGKILL)
            end
            wait(process)
            @test completed
            @test success(process)
            @test isfile(mesh) && filesize(mesh)>0
        end
    end
end

@testitem "FEM / parallel mesh launches preserve sequential mesh snapshots" tags=[:extension] setup=[NativeFEMFixtures] begin
    using SHA, JSON3
    N=NativeFEMFixtures;E=N.FEM
    seed=N.problem(;frequency=.1,rho=.1,eps_r=1.,positions=[(0.,1.),(1.,-1.)])
    form=Formulation(:LineCableModelsFEM;options=(reduce_bundle=false,
        kron_reduction=false,ideal_transposition=false))
    mktempdir() do root
        for (label,frequencies,policy) in (("distinct",[.1,.2],:remesh),("shared",[.1,.1],:reuse),
                ("shared-remesh",[.1,.1],:remesh))
            problem=LineParametersProblem(seed.system;frequencies,
                temperature=seed.temperature,earth_props=seed.earth_props)
            records=[]
            for workers in (1,2)
                directory=joinpath(root,label,string(workers));run=E._create_run(directory,problem.system.system_id)
                execution=E.computation_options(E.LineCableModelsFEM,ComputationOptions(;
                    frequency_workers=workers,mesh_policy=policy,gmsh_verbosity=0))
                model=E._resolved_fem_model(problem,form);session=E._start_gmsh(0)
                try
                    E._prepare_run_inputs!(run,model,execution)
                    physical=E._build_physical_geometry!(model,"parallel-mesh-test")
                    E._write_physical_geometry(joinpath(run.path,"input","physical.geo"),model,physical)
                    E._write_native_mesh_entry(joinpath(run.path,"input","model.geo"),"model_data.pro","physical.geo","getdp")
                    paths=E._select_meshes!(run,model,execution)
                    push!(records,(meshes=[bytes2hex(open(sha256,path)) for path in paths],
                        sidecars=[read(replace(path,".msh"=>".json"),String) for path in paths]))
                    @test Set(readdir(joinpath(run.path,"mesh")))==Set(["frequency_$(lpad(i,4,'0')).$ext" for i in 1:2 for ext in ("msh","json")])
                    @test !ispath(joinpath(directory,"meshes"))
                    if startswith(label,"shared")
                        @test records[end].meshes[1]==records[end].meshes[2]
                        @test JSON3.read(records[end].sidecars[2]).source=="shared"
                        if policy===:reuse && workers==2
                            retained=read(paths[2])
                            log=joinpath(run.path,"logs","frequency_0001-gmsh.log")
                            previous_log=stat(log).mtime
                            rm(paths[1]); rm(replace(paths[1],".msh"=>".json"))
                            reopened=E._select_meshes!(run,model,execution)
                            @test read(reopened[1])==retained
                            @test stat(log).mtime==previous_log
                            @test JSON3.read(read(replace(reopened[1],".msh"=>".json"),String)).source=="shared"
                            @test JSON3.read(read(replace(reopened[2],".msh"=>".json"),String)).source=="resume"
                        end
                    end
                finally
                    E._finish_gmsh(session)
                end
            end
            @test records[1].meshes==records[2].meshes
            @test records[1].sidecars==records[2].sidecars
        end
    end
end

@testitem "FEM / supplied run meshes match fingerprints independently of frequency order" tags=[:extension] setup=[NativeFEMFixtures,TemporaryFEMRuntime] begin
    using Gmsh, JSON3
    N=NativeFEMFixtures
    cd(fem_test_runtime_directory)
    try
        seed=N.problem(;frequency=50.,rho=100.,eps_r=1.,radius=.005)
        problem(frequencies)=LineParametersProblem(seed.system;temperature=20.,frequencies,
            earth_props=seed.earth_props)
        form=Formulation(:LineCableModelsFEM;options=(reduce_bundle=false,
            kron_reduction=false,ideal_transposition=false))
        function run(frequencies;kwargs...)
            result=compute(problem(frequencies),form;options=(keep_run_directory=true,
                frequency_workers=2,verbosity=(default=0,),gmsh_verbosity=0,getdp_verbosity=0,kwargs...))
            record=details(result).data.fem
            path=record.run.run_directory
            sidecars=[JSON3.read(read(joinpath(path,"mesh","frequency_$(lpad(i,4,'0')).json"),String))
                for i in eachindex(frequencies)]
            meshes=[read(joinpath(path,"mesh","frequency_$(lpad(i,4,'0')).msh"))
                for i in eachindex(frequencies)]
            (;result,record,path,sidecars,meshes)
        end
        original=run([50.,1000.])
        for extension in ("msh","json")
            first=joinpath(original.path,"mesh","frequency_0001.$extension")
            second=joinpath(original.path,"mesh","frequency_0002.$extension")
            first_bytes=read(first); second_bytes=read(second)
            write(first,second_bytes); write(second,first_bytes)
        end
        reused=run([50.,1000.];mesh_path=original.path,overrides=(PetscPrealloc=128,))
        @test all(s->s.source=="previous_run",reused.sidecars)
        @test reused.meshes==original.meshes
        @test Y(reused.result)==Y(original.result)
        @test all(s->dirname(dirname(String(s.source_mesh)))==original.path,reused.sidecars)
        coarser=run([50.,1000.];mesh_path=original.path,mesh_size_factor=2.)
        @test all(s->s.source=="generated",coarser.sidecars)
        @test all(coarser.meshes[i]!=original.meshes[i] for i in 1:2)
        reordered=run([10.,50.,1000.];mesh_path=original.path)
        @test [String(s.source) for s in reordered.sidecars]==["generated","previous_run","previous_run"]
        @test reordered.meshes[2:3]==original.meshes
        @test Y(reordered.result)[:,:,2:3]==Y(original.result)
        subset=run([1000.];mesh_path=original.path)
        @test only(subset.sidecars).source=="previous_run"
        @test only(subset.meshes)==original.meshes[2]
        rm(original.path;recursive=true)
        for adopted in (reused,reordered,subset),i in eachindex(adopted.meshes)
            @test read(joinpath(adopted.path,"mesh","frequency_$(lpad(i,4,'0')).msh"))==adopted.meshes[i]
        end
        empty_run=mktempdir(); mkpath(joinpath(empty_run,"mesh"))
        for path in (joinpath(fem_test_runtime_directory,"missing-run"),mktempdir(),empty_run)
            err=try
                N.FEM.computation_options(N.FEM.LineCableModelsFEM,ComputationOptions(mesh_path=path))
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(path,sprint(showerror,err))
            isdir(path) && rm(path;recursive=true)
        end
        @test_throws ArgumentError N.FEM.computation_options(N.FEM.LineCableModelsFEM,
            ComputationOptions(mesh_path=reused.path,mesh_policy=:remesh))
    finally
        cd(fem_test_working_directory)
        rm(fem_test_runtime_directory;recursive=true)
    end
end
