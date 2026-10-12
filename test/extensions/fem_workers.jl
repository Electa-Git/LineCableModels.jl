@testitem "Gmsh FEM / isolated workers, checkpoints and recovery" tags=[:extension, :fem] begin
    using Gmsh, JSON3, SHA
    const LineCableModelsFEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEM
    const LineCableModelsFEMError = Base.get_extension(LineCableModels, :LineCableModelsGmshExt).LineCableModelsFEMError
    E = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    python = Sys.which("python3")
    if Sys.isunix() && python !== nothing
        copper = Material(kind=:conductor, rho=1.72e-8)
        design = build(CableDesign, "worker-contract", terminal(:core, core(copper; r=0.005)))
        system = build(LineCableSystem, [design,design], [(0.0,-0.1),(0.1,-0.1)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
        problem = LineParametersProblem(system; frequencies=[50.0,1000.0],
            earth_props=homogeneous(rho=100.0,eps_r=10.0))
        mktempdir(;prefix="fem workers with spaces ") do root
            config = joinpath(root,"control.json")
            executable = joinpath(root,"getdp fixture")
            write(executable,"#!$python\ncontrol_path = " * JSON3.write(config) * "\n" * raw"""
import json, sys, pathlib, re, time, os
if sys.argv[1:] == ['-info']:
    print('GetDP Version 3.6.0 worker protocol fixture'); sys.exit(0)
args=sys.argv[1:]; settings={}
for i,a in enumerate(args):
    if a in ('-setstring','-setnumber'): settings[args[i+1]]=args[i+2]
control=json.loads(pathlib.Path(control_path).read_text())
if control.get('unsupported'):
    print('Unknown operation GenerateRHSGroup',flush=True);sys.exit(2)
if control.get('signal'):
    import resource
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    os.kill(os.getpid(),control['signal'])
root=pathlib.Path(settings['RunDirectory']).resolve()
f=int(settings['FrequencyIndex']); data=pathlib.Path(settings['ModelDataPath']).read_text()
hz=float(re.search(r'Frequencies\(\) = \{([^}]*)\}',data)[1].split(',')[f-1])
bases=[int(x) for x in re.search(r'\{([^}]*)\}',pathlib.Path(settings['BasisListPath']).read_text())[1].split(',')]
assert os.getcwd()==str(root)
assert all(os.environ[x]=='1' for x in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'))
(root/'observed.json').write_text(json.dumps({'frequency':f,'bases':bases,'pid':os.getpid()}))
for b in bases:
    time.sleep(control.get('delay',0.0)*(3-f))
    stem=root/'results'/'jobs'/f'getdp-f{f:04d}-b{b:04d}'
    for q in ('Z','P','Pscalar'):
        text=''.join(f'{f}\t{hz:.17g}\t{r}\t{b}\t{100*f+10*r+b}\t{r-b}\n' for r in (1,2))
        pathlib.Path(str(stem)+f'-{q}.tsv').write_text(text)
    obs=[f,hz,9.,0,0,0,0,9.,9.,9.,9.,0,48,48,48,1.,1.,0.,0,0.,0,0,0,0,1]
    (stem.parent/f'pml-f{f:04d}.tsv').write_text('\t'.join(map(str,obs))+'\n')
    pathlib.Path(str(stem)+'-timing.tsv').write_text(f'{f}\t{b}\t0\t0\t0\t0\t{int(b==bases[0])}\n')
    if control.get('fail') and f==1 and b==2:
        print('deliberate failure after partial column',flush=True); sys.exit(17)
    pathlib.Path(str(stem)+'.done').write_text(f'2\t{f}\t{hz:.17g}\t{b}\t2\t0\n')
""")
            chmod(executable,0o700)
            form = Formulation(:fem; options=(ideal_transposition=false,))
            form_controls = (getdp_executable=executable,frequency_workers=2,solver_threads=1,gmsh_verbosity=0)
            execution_options = computation_options(LineCableModelsFEM, ComputationOptions(form_controls))
            model=E._resolved_fem_model(problem,form)
            inputs=E._fem_input_record(model,form, execution_options)
            execution=inputs.execution
            meshes=map(1:2) do i
                path=joinpath(root,"mesh $i.msh")
                write(path,"mesh $i")
                path
            end
            function fresh()
                run=E._create_run(root,problem.system.system_id)
                E._write_json_atomic(joinpath(run.path,"input/computation.json"),inputs)
                E._write_run_inputs!(run,model,execution_options)
                run
            end
            write(config,JSON3.write((delay=0.02,fail=false)))
            run=fresh()
            owner=E._claim_run(run)
            @test_throws LineCableModelsFEMError E._claim_run(run)
            lock_probe=raw"""
import fcntl,sys
with open(sys.argv[1],'a+') as f:
    try: fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB);print('acquired')
    except BlockingIOError: print('busy')
"""
            lock_command=`$python -c $lock_probe $(joinpath(run.path,"coordinator.lock"))`
            @test strip(read(lock_command,String))=="busy"
            E._release_run(owner)
            @test strip(read(lock_command,String))=="acquired"
            E._release_run(E._claim_run(run))
            E._run_getdp!(run,model,form, execution_options,meshes)
            @test (run.getdp_invocations,run.completed_columns,run.completed_frequencies)==(2,4,2)
            scan=E._parse_scan(run,model,form, execution_options)
            @test real.(scan.Z[:,:,1])==[111 112;121 122]
            @test real.(scan.Z[:,:,2])==[211 212;221 222]
            digest=bytes2hex(open(sha256,meshes[1]))
            @test E._valid_column_checkpoint(run.path,1,50.0,2,2,false,digest)
            @test !E._valid_column_checkpoint(run.path,1,50.0,2,2,false,"different mesh")
            @test !E._valid_column_checkpoint(run.path,1,50.0,2,2,true,digest)
            paths=E._column_paths(run.path,1,2,false)
            before=read(paths.Z)
            write(paths.Z,"broken\n")
            @test !E._valid_column_checkpoint(run.path,1,50.0,2,2,false,digest)
            # Successful scratch is deleted; a corrupt adopted column is solved again.
            E._run_getdp!(run,model,form, execution_options,meshes)
            @test run.getdp_invocations==3
            @test read(paths.Z)==before
            # A corrupt manifest cannot prevent adopting other attempts.
            for (kind,content) in (("object","{}"),("array","[]"),("scalar","1"))
                directory=joinpath(run.path,"work","malformed_"*kind)
                mkpath(directory)
                write(joinpath(directory,"attempt.json"),content)
            end
            @test E._live_attempt(run)===nothing
            # Lose both the canonical column and its completed attempt marker.
            rm(paths.checkpoint)
            for attempt in filter(isdir,readdir(joinpath(run.path,"work");join=true))
                rm(E._column_paths(attempt,1,2,false).marker;force=true)
            end
            E._run_getdp!(run,model,form, execution_options,meshes)
            @test run.getdp_invocations==4
            @test !any(name->startswith(name,"frequency_"),readdir(joinpath(run.path,"work")))
            records=[JSON3.read(read(joinpath(d,"attempt.json"),String)) for d in
                filter(isdir,readdir(joinpath(run.path,"logs");join=true))]
            @test any(record->record.requested_bases==[2],records)
            @test E._parse_scan(run,model,form, execution_options).Z==scan.Z
            # Concurrency is scheduling metadata; thread count is numerical provenance.
            one=merge(inputs,(execution=merge(execution,(frequency_workers=1,)),))
            threads=merge(inputs,(execution=merge(execution,(solver_threads=2,)),))
            @test E._resume_inputs_match(run.path,model,one)
            @test !E._resume_inputs_match(run.path,model,threads)
            serial=fresh()
            serial_form=Formulation(:fem;options=form.options)
            serial_form_controls = merge(form_controls,(frequency_workers=1,))
            E._run_getdp!(serial,model,serial_form, computation_options(LineCableModelsFEM, ComputationOptions(serial_form_controls)),meshes)
            @test E._parse_scan(serial,model,serial_form, computation_options(LineCableModelsFEM, ComputationOptions(serial_form_controls))).Z==scan.Z
            # A failed basis does not get a checkpoint and must be retried.
            broken=fresh();write(config,JSON3.write((delay=0.01,fail=true)))
            @test_throws LineCableModelsFEMError E._run_getdp!(broken,model,form, execution_options,meshes)
            @test isfile(E._column_paths(broken.path,1,1,false).checkpoint)
            @test !isfile(E._column_paths(broken.path,1,2,false).checkpoint)
            @test E._live_attempt(broken)===nothing
            write(config,JSON3.write((delay=0.0,fail=false)))
            E._run_getdp!(broken,model,form, execution_options,meshes)
            @test E._parse_scan(broken,model,form, execution_options).Z==scan.Z
            unsupported=fresh();write(config,JSON3.write((unsupported=true,)))
            failure=try
                E._run_getdp!(unsupported,model,form, execution_options,meshes)
                nothing
            catch exception
                exception
            end
            @test failure isa LineCableModelsFEMError
            @test failure.field===:capability
            @test E._live_attempt(unsupported)===nothing
            signalled=fresh();write(config,JSON3.write((signal=6,)))
            failure=try
                E._run_getdp!(signalled,model,form,execution_options,meshes)
                nothing
            catch exception
                exception
            end
            @test failure isa LineCableModelsFEMError
            @test occursin("terminated by signal 6",sprint(showerror,failure))
            @test !occursin("(exit 0)",sprint(showerror,failure))
            attempts=[JSON3.read(read(joinpath(d,"attempt.json"),String)) for d in
                filter(isdir,readdir(joinpath(signalled.path,"work");join=true))]
            @test any(record->record.signal==6,attempts)
            @test E._live_attempt(signalled)===nothing
            # Cancellation kills and reaps workers and leaves no completion claim.
            stopped=fresh();write(config,JSON3.write((delay=2.0,fail=false)))
            task = @async try
                E._run_getdp!(stopped, model, form, execution_options, meshes)
            catch exception
                exception
            end
            deadline = time() + 30
            while stopped.getdp_invocations == 0 && !istaskdone(task)
                time() < deadline || error("worker launch timed out")
                sleep(0.01)
            end
            schedule(task, InterruptException(); error=true)
            @test fetch(task) isa InterruptException
            @test stopped.state===E.cancelled
            @test E._live_attempt(stopped)===nothing
            @test_throws LineCableModelsFEMError E._parse_scan(stopped,model,form, execution_options)
            # A surviving child from a crashed coordinator blocks retry.
            orphan=joinpath(stopped.path,"work","orphan");mkpath(orphan)
            E._write_json_atomic(joinpath(orphan,"attempt.json"),
                (state="running",pid=getpid(),process_token=E._process_token(getpid())))
            @test E._live_attempt(stopped) == (pid=getpid(),directory=orphan)
            E._write_json_atomic(joinpath(orphan,"attempt.json"),
                (state="running",pid=getpid(),process_token="old boot or reused PID"))
            Sys.islinux() && @test E._live_attempt(stopped)===nothing
        end
    else
        @test_skip "The worker process fixture requires Unix and python3"
    end
end

@testitem "FEM / passive native PML observations retain extent and cutoff flags" tags=[:extension, :fem] begin
    using Gmsh
    FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    mktempdir() do directory
        path = joinpath(directory,"pml.tsv")
        values = [1., .1, 3.45, 1., 0., 0., 0., 0., 10., 0., 10., 1.,
            16., 16., 16., 0., 0., 0., 0., 0., 0., 1., 1., 0., 0.]
        write(path,join(values,'\t')*"\n")
        observation = FEM._pml_observation(path,1,.1)
        @test observation !== nothing
        @test observation.quasi_static_cap_flags == (true,true,false)
        @test observation.net_top_exponent == 0
        @test !observation.air_excited
        @test !hasproperty(observation,:g_not_qualified)
        @test !hasproperty(observation,:attenuation_below_target)
        values[6] = 1
        values[19] = 1
        values[25] = 1
        write(path,join(values,'\t')*"\n")
        observation = FEM._pml_observation(path,1,.1)
        @test observation.cutoff_flags == (true,false)
        @test observation.earth_sizing_ceiling_active
        @test observation.air_excited
    end
end

@testitem "Gmsh FEM / deep managed paths stay below native argument limits" tags=[:extension,:integration,:fem_numerical,:fem] setup=[TemporaryFEMRuntime] begin
    using Gmsh, JSON3
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    deep=joinpath(fem_test_runtime_directory,repeat("deep-",16),repeat("folder-",14))
    mkpath(deep)
    cd(deep)
    try
        wire=build(CableDesign,"path-compliance",terminal(:core,
            core(Material(kind=:conductor,rho=1.72e-8);r=.01)))
        formulation=Formulation(:fem;options=(
            reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        options=(overrides=(PmlSideLayers=8,PmlTopLayers=8,PmlBottomLayers=8,),
            frequency_workers=1,solver_threads=1,mesh_policy=:remesh,
            keep_run_directory=true,plot_field_maps=true,gmsh_verbosity=0,getdp_verbosity=4)
        mesh_inputs=Dict{String,String}[]
        for system_id in ("short",repeat("long-system-",12))
            system=build(LineCableSystem,[wire,wire],[(0.,1.),(1.,1.)];
                connections=[Dict(:core=>1),Dict(:core=>2)],system_id)
            problem=LineParametersProblem(system;frequencies=[1e6],
                earth_props=homogeneous(rho=1000.,eps_r=12.))
            result=compute(problem,formulation;options)
            @test pwd()==deep
            @test all(isfinite,Z(result)) && all(isfinite,Y(result))
            run_path=details(result).data.fem.run.run_directory
            @test ncodeunits(run_path)>300
            records=[JSON3.read(read(joinpath(directory,"attempt.json"),String))
                for directory in filter(isdir,readdir(joinpath(run_path,"logs");join=true))]
            @test length(records)==1
            arguments=String.(only(records).command)
            @test all(argument->ncodeunits(argument)<256,arguments)
            @test all(!isabspath,arguments[2:end])
            strings=Dict(arguments[i+1]=>arguments[i+2] for i in eachindex(arguments)
                if arguments[i]=="-setstring")
            @test strings==Dict("ModelDataPath"=>"model_data.pro","RunDirectory"=>".",
                "RawDirectory"=>"results","BasisListPath"=>"bases.pro")
            @test ncodeunits(joinpath(run_path,"work","frequency_0001","bases.pro"))>300
            mesh_record=JSON3.read(first(readlines(joinpath(run_path,"logs","frequency_0001-gmsh.log"))))
            @test all(argument->ncodeunits(argument)<256,mesh_record.command)
            metadata_index=findfirst(==("MeshMetadataPath"),mesh_record.command)
            @test mesh_record.command[metadata_index+1]=="mesh-metadata.txt"
            push!(mesh_inputs,Dict(name=>read(joinpath(run_path,"input",name),String)
                for name in ("physical.geo","model_data.pro","model.geo")))
            @info "Native path compliance" absolute_run_bytes=ncodeunits(run_path) maximum_argument_bytes=maximum(ncodeunits,arguments)
        end
        @test mesh_inputs[1]==mesh_inputs[2]
        @test FEM._check_native_path("BasisListPath",repeat("x",255))==repeat("x",255)
        failure=try
            FEM._check_native_path("BasisListPath",repeat("x",256))
            nothing
        catch exception
            exception
        end
        @test failure isa FEM.LineCableModelsFEMError
        @test occursin("BasisListPath is 256 bytes long",sprint(showerror,failure))
        @test_throws FEM.LineCableModelsFEMError FEM._check_native_path("BasisListPath",repeat("é",128))
    finally
        cd(fem_test_working_directory)
        rm(fem_test_runtime_directory;recursive=true)
    end
end
