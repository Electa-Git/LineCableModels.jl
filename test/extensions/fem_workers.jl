@testitem "Gmsh FEM / isolated workers, checkpoints and recovery" tags=[:extension] begin
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
root=pathlib.Path(settings['RunDirectory'])
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
            form = Formulation(:LineCableModelsFEM; options=(ideal_transposition=false,))
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
                E._prepare_run_inputs!(run,model,execution_options)
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
            @test E._assert_no_live_attempts(run)===nothing
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
            serial_form=Formulation(:LineCableModelsFEM;options=form.options)
            serial_form_controls = merge(form_controls,(frequency_workers=1,))
            E._run_getdp!(serial,model,serial_form, computation_options(LineCableModelsFEM, ComputationOptions(serial_form_controls)),meshes)
            @test E._parse_scan(serial,model,serial_form, computation_options(LineCableModelsFEM, ComputationOptions(serial_form_controls))).Z==scan.Z
            # A failed basis does not get a checkpoint and must be retried.
            broken=fresh();write(config,JSON3.write((delay=0.01,fail=true)))
            @test_throws LineCableModelsFEMError E._run_getdp!(broken,model,form, execution_options,meshes)
            @test isfile(E._column_paths(broken.path,1,1,false).checkpoint)
            @test !isfile(E._column_paths(broken.path,1,2,false).checkpoint)
            @test E._assert_no_live_attempts(broken)===nothing
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
            @test E._assert_no_live_attempts(unsupported)===nothing
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
            @test E._assert_no_live_attempts(stopped)===nothing
            @test_throws LineCableModelsFEMError E._parse_scan(stopped,model,form, execution_options)
            # A surviving child from a crashed coordinator blocks retry.
            orphan=joinpath(stopped.path,"work","orphan");mkpath(orphan)
            E._write_json_atomic(joinpath(orphan,"attempt.json"),
                (state="running",pid=getpid(),process_token=E._process_token(getpid())))
            @test_throws LineCableModelsFEMError E._assert_no_live_attempts(stopped)
            E._write_json_atomic(joinpath(orphan,"attempt.json"),
                (state="running",pid=getpid(),process_token="old boot or reused PID"))
            Sys.islinux() && @test E._assert_no_live_attempts(stopped)===nothing
        end
    else
        @test_skip "The worker process fixture requires Unix and python3"
    end
end

@testitem "FEM / passive native PML observations retain extent and cutoff flags" tags=[:extension] begin
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
