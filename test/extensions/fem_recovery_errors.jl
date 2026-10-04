@testmodule FEMRecoveryFaults begin
    import LineCableModels, Gmsh, JSON3
    using SHA: sha256
    const E = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const ImportExport = LineCableModels.ImportExport
    const FEMResolvedModel = NamedTuple
    const FEMScan = NamedTuple
    const FEMRun = E.FEMRun
    const FEMActiveWorker = E.FEMActiveWorker
    const created, completed, running = E.created, E.completed, E.running
    const _fem_error = E._fem_error
    const _create_run = E._create_run
    const _resume_value_matches = E._resume_value_matches
    const _column_paths = E._column_paths
    const _column_files = E._column_files
    const _column_timing = E._column_timing
    const _write_json_atomic = E._write_json_atomic
    const _attempt_record = E._attempt_record
    const _transition! = E._transition!
    const read_fault = Ref{Any}(nothing)
    const fault_path = Ref("")
    const fault_after = Ref(1)
    const read_count = Ref(0)
    const raw_fault = Ref{Any}(nothing)
    const pid_fault = Ref{Any}(nothing)
    const wait_for_exit = Ref(false)
    const child = Ref{Base.Process}()
    const log = Ref{IOStream}()

    function read(path::String, ::Type{String})
        if path == fault_path[]
            read_count[] += 1
            if read_count[] >= fault_after[] && read_fault[] !== nothing
                read_fault[] isa String && return read_fault[]
                throw(read_fault[])
            end
        end
        return Base.read(path, String)
    end
    function getpid(process::Base.Process)
        child[] = process
        wait_for_exit[] && wait(process)
        pid_fault[] === nothing || throw(pid_fault[])
        return Base.getpid(process)
    end
    function open(path::String, mode::String)
        log[] = Base.open(path, mode)
    end
    open(f::Function, path::String) = Base.open(f, path)
    function _valid_job_raw(arguments...)
        raw_fault[] === nothing || throw(raw_fault[])
        return E._valid_job_raw(arguments...)
    end

    # Isolate file and process faults without changing package or Base methods. A name
    # selects every definition of that function; a signature selects one method.
    for (file, selected) in (("compute.jl", (:_resume_inputs_match, :_resume_run)),
            ("results.jl", (:(validate(scan::FEMScan, run::FEMRun)),)),
            ("workers.jl", (:_valid_column_checkpoint, :_process_token, :_start_worker!)))
        Base.include(@__MODULE__, joinpath(@__DIR__, "../../ext/LineCableModelsGmshExt", file)) do expression
            expression isa Expr && expression.head === :function &&
                expression.args[1] isa Expr && expression.args[1].head === :call &&
                (expression.args[1].args[1] in selected || expression.args[1] in selected) ?
                expression : nothing
        end
    end
end

@testitem "Gmsh FEM / recovery rejects invalid files and propagates unexpected failures" tags=[:extension, :fem] setup=[FEMRecoveryFaults] begin
    using JSON3
    const F = FEMRecoveryFaults
    model = (problem=(system=(system_id="recovery-fixture",),),)
    inputs = (owned_gmsh=true, schema_version=7, solver_protocol=3,
        getdp_identity=nothing, execution=(ui=false, mesh_mode=:reuse))
    unexpected = (InterruptException(), ErrorException("unexpected reader bug"),
        MethodError(sin, (nothing,)))
    caught(f) = try
        f()
    catch exception
        exception
    end

    mktempdir() do root
        run = F.E._create_run(root)
        F.E._write_json_atomic(joinpath(run.path, "input", "problem.json"),
            LineCableModels.ImportExport.serialize_value(model.problem))
        F.E._write_json_atomic(joinpath(run.path, "input", "computation.json"), inputs)
        state = joinpath(run.path, "run.json")
        original = read(state, String)
        @test F._resume_inputs_match(run.path, model, inputs)
        # JSON snapshots keep their lossless integer view. Input records on
        # both sides are parsed alike and remain insensitive to object order.
        nested = merge(inputs, (extra=(large=Int64(9007199254740993),
            values=[(a=1, b=2), (a=3, b=4)]),))
        computation = joinpath(run.path, "input", "computation.json")
        reordered = merge(nested, (extra=(values=[(b=2, a=1), (b=4, a=3)],
            large=Int64(9007199254740993)),))
        F.E._write_json_atomic(computation, reordered)
        @test F._resume_inputs_match(run.path, model, nested)
        changed = merge(nested, (extra=merge(nested.extra, (large=Int64(9007199254740992),)),))
        @test !F._resume_inputs_match(run.path, model, changed)
        F.E._write_json_atomic(computation, inputs)
        @test F._resume_run(root, run.path, model, inputs).getdp_invocations == 0
        for malformed in ("{", "[]", "null", "{}", "{\"state\":42}")
            write(state, malformed)
            @test !F._resume_inputs_match(run.path, model, inputs)
        end
        write(state, original)
        F.fault_path[] = state
        for exception in unexpected
            F.read_fault[] = exception
            @test caught(() -> F._resume_inputs_match(run.path, model, inputs)) === exception
        end
        F.read_fault[] = SystemError("read fixture", 2)
        @test !F._resume_inputs_match(run.path, model, inputs)

        # The first metadata read selects the run. A later read must not invent
        # fresh metadata if the selected file disappears or becomes malformed.
        F.fault_after[] = 2
        for fault in (SystemError("selected file disappeared", 2), "{")
            F.read_count[] = 0
            F.read_fault[] = fault
            exception = caught(() -> F._resume_run(root, run.path, model, inputs))
            @test exception isa LineCableModelsFEMError
            @test exception.run_directory == run.path
            @test occursin(state, exception.message)
            @test occursin(fault isa String ? "invalid JSON" : "selected file disappeared", exception.message)
            @test read(state, String) == original
        end
        for exception in unexpected
            F.read_count[] = 0
            F.read_fault[] = exception
            @test caught(() -> F._resume_run(root, run.path, model, inputs)) === exception
        end
        F.read_count[] = 0
        F.read_fault[] = "{}"
        @test caught(() -> F._resume_run(root, run.path, model, inputs)) isa LineCableModelsFEMError

        F.fault_after[] = 1
        F.fault_path[] = joinpath(run.path, "raw", "checksums.json")
        for fault in (SystemError("checksum file disappeared", 2), "{")
            F.read_fault[] = fault
            exception = caught(() -> F.validate((map_paths=String[],), run))
            @test exception isa LineCableModelsFEMError
            @test exception.field === :checksum
            @test occursin(F.fault_path[], exception.message)
            @test occursin(fault isa String ? "invalid JSON" : "checksum file disappeared", exception.message)
        end
        for exception in unexpected
            F.read_fault[] = exception
            @test caught(() -> F.validate((map_paths=String[],), run)) === exception
        end

        paths = F._column_paths(run.path, 1, 1, false)
        mkpath(dirname(paths.checkpoint))
        write(paths.checkpoint, JSON3.write((protocol=2, frequency_index=1, basis=1,
            terminals=1, frequency_hz=50.0, plot_field_maps=false,
            mesh_digest="mesh", checksums=Dict())))
        F.fault_path[] = paths.checkpoint
        for exception in unexpected
            F.read_fault[] = exception
            @test caught(() -> F._valid_column_checkpoint(run.path, 1, 50.0, 1, 1, false, "mesh")) === exception
        end
        F.read_fault[] = SystemError("checkpoint file disappeared", 2)
        @test !F._valid_column_checkpoint(run.path, 1, 50.0, 1, 1, false, "mesh")
        F.read_fault[] = nothing
        for exception in (unexpected..., ArgumentError("raw validator bug"))
            F.raw_fault[] = exception
            @test caught(() -> F._valid_column_checkpoint(run.path, 1, 50.0, 1, 1, false, "mesh")) === exception
        end
        F.raw_fault[] = SystemError("raw file disappeared", 2)
        @test !F._valid_column_checkpoint(run.path, 1, 50.0, 1, 1, false, "mesh")
        F.raw_fault[] = nothing

        if Sys.islinux()
            F.fault_path[] = "/proc/$(getpid())/stat"
            for exception in unexpected
                F.read_fault[] = exception
                @test caught(() -> F._process_token(getpid())) === exception
            end
            F.read_fault[] = SystemError("process disappeared", 2)
            @test F._process_token(getpid()) === nothing
            F.read_fault[] = "incomplete stat"
            @test F._process_token(getpid()) === nothing
            F.read_fault[] = nothing
            @test F._process_token(getpid()) isa String
        end
        F.read_fault[] = nothing
    end
end

@testitem "Gmsh FEM / PID lookup only recovers exited children and cleans up failures" tags=[:extension, :fem] setup=[FEMRecoveryFaults] begin
    const F = FEMRecoveryFaults
    mktempdir() do root
        run = F.E._create_run(root)
        execution = (data=(solver_threads=1,),)
        exited = `$(Base.julia_cmd()) --startup-file=no --handle-signals=no --project=@stdlib -e 'exit()'`
        job = F.E.FEMFrequencyJob(1, 50.0, [1], root, "mesh", exited)
        F.wait_for_exit[] = true
        worker = F._start_worker!(run, job, execution)
        @test worker.pid == 0
        @test worker.process_token === nothing
        @test process_exited(worker.process)
        close(worker.log)
        F.wait_for_exit[] = false

        active = `$(Base.julia_cmd()) --startup-file=no --handle-signals=no --project=@stdlib -e 'sleep(30)'`
        job = F.E.FEMFrequencyJob(1, 50.0, [1], root, "mesh", active)
        for exception in (InterruptException(), ErrorException("PID lookup bug"),
                Base.IOError("unexpected PID error", 0))
            F.pid_fault[] = exception
            caught = try
                F._start_worker!(run, job, execution)
            catch error
                error
            end
            @test caught === exception
            @test process_exited(F.child[])
            @test !isopen(F.log[])
        end
        F.pid_fault[] = nothing
    end
end
