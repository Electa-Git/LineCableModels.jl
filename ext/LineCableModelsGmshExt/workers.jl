"""
    FEMFrequencyJob

Immutable process inputs for one mesh, frequency and requested terminal list.
"""
struct FEMFrequencyJob
    "Position in the input frequency vector."
    frequency_index::Int
    "Physical frequency \\[Hz\\]."
    frequency::Float64
    "Terminal columns requested from this process."
    bases::Vector{Int}
    "Exclusive attempt directory."
    directory::String
    "Digest of the mesh used by this attempt."
    mesh_digest::String
    "Standalone GetDP command with its own environment and working directory."
    command::Cmd
end

struct FEMActiveWorker{P}
    job::FEMFrequencyJob
    process::P
    pid::Int
    process_token::Union{Nothing, String}
    log::IOStream
    started::Float64
    started_ns::UInt64
    pending::Set{Int}
end

_column_stem(frequency::Int, basis::Int) = @sprintf("getdp-f%04d-b%04d", frequency, basis)

function _column_paths(root::String, frequency::Int, basis::Int, maps::Bool;
        physics::Symbol=Symbol("quasi-tem"))
    stem = _column_stem(frequency, basis)
    raw = _job_raw_paths(root, stem)
    return (; raw...,
        diagnostics = _quasi_full(physics) ?
            [joinpath(root, "raw", "jobs", "$stem-Pscalar.tsv")] : String[],
        timing = joinpath(root, "raw", "jobs", "$stem-timing.tsv"),
        marker = joinpath(root, "raw", "jobs", "$stem.done"),
        checkpoint = joinpath(root, "raw", "jobs", "$stem.json"),
        maps = maps ?
               [joinpath(root, "maps",
                    @sprintf("%s_f%04d_b%04d.pos", quantity, frequency, basis))
                for quantity in _field_quantities(physics)] : String[])
end

_column_files(paths) = [paths.Z, paths.P, paths.timing, paths.maps..., paths.diagnostics...]

function _valid_column_marker(path, frequency_index, frequency, basis, terminals, maps)
    isfile(path) || return false
    text = read(path, String)
    endswith(text, '\n') || return false
    rows = split(strip(text), '\n')
    length(rows) == 1 || return false
    fields = split(only(rows), '\t')
    length(fields) == 6 || return false
    values = tryparse.(Float64, fields)
    any(isnothing, values) && return false
    return values[[1, 2, 4, 5, 6]] == [2, frequency_index, basis, terminals, Int(maps)] &&
           isapprox(values[3], frequency; rtol = 16eps(Float64), atol = 0)
end

function _column_timing(path, frequency_index, basis)
    isfile(path) || return nothing
    rows = readlines(path)
    length(rows) == 1 || return nothing
    values = tryparse.(Float64, split(only(rows), '\t'))
    length(values) == 7 && all(v -> v !== nothing && isfinite(v), values) || return nothing
    values[1:2] == [frequency_index, basis] || return nothing
    all(>=(0), values[3:6]) && values[7] in (0, 1) || return nothing
    return (constraint_seconds = values[3], assembly_seconds = values[4],
        solve_seconds = values[5], output_seconds = values[6], factorized = Bool(values[7]))
end

function _complete_map(path)
    isfile(path) && filesize(path) > 3 || return false
    return open(path) do io
        seek(io, max(0, filesize(path) - 128))
        endswith(strip(read(io, String)), "};")
    end
end

function _valid_column_checkpoint(
        root, frequency_index, frequency, basis, terminals, maps, mesh_digest;
        physics::Symbol=Symbol("quasi-tem"))
    paths = _column_paths(root, frequency_index, basis, maps; physics)
    isfile(paths.checkpoint) || return false
    return try
        record = JSON3.read(read(paths.checkpoint, String))
        record.protocol == 2 && record.frequency_index == frequency_index &&
        record.basis == basis && record.terminals == terminals &&
        get(record, :physics, "quasi-tem") == String(physics) &&
        record.plot_field_maps == maps && record.mesh_digest == mesh_digest &&
        isapprox(record.frequency_hz, frequency; rtol = 16eps(Float64), atol = 0) ||
            return false
        _valid_job_raw(paths.Z, terminals, frequency_index, frequency, basis) &&
        _valid_job_raw(paths.P, terminals, frequency_index, frequency, basis) ||
            return false
        _column_timing(paths.timing, frequency_index, basis) === nothing && return false
        files = _column_files(paths)
        length(record.checksums) == length(files) || return false
        all(files) do file
            key = relpath(file, root)
            isfile(file) && haskey(record.checksums, key) &&
                record.checksums[key] == bytes2hex(open(sha256, file))
        end
    catch
        false
    end
end

function _copy_column_file(source, destination)
    mkpath(dirname(destination))
    temporary = tempname(dirname(destination))
    try
        cp(source, temporary)
        mv(temporary, destination; force = true)
    finally
        rm(temporary; force = true)
    end
end

function _adopt_column!(
        run, source, frequency_index, frequency, basis, terminals, maps, mesh_digest;
        physics::Symbol=Symbol("quasi-tem"))
    paths = _column_paths(source, frequency_index, basis, maps; physics)
    _valid_column_marker(
        paths.marker, frequency_index, frequency, basis, terminals, maps) || return false
    _valid_job_raw(paths.Z, terminals, frequency_index, frequency, basis) &&
    _valid_job_raw(paths.P, terminals, frequency_index, frequency, basis) || return false
    all(path -> _valid_job_raw(path, terminals, frequency_index, frequency, basis),
        paths.diagnostics) || return false
    timing = _column_timing(paths.timing, frequency_index, basis)
    timing === nothing && return false
    all(_complete_map, paths.maps) || return false
    destination = _column_paths(run.path, frequency_index, basis, maps; physics)
    for (input, output) in zip(_column_files(paths), _column_files(destination))
        _copy_column_file(input, output)
    end
    checksums = Dict(relpath(file, run.path) => bytes2hex(open(sha256, file))
    for file in _column_files(destination))
    _write_json_atomic(destination.checkpoint,
        (; protocol = 2, frequency_index, frequency_hz = frequency, basis, terminals,
            plot_field_maps = maps, physics, mesh_digest, timing, checksums,
            attempt = relpath(source, run.path)))
    return true
end

function _getdp_command(executable, model_path, mesh_path, run, formulation, execution, mesh_plan,
        bases, directory; reuse_factorization = true)
    basis_path = joinpath(directory, "bases.pro")
    write(basis_path, "RequestedBases() = $(_pro_array(bases));\n")
    verbosity = LineCableModels.performance_sample_active() ? 0 : execution.getdp_verbosity
    arguments = [executable, model_path, "-solve", "LineCableModelsFEMScan",
        "-setnumber", "Physics", string(_fem_physics_code(formulation)),
        "-msh", abspath(mesh_path), "-name", joinpath(directory, "solver"),
        "-v", string(verbosity),
        "-setstring", "ModelDataPath", joinpath(run.path, "input", "model_data.pro"),
        "-setstring", "RunDirectory", directory,
        "-setstring", "BasisListPath", basis_path,
        "-setnumber", "FrequencyIndex", string(mesh_plan.frequency_index),
        "-setnumber", "FrequencyHz", _pro_number(mesh_plan.frequency),
        "-setnumber", "Val_Rint", _pro_number(mesh_plan.domain_radius),
        "-setnumber", "Val_Rext", _pro_number(mesh_plan.shell_outer_radius),
        "-setnumber", "PlotFieldMaps", string(Int(execution.plot_field_maps)),
        "-setnumber", "ReuseFactorization", string(Int(reuse_factorization))]
    if _quasi_full(formulation.options.physics)
        append!(arguments, ["-setstring", "PathDataPath",
            _voltage_path_file(run, mesh_plan.frequency_index)])
    end
    if verbosity >= 4
        append!(arguments, [
            "-cpu", "-ksp_view", "-log_view", ":" * joinpath(directory, "petsc.log")])
    end
    threads = string(execution.solver_threads)
    return addenv(Cmd(Cmd(arguments); dir = directory), "OPENBLAS_NUM_THREADS" => threads,
        "OPENBLAS_DEFAULT_NUM_THREADS" => threads,
        "OMP_NUM_THREADS" => threads, "MKL_NUM_THREADS" => threads)
end

function _frequency_job(
        run, model, formulation, execution, executable, mesh, mesh_digest, frequency_index,
        bases; reuse_factorization = true)
    isempty(bases) && throw(ArgumentError("a frequency job requires at least one terminal"))
    all(b -> b in eachindex(model.terminal_ids), bases) && allunique(bases) ||
        throw(ArgumentError("frequency job terminal indices must be distinct and in range"))
    root = joinpath(run.path, "attempts")
    mkpath(root)
    directory = mktempdir(root; prefix = @sprintf("f%04d-", frequency_index), cleanup = false)
    for subdirectory in ("raw/jobs", "maps")
        mkpath(joinpath(directory, subdirectory))
    end
    plan = model.mesh_plans[frequency_index]
    command = _getdp_command(
        executable, _getdp_assets(joinpath(run.path, "input", "getdp")).model,
        mesh, run, formulation, execution, plan, bases, directory; reuse_factorization)
    return FEMFrequencyJob(frequency_index, Float64(plan.frequency), collect(bases),
        directory, mesh_digest, command)
end

function _attempt_record(job; kwargs...)
    return (;
        protocol = 2, frequency_index = job.frequency_index, frequency_hz = job.frequency,
        requested_bases = job.bases, mesh_digest = job.mesh_digest, command = collect(job.command), kwargs...)
end

function _start_worker!(run, job, execution)
    log = open(joinpath(job.directory, "getdp.log"), "w")
    started, started_ns = time(), time_ns()
    process = try
        Base.run(pipeline(job.command; stdout = log, stderr = log); wait = false)
    catch
        close(log)
        rethrow()
    end
    pid = try
        Int(getpid(process))
    catch
        0 # A process can exit before its handle's PID is queried.
    end
    process_token = _process_token(pid)
    worker = FEMActiveWorker(
        job, process, pid, process_token, log, started, started_ns, Set(job.bases))
    try
        _write_json_atomic(joinpath(job.directory, "attempt.json"),
            _attempt_record(
                job; state = "running", pid, process_token, started_unix_seconds = started,
                solver_threads = execution.solver_threads))
        run.getdp_invocations += 1
        _transition!(run, running, "frequency $(job.frequency_index) launched")
    catch
        process_running(process) && kill(process, 9) # SIGKILL, accepted by the public kill(process, signal) API
        wait(process)
        close(log)
        rethrow()
    end
    return worker
end

function _record_progress!(run, valid)
    run.completed_columns = count(valid)
    run.completed_frequencies = count(all, eachcol(valid))
    _transition!(run, running,
        "$(run.completed_columns)/$(length(valid)) terminal columns validated")
end

function _collect_worker_columns!(run, worker, valid, maps; physics::Symbol=Symbol("quasi-tem"))
    job = worker.job
    for basis in sort!(collect(worker.pending))
        if _adopt_column!(run, job.directory, job.frequency_index, job.frequency,
            basis, size(valid, 1), maps, job.mesh_digest; physics)
            valid[basis, job.frequency_index] = true
            delete!(worker.pending, basis)
            _record_progress!(run, valid)
        end
    end
end

function _finish_worker!(run, worker, execution; stopped = false)
    wait(worker.process)
    isopen(worker.log) && close(worker.log)
    elapsed = (time_ns() - worker.started_ns) / 1e9
    status = stopped ? "stopped" :
             success(worker.process) && isempty(worker.pending) ? "complete" : "failed"
    _write_json_atomic(joinpath(worker.job.directory, "attempt.json"),
        _attempt_record(worker.job; state = status, pid = worker.pid,
            process_token = worker.process_token,
            started_unix_seconds = worker.started, elapsed_seconds = elapsed,
            solver_threads = execution.solver_threads,
            exit_code = worker.process.exitcode, signal = worker.process.termsignal,
            completed_bases = setdiff(worker.job.bases, collect(worker.pending))))
    open(joinpath(run.path, "logs", "getdp.log"), "a") do io
        println(io,
            "frequency=$(worker.job.frequency_index) state=$status elapsed_seconds=$elapsed " *
            "exit_code=$(worker.process.exitcode) attempt=$(relpath(worker.job.directory, run.path))")
        if status != "complete"
            println(io, _log_tail(readlines(joinpath(worker.job.directory, "getdp.log"))))
        end
    end
    return elapsed
end

function _stop_workers!(run, active, valid, formulation, execution)
    errors = Any[]
    for worker in active
        try
            process_running(worker.process) && kill(worker.process)
        catch exception
            push!(errors, exception)
        end
    end
    for worker in active
        try
            if timedwait(() -> process_exited(worker.process), 2.0; pollint = 0.02) ===
               :timed_out
                kill(worker.process, 9)
            end
            wait(worker.process)
            _collect_worker_columns!(run, worker, valid, execution.plot_field_maps;
                physics=formulation.options.physics)
            _finish_worker!(run, worker, execution; stopped = true)
        catch exception
            push!(errors, exception)
        finally
            isopen(worker.log) && close(worker.log)
        end
    end
    empty!(active)
    for exception in errors
        @warn "Error while stopping a FEM worker" exception
    end
end

# A Linux PID can be reused after a crash. Include boot and process start time
# so a later unrelated process does not prevent recovery of this attempt.
function _process_token(pid::Integer)
    pid > 0 && Sys.islinux() || return nothing
    return try
        fields = split(last(split(read("/proc/$pid/stat", String), ") "; limit = 2)))
        strip(read("/proc/sys/kernel/random/boot_id", String)) * ":" * fields[20]
    catch
        nothing
    end
end

function _assert_no_live_attempts(run)
    root = joinpath(run.path, "attempts")
    isdir(root) || return nothing
    for directory in filter(isdir, readdir(root; join = true))
        record = try
            JSON3.read(read(joinpath(directory, "attempt.json"), String))
        catch
            continue
        end
        record isa AbstractDict || continue
        get(record, :state, "") == "running" || continue
        pid = get(record, :pid, 0)
        pid isa Integer && pid > 0 || continue
        status = ccall(:uv_kill, Cint, (Cint, Cint), pid, 0)
        # Query libuv's documented error name instead of Julia's private UV constant.
        alive = iszero(status) ||
                unsafe_string(ccall(:uv_err_name, Cstring, (Cint,), status)) != "ESRCH"
        token = get(record, :process_token, nothing)
        current = _process_token(Int(pid))
        alive && (token === nothing || current === nothing || token == current) || continue
        _fem_error(:execution, "GetDP", :ownership,
            "cannot resume while a previous GetDP process is still alive (PID $pid, attempt $directory)";
            run_directory = run.path)
    end
    return nothing
end

function _recover_columns!(run, model, maps, mesh_digests; physics::Symbol=Symbol("quasi-tem"))
    valid = falses(length(model.terminal_ids), length(model.problem.frequencies))
    for frequency in axes(valid, 2), basis in axes(valid, 1)

        valid[basis, frequency] = _valid_column_checkpoint(run.path, frequency,
            model.problem.frequencies[frequency], basis, size(valid, 1), maps, mesh_digests[frequency]; physics)
    end
    root = joinpath(run.path, "attempts")
    if isdir(root)
        for directory in sort!(filter(isdir, readdir(root; join = true)))
            manifest = joinpath(directory, "attempt.json")
            record = try
                JSON3.read(read(manifest, String))
            catch
                continue
            end
            record isa AbstractDict || continue
            frequency = get(record, :frequency_index, nothing)
            frequency isa Integer && frequency in axes(valid, 2) &&
            get(record, :protocol, 0) == 2 &&
            get(record, :mesh_digest, "") == mesh_digests[frequency] || continue
            for basis in axes(valid, 1)
                valid[basis, frequency] && continue
                valid[basis, frequency] = _adopt_column!(run, directory, frequency,
                    model.problem.frequencies[frequency], basis, size(valid, 1), maps, mesh_digests[frequency]; physics)
            end
        end
    end
    _record_progress!(run, valid)
    return valid
end

function _assemble_columns!(run, model)
    for quantity in (:Z, :P)
        destination = joinpath(run.path, "raw", "$quantity.tsv")
        temporary = tempname(dirname(destination))
        try
            open(temporary, "w") do io
                println(io, join(FEM_RAW_HEADER, '\t'))
                for frequency in eachindex(model.problem.frequencies),
                    basis in eachindex(model.terminal_ids)

                    paths = _job_raw_paths(run, _column_stem(frequency, basis))
                    write(io, read(getproperty(paths, quantity)))
                end
            end
            mv(temporary, destination; force = true)
        finally
            rm(temporary; force = true)
        end
    end
    _write_scan_completion!(run, model)
    return nothing
end

function _mesh_work_size(path)
    # Read only the ASCII MSH header, once per unique mesh, outside worker loops.
    # Unsupported external formats simply have no mesh-weighted estimate.
    try
        return open(path) do io
            strip(readline(io)) == "\$MeshFormat" || return nothing
            format = split(readline(io))
            length(format) == 3 && format[1] == "4.1" && format[2] == "0" || return nothing
            for line in eachline(io)
                strip(line) == "\$Nodes" || continue
                fields = split(readline(io))
                length(fields) == 4 || return nothing
                count = tryparse(Int, fields[2])
                return count !== nothing && count > 0 ? Float64(count) : nothing
            end
            return nothing
        end
    catch error
        error isa InterruptException && rethrow()
        return nothing
    end
end

function _fem_remaining_seconds(work, rates, active, next_job, workers, now_ns)
    length(rates) >= 3 || return -1.0
    # A median of recent seconds per mesh-work unit limits startup outliers.
    ordered = sort(rates)
    n = length(ordered)
    rate = (ordered[div(n+1,2)] + ordered[div(n+2,2)]) / 2
    # Simulate the existing FIFO scheduler, including its serial tail. Summing
    # worker wall durations would overestimate concurrent work.
    lanes = zeros(Float64, workers)
    for (index, (job_index, started_ns)) in enumerate(active)
        elapsed = max(0.0, Float64(now_ns-started_ns)*1e-9)
        predicted = rate*work[job_index]
        elapsed >= predicted && return -1.0 # Overdue work is still unfinished.
        lanes[index] = predicted-elapsed
    end
    for index in next_job:length(work)
        lane = argmin(lanes)
        lanes[lane] += rate*work[index]
    end
    return maximum(lanes)
end

function _run_getdp!(run::FEMRun, model::FEMResolvedModel, formulation::LineCableModelsFEM,
        execution::ComputationOptions, mesh_paths::AbstractVector{<:AbstractString}; pump = () -> true,
        reuse_factorization::Bool = true, batch_terminals::Bool = true)
    length(mesh_paths) == length(model.problem.frequencies) ||
        throw(DimensionMismatch("one FEM mesh path is required per frequency"))
    executable = _resolve_getdp(execution, run)
    assets = _getdp_assets(joinpath(run.path, "input", "getdp"))
    all(isfile, assets) || _fem_error(:getdp, "GetDP", :assets,
        "one or more retained solver sources are missing"; run_directory = run.path)
    mesh_digests = [bytes2hex(open(sha256, path)) for path in mesh_paths]
    valid = _recover_columns!(run, model, execution.plot_field_maps, mesh_digests;
        physics=formulation.options.physics)
    recovered_columns=count(valid)
    recovered_jobs=batch_terminals ? count(all,eachcol(valid)) : recovered_columns
    receiver=LineCableModels.progress_receiver()
    pending = Tuple{Int, Vector{Int}}[]
    for frequency in axes(valid, 2)
        bases = findall(!, valid[:, frequency])
        isempty(bases) && continue
        if batch_terminals
            push!(pending, (frequency, bases))
        else
            append!(pending, [(frequency, [basis]) for basis in bases])
        end
    end
    active = FEMActiveWorker[]
    next_job = 1
    worker_wall_seconds=0.0
    previous_progress=nothing
    # Node count is a proxy for sparse factorization work, scaled by remaining
    # terminal columns. Calibrate it against measured process wall durations.
    work = Float64[]
    work_index = Dict{Int, Int}()
    rates = Float64[]
    if receiver !== nothing
        sizes = Dict(path => _mesh_work_size(path) for path in unique(mesh_paths))
        if all(size -> size !== nothing, values(sizes))
            work = [sizes[mesh_paths[frequency]]^1.5 * length(bases)
                for (frequency, bases) in pending]
        end
    end
    # Mixed meshes or partially recovered terminal batches have unequal costs.
    # Their job counts remain useful, but do not imply equal-cost throughput.
    estimate_throughput=all(==(first(mesh_digests)),mesh_digests) &&
        (isempty(pending) || all(job->length(job[2])==length(first(pending)[2]),pending))
    receiver === nothing || LineCableModels.report_progress(receiver,
        (stage=:solving,unit=:jobs,completed=recovered_jobs,
            total=length(pending)+recovered_jobs,workers=0,queued=length(pending),
            recovered=recovered_jobs,estimate_throughput))
    try
        while next_job <= length(pending) || !isempty(active)
            pump() || _fem_error(:cancelled, "GetDP", :ui,
                "FEM solve cancelled; completed terminal columns are retained"; run_directory = run.path)
            while next_job <= length(pending) &&
                length(active) < execution.frequency_workers
                frequency, bases = pending[next_job]
                job = _frequency_job(
                    run, model, formulation, execution, executable, mesh_paths[frequency],
                    mesh_digests[frequency], frequency, bases; reuse_factorization)
                push!(active, _start_worker!(run, job, execution))
                isempty(work) || (work_index[active[end].pid] = next_job)
                next_job += 1
            end
            for index in reverse(eachindex(active))
                worker = active[index]
                _collect_worker_columns!(run, worker, valid, execution.plot_field_maps;
                    physics=formulation.options.physics)
                process_exited(worker.process) || continue
                # The last marker can arrive between the poll above and exit.
                _collect_worker_columns!(run, worker, valid, execution.plot_field_maps;
                    physics=formulation.options.physics)
                duration = _finish_worker!(run, worker, execution)
                worker_wall_seconds += duration
                if !isempty(work) && success(worker.process) && isempty(worker.pending)
                    push!(rates, duration/work[pop!(work_index, worker.pid)])
                    length(rates) > 32 && popfirst!(rates)
                end
                deleteat!(active, index)
                if !success(worker.process) || !isempty(worker.pending)
                    tail = _log_tail(readlines(joinpath(worker.job.directory, "getdp.log")))
                    capability = occursin(r"(?i)(unknown|syntax error).*?(GenerateRHS|SolveAgain)", tail)
                    requirement = capability ?
                                  "This backend requires GetDP with GenerateRHSGroup and SolveAgain support. " :
                                  ""
                    _fem_error(:getdp, "GetDP", capability ? :capability : :client,
                        requirement *
                        "GetDP frequency $(worker.job.frequency_index) failed " *
                        "(exit $(worker.process.exitcode)); missing columns $(sort!(collect(worker.pending))). " *
                        "Attempt: $(worker.job.directory)\nGetDP log tail:\n$tail"; run_directory = run.path)
                end
            end
            current=(next_job,length(active),run.completed_columns)
            if receiver !== nothing && current != previous_progress
                remaining_seconds = isempty(work) ? -1.0 : _fem_remaining_seconds(
                    work, rates, [(work_index[w.pid], w.started_ns) for w in active],
                    next_job, execution.frequency_workers, time_ns())
                LineCableModels.report_progress(receiver,
                    (stage=:solving,unit=:jobs,
                        completed=recovered_jobs+next_job-1-length(active),
                        total=length(pending)+recovered_jobs,workers=length(active),
                        queued=length(pending)-next_job+1,recovered=recovered_jobs,
                        remaining_seconds,eta_source=:mesh_work))
                previous_progress=current
            end
            isempty(active) || sleep(0.02)
        end
    catch exception
        _stop_workers!(run, active, valid, formulation, execution)
        if exception isa InterruptException ||
           exception isa LineCableModelsFEMError && exception.category === :cancelled
            _transition!(run, cancelled, "solver processes stopped; completed columns retained")
        end
        rethrow()
    finally
        isempty(active) || _stop_workers!(run, active, valid, formulation, execution)
    end
    all(valid) || _fem_error(:getdp, "GetDP", :raw_output,
        "frequency scan returned without all terminal columns"; run_directory = run.path)
    _assemble_columns!(run, model)
    timings=[_column_timing(_column_paths(run.path,frequency,basis,
            execution.plot_field_maps).timing,frequency,basis)
        for frequency in axes(valid,2) for basis in axes(valid,1)]
    totals=map((:constraint_seconds,:assembly_seconds,:solve_seconds,:output_seconds)) do key
        sum(getproperty(timing,key) for timing in timings)
    end
    _write_json_atomic(joinpath(run.path,"timing-summary.json"),
        (;schema=1,backend="getdp",scope="accumulated native worker wall time; not elapsed scan time",
            constraint_seconds=totals[1],assembly_seconds=totals[2],solve_seconds=totals[3],
            output_seconds=totals[4],columns=length(timings),recovered_columns,
            worker_wall_seconds,
            worker_wall_scope="sum of newly executed process wall durations in this invocation",
            factorized_columns=count(timing->timing.factorized,timings)))
    return nothing
end

function _run_getdp!(run::FEMRun, model::FEMResolvedModel, formulation::LineCableModelsFEM,
        execution::ComputationOptions, mesh_path::String; kwargs...)
    length(model.problem.frequencies) == 1 || throw(DimensionMismatch(
        "a single FEM mesh path is valid only for a one-frequency problem"))
    return _run_getdp!(run, model, formulation, execution, [mesh_path]; kwargs...)
end
