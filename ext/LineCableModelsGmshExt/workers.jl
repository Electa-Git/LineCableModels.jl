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
        physics::Symbol=:helmholtz)
    stem = _column_stem(frequency, basis)
    raw = _job_raw_paths(root, stem)
    return (; raw...,
        diagnostics = [joinpath(root, "results", "jobs", "$stem-Pscalar.tsv")],
        timing = joinpath(root, "results", "jobs", "$stem-timing.tsv"),
        marker = joinpath(root, "results", "jobs", "$stem.done"),
        checkpoint = joinpath(root, "results", "jobs", "$stem.json"),
        pml = joinpath(root, "results", "jobs", @sprintf("pml-f%04d.tsv", frequency)),
        maps = maps ?
               [joinpath(root, "maps",
                    @sprintf("%s_f%04d_b%04d.pos", quantity, frequency, basis))
                for quantity in _field_quantities(physics)] : String[])
end

_column_files(paths) = [paths.Z, paths.P, paths.timing, paths.pml, paths.maps..., paths.diagnostics...]

function _pml_observation(path, frequency_index, frequency)
    isfile(path) || return nothing
    rows = readlines(path)
    length(rows) == 1 || return nothing
    values = tryparse.(Float64, split(only(rows), '\t'))
    length(values) == 24 && all(v -> v !== nothing && isfinite(v), values) || return nothing
    values[1] == frequency_index && isapprox(values[2], frequency; rtol=16eps(Float64), atol=0) || return nothing
    values[3] > 0 && all(v -> v in (0,1), values[[4,5,6,7,12,19,21,22,23,24]]) || return nothing
    all(v -> 1 <= v <= typemax(Int) && isinteger(v), values[13:15]) || return nothing
    all(v -> 0 <= v <= 1, values[16:18]) && values[20] >= 0 || return nothing
    return (frequency_index=frequency_index, frequency_hz=values[2],
        target_exponent=values[3], sizing_floor_flags=(Bool(values[4]),Bool(values[5])),
        cutoff_flags=(Bool(values[6]),Bool(values[7])),
        net_side_exponents=(values[8],values[9]), net_top_exponent=values[10],
        net_bottom_exponent=values[11], sizing_floor_active=Bool(values[12]),
        earth_sizing_ceiling_active=Bool(values[19]),
        earth_layer_thickness_m=values[20],
        earth_layer_clipped_or_omitted=Bool(values[21]),
        pml_eta=(values[16],values[17],values[18]),
        effective_pml_layers=(Int(values[13]),Int(values[14]),Int(values[15])),
        quasi_static_cap_flags=(Bool(values[22]),Bool(values[23]),Bool(values[24])))
end

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
        physics::Symbol=:helmholtz)
    paths = _column_paths(root, frequency_index, basis, maps; physics)
    isfile(paths.checkpoint) || return false
    return try
        record = JSON3.read(read(paths.checkpoint, String))
        record.protocol == 2 && record.frequency_index == frequency_index &&
        record.basis == basis && record.terminals == terminals &&
        get(record, :physics, nothing) == String(physics) &&
        record.plot_field_maps == maps && record.mesh_digest == mesh_digest &&
        isapprox(record.frequency_hz, frequency; rtol = 16eps(Float64), atol = 0) ||
            return false
        _valid_job_raw(paths.Z, terminals, frequency_index, frequency, basis) &&
        _valid_job_raw(paths.P, terminals, frequency_index, frequency, basis) ||
            return false
        _column_timing(paths.timing, frequency_index, basis) === nothing && return false
        _pml_observation(paths.pml, frequency_index, frequency) === nothing && return false
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
        physics::Symbol=:helmholtz)
    paths = _column_paths(source, frequency_index, basis, maps; physics)
    _valid_column_marker(
        paths.marker, frequency_index, frequency, basis, terminals, maps) || return false
    _valid_job_raw(paths.Z, terminals, frequency_index, frequency, basis) &&
    _valid_job_raw(paths.P, terminals, frequency_index, frequency, basis) || return false
    all(path -> _valid_job_raw(path, terminals, frequency_index, frequency, basis),
        paths.diagnostics) || return false
    timing = _column_timing(paths.timing, frequency_index, basis)
    timing === nothing && return false
    _pml_observation(paths.pml, frequency_index, frequency) === nothing && return false
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

function _getdp_command(executable, model_path, mesh_path, run, formulation, execution, frequency_index,
        bases, directory; reuse_factorization = true)
    basis_path = joinpath(directory, "bases.pro")
    write(basis_path, "RequestedBases() = $(_pro_array(bases));\n")
    verbosity = execution.data.getdp_verbosity
    arguments = [executable, model_path, "-solve", "LineCableModelsFEMScan",
        "-setnumber", "Physics", string(_fem_physics_code(formulation)),
        "-msh", abspath(mesh_path), "-name", joinpath(directory, "solver"),
        "-v", string(verbosity),
        "-setstring", "ModelDataPath", joinpath(run.path, "input", "model_data.pro"),
        "-setstring", "RunDirectory", directory,
        "-setstring", "RawDirectory", joinpath(directory,"results"),
        "-setstring", "BasisListPath", basis_path,
        "-setnumber", "FrequencyIndex", string(frequency_index),
        "-setnumber", "PlotFieldMaps", string(Int(execution.data.plot_field_maps)),
        "-setnumber", "ReuseFactorization", string(Int(reuse_factorization)),
        "-setnumber", "LinearSolver", string(Int(execution.data.linear_solver === :gmres)),
        "-setnumber", "GmresIterationsMax", string(execution.data.gmres_iterations_max),
        "-setnumber", "GmresRelativeTolerance", _pro_number(execution.data.gmres_relative_tolerance),
        "-setnumber", "GmresAbsoluteTolerance", _pro_number(execution.data.gmres_absolute_tolerance),
        "-setnumber", "MumpsOrdering", string(execution.data.mumps_ordering),
        "-setnumber", "PetscPrealloc", string(execution.data.petsc_prealloc),
        "-setnumber", "MumpsErrorAnalysis", string(execution.data.mumps_error_analysis),
        "-setnumber", "MumpsRefinementMax", string(execution.data.mumps_refinement_max),
        "-setnumber", "MumpsBackwardErrorTolerance", _pro_number(execution.data.mumps_backward_error_tolerance)]
    if verbosity >= 4
        append!(arguments, [
            "-cpu", "-ksp_view", "-log_view", ":" * joinpath(directory, "petsc.log")])
    end
    threads = string(execution.data.solver_threads)
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
    root = joinpath(run.path, "work")
    mkpath(root)
    directory = joinpath(root,@sprintf("frequency_%04d",frequency_index))
    if ispath(directory) || ispath(joinpath(run.path,"logs",basename(directory)))
        directory = mktempdir(root;prefix=@sprintf("frequency_%04d-",frequency_index),cleanup=false)
    else
        mkdir(directory)
    end
    mkpath(joinpath(directory,"results","jobs"))
    execution.data.plot_field_maps && mkpath(joinpath(directory,"maps"))
    command = _getdp_command(
        executable, _getdp_assets(joinpath(run.path, "input", "getdp")).model,
        mesh, run, formulation, execution, frequency_index, bases, directory; reuse_factorization)
    return FEMFrequencyJob(frequency_index, Float64(model.problem.frequencies[frequency_index]), collect(bases),
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
                solver_threads = execution.data.solver_threads))
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

function _collect_worker_columns!(run, worker, valid, maps; physics::Symbol=:helmholtz)
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

# Read PETSc's native post-solve MUMPS view, not GetDP's original-coordinate
# residual. Missing/changed diagnostic output never rejects a computed column.
function _solver_diagnostics(text::AbstractString, ::Val{:mumps})
    return map(eachmatch(r"(?ms)^KSP Object:.*?(?=^KSP Object:|\z)", text)) do block
        number(pattern) = begin
            found = match(pattern, block.match)
            found === nothing ? nothing : tryparse(Float64, found[1])
        end
        analysis = number(r"ICNTL\(11\)\s*\(error analysis\):\s*(\S+)")
        backward = match(r"RINFOG\(7\),\s*RINFOG\(8\)\s*\(backward error est\):\s*([^,\s]+),\s*(\S+)", block.match)
        omega1 = backward === nothing ? nothing : tryparse(Float64, backward[1])
        omega2 = backward === nothing ? nothing : tryparse(Float64, backward[2])
        conditions = analysis == 1 ? match(
            r"RINFOG\(10\),\s*RINFOG\(11\)\s*\(condition numbers\):\s*([^,\s]+),\s*(\S+)", block.match) : nothing
        return (linear_solver=:mumps, error_analysis=analysis,
            refinement_steps=number(r"INFOG\(15\)\s*\(number of steps of iterative refinement after solution\):\s*(\S+)"),
            omega1, omega2,
            backward_error=omega1 === nothing || omega2 === nothing ? nothing : omega1 + omega2,
            scaled_residual=number(r"RINFOG\(6\)\s*\(inf norm of residual\):\s*(\S+)"),
            forward_error=analysis == 1 ? number(r"RINFOG\(9\)\s*\(error estimate\):\s*(\S+)") : nothing,
            cond1=conditions === nothing ? nothing : tryparse(Float64, conditions[1]),
            cond2=conditions === nothing ? nothing : tryparse(Float64, conditions[2]))
    end
end

function _solver_diagnostics(text::AbstractString, ::Val{:gmres})
    # PETSc's iteration-zero monitor begins each RHS even at GetDP verbosity 0.
    # GetDP's original-coordinate residual is optional at higher verbosity.
    # MUMPS RINFOG values here would describe a preconditioner application.
    return map(eachmatch(r"(?ms)^[ \t]*0 KSP unpreconditioned resid norm.*?(?=^[ \t]*0 KSP unpreconditioned resid norm|\z)", text)) do block
        reason = match(r"Linear solve (?:did not converge|converged) due to (\S+) iterations (\d+)", block.match)
        original = match(r"FEM algebraic residual: frequency=\d+ basis=\d+ residual=(\S+) rhs=(\S+)", block.match)
        monitors = collect(eachmatch(r"KSP unpreconditioned resid norm\s+(\S+) true resid norm\s+(\S+) \|\|r\(i\)\|\|/\|\|b\|\|\s+(\S+)", block.match))
        last_monitor = isempty(monitors) ? nothing : last(monitors)
        return (linear_solver=:gmres,
            iterations=reason === nothing ? nothing : tryparse(Int, reason[2]),
            convergence_reason=reason === nothing ? nothing : String(reason[1]),
            estimated_residual_norm=last_monitor === nothing ? nothing : tryparse(Float64, last_monitor[1]),
            scaled_residual_norm=last_monitor === nothing ? nothing : tryparse(Float64, last_monitor[2]),
            scaled_relative_residual=last_monitor === nothing ? nothing : tryparse(Float64, last_monitor[3]),
            original_residual_norm=original === nothing ? nothing : tryparse(Float64, original[1]),
            original_rhs_norm=original === nothing ? nothing : tryparse(Float64, original[2]))
    end
end

function _warn_solver_diagnostics(diagnostic, controls, ::Val{:mumps}; frequency_hz, basis, log)
    if diagnostic.error_analysis != controls.mumps_error_analysis ||
            diagnostic.backward_error === nothing || diagnostic.refinement_steps === nothing
        @warn "MUMPS diagnostics unavailable or inconsistent with requested analysis; inspect native log" frequency_hz basis log
        return nothing
    end
    error = diagnostic.backward_error
    target = controls.mumps_backward_error_tolerance
    steps = diagnostic.refinement_steps
    if !isfinite(error) || error > target
        @warn "MUMPS backward-error target not met; result retained" frequency_hz basis backward_error=error target refinement_steps=steps log
    end
    if controls.mumps_error_analysis == 1
        estimate = diagnostic.forward_error
        target = FEM_MUMPS_FORWARD_ERROR_BUDGET
        if estimate === nothing
            @warn "MUMPS forward-error estimate unavailable; inspect native log" frequency_hz basis log
        elseif !isfinite(estimate) || estimate > target
            @warn "MUMPS estimated scaled-solution sensitivity exceeds budget; not an error estimate for terminal G/Y; result retained" frequency_hz basis forward_error=estimate target refinement_steps=steps log
        end
    end
    return nothing
end

function _warn_solver_diagnostics(diagnostic, controls, ::Val{:gmres}; frequency_hz, basis, log)
    reason, iterations = diagnostic.convergence_reason, diagnostic.iterations
    residual, relative = diagnostic.scaled_residual_norm, diagnostic.scaled_relative_residual
    if reason === nothing || iterations === nothing || residual === nothing || relative === nothing
        @warn "GMRES diagnostics unavailable; inspect native log" frequency_hz basis log
        return nothing
    end
    rtol, atol = controls.gmres_relative_tolerance, controls.gmres_absolute_tolerance
    target_met = isfinite(residual) &&
        (residual <= atol || isfinite(relative) && relative <= rtol)
    if !startswith(reason, "CONVERGED_") || !target_met
        @warn "GMRES convergence or recomputed residual target not met; result retained" frequency_hz basis iterations reason scaled_residual_norm=residual scaled_relative_residual=relative relative_tolerance=rtol absolute_tolerance=atol log
    end
    return nothing
end

function _finish_worker!(run, worker, execution; stopped = false)
    wait(worker.process)
    isopen(worker.log) && close(worker.log)
    elapsed = (time_ns() - worker.started_ns) / 1e9
    status = stopped ? "stopped" :
             success(worker.process) && isempty(worker.pending) ? "complete" : "failed"
    if !stopped && (execution.data.linear_solver === :gmres || execution.data.mumps_error_analysis != 0)
        solver = Val(execution.data.linear_solver)
        log = joinpath(worker.job.directory, "getdp.log")
        records = _solver_diagnostics(read(log, String), solver)
        if length(records) == length(worker.job.bases)
            for (basis, diagnostic) in zip(worker.job.bases, records)
                frequency_hz = worker.job.frequency
                _warn_solver_diagnostics(diagnostic, execution.data, solver; frequency_hz, basis, log)
            end
        else
            @warn "FEM solver diagnostics incomplete; inspect native log" linear_solver=execution.data.linear_solver frequency_hz=worker.job.frequency expected=length(worker.job.bases) reported=length(records) log
        end
    end
    _write_json_atomic(joinpath(worker.job.directory, "attempt.json"),
        _attempt_record(worker.job; state = status, pid = worker.pid,
            process_token = worker.process_token,
            started_unix_seconds = worker.started, elapsed_seconds = elapsed,
            solver_threads = execution.data.solver_threads,
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
            _collect_worker_columns!(run, worker, valid, execution.data.plot_field_maps;
                physics=formulation.options.data.physics)
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
    root = joinpath(run.path, "work")
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

function _recover_columns!(run, model, maps, mesh_digests; physics::Symbol=:helmholtz)
    valid = falses(length(model.terminal_ids), length(model.problem.frequencies))
    for frequency in axes(valid, 2), basis in axes(valid, 1)

        valid[basis, frequency] = _valid_column_checkpoint(run.path, frequency,
            model.problem.frequencies[frequency], basis, size(valid, 1), maps, mesh_digests[frequency]; physics)
    end
    root = joinpath(run.path, "work")
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
        destination = joinpath(run.path, "results", "$quantity.tsv")
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

# Keep diagnostic evidence after adopted solver scratch is removed.
function _retain_worker_logs!(run, directory)
    destination=joinpath(run.path,"logs",basename(directory))
    for name in ("attempt.json","getdp.log","petsc.log")
        source=joinpath(directory,name)
        isfile(source) && _copy_column_file(source,joinpath(destination,name))
    end
    previous=relpath(directory,run.path)
    manifest=try
        JSON3.read(read(joinpath(directory,"attempt.json"),String))
    catch
        nothing
    end
    if manifest isa AbstractDict && get(manifest,:frequency_index,nothing) isa Integer
        frequency=manifest.frequency_index
        for basis in get(manifest,:requested_bases,Int[])
            path=_column_paths(run.path,frequency,Int(basis),false).checkpoint
            isfile(path) || continue
            record=JSON3.read(read(path,String),Dict{String,Any})
            get(record,"attempt",nothing)==previous || continue
            record["attempt"]=relpath(destination,run.path)
            _write_json_atomic(path,record)
        end
    end
    rm(directory;recursive=true,force=true)
    return nothing
end

function _run_getdp!(run::FEMRun, model::FEMResolvedModel, formulation::LineCableModelsFEM,
        execution::ComputationOptions, mesh_paths::AbstractVector{<:AbstractString};
        reuse_factorization::Bool = true, batch_terminals::Bool = true)
    length(mesh_paths) == length(model.problem.frequencies) ||
        throw(DimensionMismatch("one FEM mesh path is required per frequency"))
    executable = _resolve_getdp(execution, run)
    assets = _getdp_assets(joinpath(run.path, "input", "getdp"))
    all(isfile, assets) || _fem_error(:getdp, "GetDP", :assets,
        "one or more retained solver sources are missing"; run_directory = run.path)
    mesh_digests = [bytes2hex(open(sha256, path)) for path in mesh_paths]
    valid = _recover_columns!(run, model, execution.data.plot_field_maps, mesh_digests;
        physics=formulation.options.data.physics)
    recovered_columns=count(valid)
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
    progress = LineCableModels.verbosity(execution, :progress) > 0
    started = progress ? time_ns() : UInt64(0)
    previous = started
    last_log = started
    previous_columns = recovered_columns
    average_seconds = 0.0
    try
        while next_job <= length(pending) || !isempty(active)
            while next_job <= length(pending) &&
                length(active) < execution.data.frequency_workers
                frequency, bases = pending[next_job]
                job = _frequency_job(
                    run, model, formulation, execution, executable, mesh_paths[frequency],
                    mesh_digests[frequency], frequency, bases; reuse_factorization)
                push!(active, _start_worker!(run, job, execution))
                next_job += 1
            end
            for index in reverse(eachindex(active))
                worker = active[index]
                _collect_worker_columns!(run, worker, valid, execution.data.plot_field_maps;
                    physics=formulation.options.data.physics)
                process_exited(worker.process) || continue
                # The last marker can arrive between the poll above and exit.
                _collect_worker_columns!(run, worker, valid, execution.data.plot_field_maps;
                    physics=formulation.options.data.physics)
                duration = _finish_worker!(run, worker, execution)
                worker_wall_seconds += duration
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
                _retain_worker_logs!(run,worker.job.directory)
            end
            if progress && run.completed_columns > previous_columns
                now = time_ns()
                interval = (now - previous) * 1e-9 / (run.completed_columns - previous_columns)
                average_seconds = previous_columns == recovered_columns ? interval :
                    0.2 * interval + 0.8 * average_seconds
                previous = now
                previous_columns = run.completed_columns
                if now - last_log >= 5_000_000_000
                    @info "FEM frequency sweep progress" _group=:progress completed_columns=run.completed_columns total_columns=length(valid) completed_frequencies=run.completed_frequencies total_frequencies=size(valid,2) elapsed_seconds=(now-started)*1e-9 eta_hours=(length(valid)-run.completed_columns)*average_seconds/3600
                    last_log = now
                end
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
            execution.data.plot_field_maps).timing,frequency,basis)
        for frequency in axes(valid,2) for basis in axes(valid,1)]
    totals=map((:constraint_seconds,:assembly_seconds,:solve_seconds,:output_seconds)) do key
        sum(getproperty(timing,key) for timing in timings)
    end
    _write_json_atomic(joinpath(run.path,"results","timing-summary.json"),
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
