"""
    UIHostHandle

Track one owned operating-system process. Credentials and process launch
environment are intentionally absent from display and diagnostics.
"""
mutable struct UIHostHandle
    "Durable run identity."
    id::UUID
    "Original authenticated owner."
    owner::Principal
    "Owned child process; never reconstructed from a stored PID."
    process::Base.Process
    "Private directory containing this child's readiness receipt."
    directory::String
    "Per-launch proxy credential."
    key::String
    "Verified loopback listener port, or nothing before readiness."
    port::Union{Nothing,Int}
    "Whether an explicit stop was requested."
    stopping::Bool
    "Whether terminal bookkeeping and cleanup completed."
    finished::Bool
    "Owned lifecycle monitor."
    monitor::Union{Nothing,Task}
    "Number of currently authorized browser WebSockets."
    connections::Int
    "Monotonic time of the last browser activity or completed startup."
    last_seen::UInt64
end
Base.show(io::IO, handle::UIHostHandle) = print(io, "UIHostHandle(", handle.id, ", <private>)")
Base.show(io::IO, ::MIME"text/plain", handle::UIHostHandle) = show(io, handle)

"""
    UIHostSupervisor(store, registry, scratch_root; limits=RunLimits())

Supervise approved local UI processes without loading their UI implementations.
An exclusive Linux file lock permits only one supervisor per database. Existing
active records become failed on recovery; they never imply restored UI memory.
"""
mutable struct UIHostSupervisor
    "Durable run bookkeeping."
    store::RuntimeStore
    "Trusted application declarations and launch hooks."
    registry::ApplicationRegistry
    "Private directory reserved for owned UI-host records."
    root::String
    "Admission and lifetime bounds."
    limits::RunLimits
    "Live process handles owned by this supervisor."
    handles::Dict{UUID,UIHostHandle}
    "Serialize launch, completion and shutdown state."
    lock::ReentrantLock
    "Kernel-held single-supervisor lock."
    ownership_lock::IOStream
    "Whether new launches are forbidden."
    closed::Bool
end

function UIHostSupervisor(store::RuntimeStore, registry::ApplicationRegistry,
        scratch_root::AbstractString; limits::RunLimits=RunLimits())
    Sys.islinux() || throw(ArgumentError("v1 UI supervision requires Linux"))
    isnothing(Sys.which("setpriv")) &&
        throw(ArgumentError("setpriv is required for parent-death cleanup"))
    root = String(rstrip(abspath(scratch_root), '/'))
    isempty(root) && (root = "/")
    root in ("/", homedir(), pwd(), dirname(store.path)) &&
        throw(ArgumentError("UI hosts require a dedicated scratch directory"))
    islink(root) && throw(ArgumentError("UI scratch root cannot be a symbolic link"))
    isdir(root) || mkpath(root; mode=0o700)
    lockpath = store.path * ".supervisor.lock"
    islink(lockpath) && throw(ArgumentError("supervisor lock cannot be a symbolic link"))
    ownership = open(lockpath, "a+")
    # flock(LOCK_EX | LOCK_NB): released by the kernel even on SIGKILL. No
    # age-based lock stealing may create two live allocation authorities.
    if ccall(:flock, Cint, (Cint, Cint), fd(ownership), 6) != 0
        close(ownership)
        throw(ArgumentError("a supervisor already owns this database"))
    end
    supervisor = UIHostSupervisor(store, registry, realpath(root), limits,
        Dict{UUID,UIHostHandle}(), ReentrantLock(), ownership, false)
    try
        recovery = Principal("runtime-recovery"; administrator=true)
        for run in list_runs(store, recovery)
            if run.state in ACTIVE_RUN_STATES
                transition_run!(store, recovery, run.id, :failed;
                    reason="UI host lost during coordinator restart")
            end
            try_remove_ui_directory!(supervisor, run.id)
        end
        return supervisor
    catch
        close(ownership)
        rethrow()
    end
end

function remove_ui_directory!(supervisor::UIHostSupervisor, id::UUID)
    directory = joinpath(supervisor.root, string(id))
    ispath(directory) || return
    islink(directory) && throw(ArgumentError("refusing a linked UI-host directory"))
    marker = joinpath(directory, "ownership.json")
    isfile(marker) && !islink(marker) && filesize(marker) <= 1024 ||
        throw(ArgumentError("UI-host ownership receipt missing"))
    receipt = JSON3.read(read(marker, String))
    receipt.kind == "lcm-ui-host-v1" && receipt.run_id == string(id) ||
        throw(ArgumentError("UI-host ownership receipt does not match"))
    rm(directory; recursive=true)
end

function try_remove_ui_directory!(supervisor::UIHostSupervisor, id::UUID)
    try
        remove_ui_directory!(supervisor, id)
        return true
    catch
        # A missing/tampered receipt never permits broader deletion and must
        # not prevent cleanup of other owned processes. Do not log paths.
        @warn "UI scratch cleanup requires operator attention" run_id=id
        return false
    end
end

function child_environment(context::HostContext, key::String)
    # Do not pass publisher/proxy/broker/storage credentials to UI code.
    environment = Dict(name => ENV[name] for name in
        ("PATH", "HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TMPDIR",
         "JULIA_DEPOT_PATH", "JULIA_CPU_TARGET") if haskey(ENV, name))
    merge!(environment, Dict(
        "JULIA_NUM_THREADS" => "2",
        "LCM_RUN_ID" => string(context.run_id),
        "LCM_RUN_PREFIX" => context.prefix,
        "LCM_UI_READY_FILE" => context.ready_file,
        "LCM_UI_HOST_KEY" => key,
        "LCM_UI_PARENT_PID" => string(getpid()),
        "TMPDIR" => dirname(context.ready_file),
    ))
    return environment
end

function finish_ui!(supervisor::UIHostSupervisor, handle::UIHostHandle, reason::String)
    lock(supervisor.lock) do
        handle.finished && return
        Base.process_running(handle.process) &&
            throw(ArgumentError("cannot release a live UI process"))
        state = handle.stopping ? :stopped : :failed
        cleaned = try_remove_ui_directory!(supervisor, handle.id)
        final_reason = cleaned ? reason : reason * "; owned scratch cleanup requires operator attention"
        transition_run!(supervisor.store, handle.owner, handle.id, state; reason=final_reason)
        handle.finished = true
        delete!(supervisor.handles, handle.id)
    end
end

function terminate_ui_process!(handle::UIHostHandle, grace::Real)
    if Base.process_running(handle.process)
        try
            kill(handle.process, Base.SIGTERM)
        catch error
            error isa Base.IOError || rethrow()
        end
        timedwait(() -> !Base.process_running(handle.process), grace; pollint=0.02)
        Base.process_running(handle.process) && kill(handle.process, Base.SIGKILL)
    end
    wait(handle.process)
end

function ui_readiness(handle::UIHostHandle)
    file = joinpath(handle.directory, "ready.json")
    isfile(file) && !islink(file) && filesize(file) <= 4096 || return nothing
    payload = try
        JSON3.read(read(file, String))
    catch
        return nothing # An atomic receipt is preferred; tolerate an incomplete write.
    end
    get(payload, :schema_version, nothing) === 1 &&
        get(payload, :run_id, "") == string(handle.id) &&
        get(payload, :pid, 0) == getpid(handle.process) || return nothing
    port = get(payload, :port, 0)
    port isa Integer && !(port isa Bool) && 1 <= port <= 65535 || return nothing
    response = try
        HTTP.get("http://127.0.0.1:$port/health";
            headers=["X-LCM-Host-Key"=>handle.key], redirect=false, proxy=nothing,
            retry=false, status_exception=false, request_timeout=1)
    catch
        return nothing
    end
    response.status == 200 && String(response.body) == string(handle.id) || return nothing
    return Int(port)
end

function monitor_ui!(supervisor::UIHostSupervisor, handle::UIHostHandle)
    started = time_ns()
    try
        while Base.process_running(handle.process)
            if handle.port === nothing && !handle.stopping
                port = ui_readiness(handle)
                if port !== nothing
                    lock(supervisor.lock) do
                        if !handle.stopping && !handle.finished
                            handle.port = port
                            handle.last_seen = time_ns()
                            transition_run!(supervisor.store, handle.owner, handle.id, :running)
                        end
                    end
                elseif (time_ns() - started) / 1e9 > supervisor.limits.startup_seconds
                    terminate_ui_process!(handle, supervisor.limits.shutdown_seconds)
                    finish_ui!(supervisor, handle, "UI host startup timed out")
                    return
                end
            end
            expired = lock(supervisor.lock) do
                if handle.port !== nothing && !handle.stopping && handle.connections == 0 &&
                        (time_ns() - handle.last_seen) / 1e9 > supervisor.limits.disconnect_grace_seconds
                    handle.stopping = true
                    transition_run!(supervisor.store, handle.owner, handle.id, :stopping;
                        reason="Disconnected application grace expired")
                    return true
                end
                return false
            end
            if expired
                terminate_ui_process!(handle, supervisor.limits.shutdown_seconds)
                finish_ui!(supervisor, handle, "Disconnected application grace expired; volatile state was lost")
                return
            end
            sleep(0.05)
        end
        finish_ui!(supervisor, handle,
            handle.stopping ? "UI host stopped" : "UI host exited; volatile state was lost")
    catch
        # Never expose a Cmd/environment or raw child error through diagnostics.
        terminate_ui_process!(handle, supervisor.limits.shutdown_seconds)
        finish_ui!(supervisor, handle, "UI host supervision failed")
    end
end

"""
    start_ui!(supervisor, principal, application_id; request_id=uuid4()) -> RunRecord

Reserve capacity and launch one approved process. Return without waiting for UI
startup; inspect its durable state for completion. Repeated requests return the
same run and never start a second process.
"""
function start_ui!(supervisor::UIHostSupervisor, principal::Principal,
        application_id::AbstractString; request_id::UUID=uuid4())
    lock(supervisor.lock) do
        supervisor.closed && throw(ArgumentError("supervisor is closed"))
        haskey(supervisor.registry.definitions, application_id) ||
            throw(AccessDenied(404, "Application not found"))
        haskey(supervisor.registry.applications, application_id) ||
            throw(AccessDenied(409, "Application live implementation is not installed"))
        description = supervisor.registry.definitions[application_id]
        run = reserve_run!(supervisor.store, principal, description;
            limits=supervisor.limits, request_id)
        run.state == :reserved || return run
        directory = joinpath(supervisor.root, string(run.id))
        context = HostContext(run.id, "/applications/runs/$(run.id)/",
            joinpath(directory, "ready.json"))
        transition_run!(supervisor.store, principal, run.id, :starting)
        created = false
        process = nothing
        try
            mkdir(directory; mode=0o700)
            created = true
            write(joinpath(directory, "ownership.json"),
                JSON3.write((kind="lcm-ui-host-v1", run_id=string(run.id))))
            key = bytes2hex(rand(Random.RandomDevice(), UInt8, 32))
            command = Base.invokelatest(ui_command,
                supervisor.registry.applications[application_id], context)
            command isa Cmd || throw(ArgumentError("ui_command must return Cmd"))
            launcher = Sys.which("setpriv")
            guarded = `$launcher --pdeathsig KILL -- $command`
            process = Base.run(pipeline(ignorestatus(setenv(guarded,
                child_environment(context, key))); stdout=devnull, stderr=devnull); wait=false)
            handle = UIHostHandle(run.id, principal, process, directory, key,
                nothing, false, false, nothing, 0, time_ns())
            supervisor.handles[run.id] = handle
            handle.monitor = @async monitor_ui!(supervisor, handle)
            return get_run(supervisor.store, principal, run.id)
        catch
            if process !== nothing
                Base.process_running(process) && kill(process, Base.SIGKILL)
                wait(process)
            end
            # Only this invocation's fresh UUID directory is eligible here.
            # It may have failed before an ownership receipt could be written.
            cleanup_failed = false
            if created
                try
                    islink(directory) && throw(ArgumentError("linked launch directory"))
                    rm(directory; recursive=true)
                catch
                    cleanup_failed = true
                end
            end
            transition_run!(supervisor.store, principal, run.id, :failed;
                reason=cleanup_failed ? "UI launch failed; owned scratch cleanup requires operator attention" :
                    "UI host could not be launched")
            throw(ArgumentError("UI host could not be launched"))
        end
    end
end

"""
    stop_ui!(supervisor, principal, id) -> RunRecord

Stop an owned UI process within the configured grace, then forcibly terminate
it if necessary. Release only the owned process and receipt-validated directory.
"""
function stop_ui!(supervisor::UIHostSupervisor, principal::Principal, id::UUID)
    handle = lock(supervisor.lock) do
        run = get_run(supervisor.store, principal, id)
        run.state in (:stopped, :failed) && return nothing
        handle = get(supervisor.handles, id, nothing)
        isnothing(handle) && throw(ArgumentError("active UI process is not owned by this supervisor"))
        handle.stopping = true
        transition_run!(supervisor.store, principal, id, :stopping)
        return handle
    end
    if handle !== nothing
        terminate_ui_process!(handle, supervisor.limits.shutdown_seconds)
        finish_ui!(supervisor, handle, "UI host stopped")
    end
    return get_run(supervisor.store, principal, id)
end

function Base.close(supervisor::UIHostSupervisor)
    handles = lock(supervisor.lock) do
        supervisor.closed && return nothing
        supervisor.closed = true
        return collect(values(supervisor.handles))
    end
    handles === nothing && return
    try
        # All owned processes receive stop promptly. One cleanup failure must
        # neither skip the remaining handles nor multiply the grace by count.
        @sync for handle in handles
            @async stop_ui!(supervisor, handle.owner, handle.id)
        end
    finally
        close(supervisor.ownership_lock)
    end
end
