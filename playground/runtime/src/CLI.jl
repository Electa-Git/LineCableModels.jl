"""
    serve_runtime(config; xray=false)

Start an explicitly enabled private runtime gateway with its approved catalogue
and existing static build. Supervised stdin shutdown releases only resources
owned by this coordinator. Direct process death is covered by child guards.
"""
function serve_runtime(config::RuntimeConfig; xray::Bool=false)
    config.enabled || throw(ArgumentError("runtime is disabled; opt in with enabled = true"))
    sitepath = something(config.site_directory, joinpath(@__DIR__, "..", "..", "_site"))
    site = PublishedSite(sitepath)
    registry = default_applications(; xray)
    store = RuntimeStore(config.database)
    supervisor = nothing
    control = nothing
    server = nothing
    closed = Ref(false)
    shutdown_lock = ReentrantLock()
    function shutdown()
        lock(shutdown_lock) do
            closed[] && return
            # Closing the listener wakes the main task. It must join this
            # complete teardown before returning, not mistake "closing" for
            # "closed" and exit while child/scratch cleanup is still yielding.
            try
                server === nothing || HTTP.forceclose(server)
            finally
                try
                    try
                        control === nothing || close(control)
                    finally
                        supervisor === nothing || close(supervisor)
                    end
                finally
                    close(store)
                    closed[] = true
                end
            end
        end
    end
    try
        supervisor = UIHostSupervisor(store, registry, config.scratch_root; limits=config.limits)
        if config.control !== nothing
            control = ControlService(config.control, store, registry)
            # Complete gateway compilation before the broker scheduler starts;
            # it must not consume a newly connected worker's liveness window.
            compile_gateway_paths(supervisor, config.identity, site, control)
            start_control!(control)
        end
        server = start_gateway(supervisor, config.identity; host=config.listen_host, port=config.port, site, control)
        atexit(shutdown)
        supervised = get(ENV, "LCM_SUPERVISED", "") == "1"
        Base.exit_on_sigint(!supervised)
        if supervised
            @async begin
                try
                    readline(stdin)
                finally
                    shutdown()
                end
            end
        end
        println("LCM runtime listening at $(config.identity.origin)/")
        println("Private identity mode: $(config.identity isa LocalIdentity ? "local development" : "trusted proxy")")
        println("UI hosts: explicit launch only. Scientific resources are not started here.")
        flush(stdout)
        wait(server)
    finally
        shutdown()
    end
    return nothing
end

"""
    serve_agent(config_file)

Start one provisioned, service-managed agent with the shared scientific and
terminal resource owner. Unsupported limits leave no eligible profiles, never an
unconfined fallback. Broker control and expiry remain independent of executors.

The root also watches its verified service incarnation for a stop request, so
normal service shutdown does not depend on which asynchronous task receives
SIGINT. Cleanup runs before process exit; failed cleanup remains an error and
the fixed post-stop command handles crash/forced-exit recovery.
"""
function serve_agent(config_file::AbstractString)
    path = realpath(config_file)
    config = read_agent_config(path)
    driver = nothing
    agent = nothing
    cleaned = false
    function cleanup()
        cleaned && return
        if agent !== nothing
            close(agent)
        elseif driver !== nothing
            close(driver)
        end
        cleaned = true
        println("LCM agent shutdown complete; owned resources retired.")
        flush(stdout)
    end
    Base.exit_on_sigint(false)
    try
        driver = ManagedScientificDriver(path,config)
        agent = AgentService(config,ManagedAgentResources(driver))
        start_agent!(agent)
        println("LCM agent $(config.worker_id) started; $(length(driver.profiles.definitions)) eligible profiles.")
        for (id,reason) in sort!(collect(driver.unavailable);by=first)
            println("Profile $id unavailable: $reason")
        end
        println("Control scheduling active; no executor prepared on startup.")
        flush(stdout)
        object = only(only(systemd_bus_command(driver.runner,["call","org.freedesktop.systemd1",
            "/org/freedesktop/systemd1","org.freedesktop.systemd1.Manager","GetUnit","s",driver.manager.unit])).data)
        # Typed D-Bus state is the existing supervisor's stop request, not a new
        # control file/socket. Monitor separately from heartbeat/lease scheduling.
        while !istaskdone(agent.task)
            managed_agent_stopping(driver.runner,object,driver.manager) && break
            sleep(1)
        end
        istaskdone(agent.task) && wait(agent.task)
    catch error
        interrupted = error isa InterruptException || (error isa TaskFailedException &&
            agent !== nothing && error.task === agent.task &&
            any(entry->entry.exception isa InterruptException,Base.current_exceptions(error.task)))
        interrupted || rethrow()
    finally
        cleanup()
    end
    return nothing
end

"""
    runtime_cli(arguments)

Handle runtime lifecycle, operator provisioning and passive configuration checks.
Configuration checks are read-only and never connect to a broker or start a UI.
"""
function runtime_cli(arguments)
    args = collect(arguments)
    !isempty(args) && first(args) == "runtime" && popfirst!(args)
    if isempty(args) || any(in(("--help", "-h")), args)
        println("Usage: lcm runtime <check|start|status|migrate|permissions|provision|check-agent> --config FILE [--xray]")
        println("       lcm runtime fingerprint --project DIRECTORY")
        println("       lcm runtime check-host [--runtime auto|podman|docker|native]")
        println("       lcm runtime agent-unit --config FILE")
        println("       lcm runtime start-agent --config FILE  (requires the rendered service)")
        println("       lcm runtime recover-agent --journal DIRECTORY --worker ID")
        println("Build first: lcm playground build")
        println("--xray enables owned-component diagnostics for newly launched UI hosts.")
        println("permissions prints server ACLs; provision explicitly creates bounded v2 streams.")
        return nothing
    end
    action = popfirst!(args)
    if action == "recover-agent"
        length(args) == 4 && args[1] == "--journal" && args[3] == "--worker" ||
            throw(ArgumentError("recover-agent requires --journal DIRECTORY --worker ID"))
        recover_agent_resources!(args[2],args[4])
        println("Owned agent recovery complete; no lease or preparation restored.")
        return nothing
    end
    if action == "check-host"
        (isempty(args) || (length(args) == 2 && args[1] == "--runtime")) ||
            throw(ArgumentError("check-host accepts only --runtime auto|podman|docker|native"))
        choice = isempty(args) ? "auto" : args[2]
        runner = CommandRunner()
        try
            check = if choice == "native"
                found = check_native_host(runner)
                println("Native backend: local user-systemd. No resources started.")
                found
            else
                found = check_container_host(runner; requested=choice)
                println("Container engine: $(found.engine.name); rootless=$(found.rootless). No resources started.")
                found
            end
            isempty(check.failures) || throw(ArgumentError("host prerequisites unavailable: " * join(string.(check.failures), ", ")))
            println("Host prerequisites present. Effective executor limits still require launch-time verification.")
        finally
            close(runner)
        end
        return nothing
    end
    if action == "fingerprint"
        length(args) == 2 && args[1] == "--project" ||
            throw(ArgumentError("fingerprint requires exactly --project DIRECTORY"))
        result = native_environment_fingerprint(args[2])
        println(result.digest)
        return nothing
    end
    action in ("check", "start", "status", "migrate", "permissions", "provision", "check-agent", "agent-unit", "start-agent") || throw(ArgumentError("unknown runtime action"))
    path = nothing
    xray = false
    while !isempty(args)
        argument = popfirst!(args)
        if argument == "--config" && !isempty(args) && path === nothing
            path = popfirst!(args)
        elseif argument == "--xray" && action == "start" && !xray
            xray = true
        else
            throw(ArgumentError("unknown or repeated runtime option"))
        end
    end
    path === nothing && throw(ArgumentError("runtime requires --config FILE"))
    if action == "agent-unit"
        print(agent_service_unit(path))
        return nothing
    end
    if action == "start-agent"
        serve_agent(path)
        return nothing
    end
    if action == "check-agent"
        agent = read_agent_config(path)
        println("Agent configuration valid for $(agent.worker_id). No profile was prepared or resource started.")
        return nothing
    end
    config = read_config(path)
    if action == "check"
        println("Runtime configuration valid; enabled=$(config.enabled). No resources started.")
        config.control === nothing || println("Worker control: $(length(config.control.workers)) provisioned identities; enrollment and preparation remain explicit.")
    elseif action in ("permissions", "provision")
        config.control === nothing && throw(ArgumentError("worker control is not configured"))
        ids = sort!(collect(keys(config.control.workers)))
        if action == "permissions"
            println("# Add these server-owned users to the broker authorization block.")
            println(broker_user_config(CoordinatorIdentity();
                password_environment="LCM_V2_COORDINATOR_PASSWORD", worker_ids=ids))
            for id in ids
                reference = "LCM_V2_WORKER_" * uppercase(bytes2hex(sha256(id))[1:16]) * "_PASSWORD"
                println(broker_user_config(WorkerIdentity(id); password_environment=reference))
            end
        else
            jobs = BrokerJobs(config.control.endpoint, CoordinatorIdentity())
            try
                for id in ids
                    ensure_worker_streams!(jobs, config.control.workers[id])
                    println("Verified bounded v2 streams for $id.")
                end
            finally
                close(jobs)
            end
        end
    elseif action == "migrate"
        backup = migrate_runtime!(config.database)
        println(backup === nothing ? "Runtime database is already current." :
            "Runtime database migrated; backup directory beside the database: " * basename(dirname(backup)))
    elseif action == "status"
        url = "http://$(config.listen_host == "::1" ? "[::1]" : config.listen_host):$(config.port)/runtime/api/capabilities"
        response = HTTP.get(url; proxy=nothing, retry=false, redirect=false, request_timeout=5)
        println(String(response.body))
    else
        serve_runtime(config; xray)
    end
    return nothing
end

function (@main)(arguments)
    try
        runtime_cli(arguments)
        return 0
    catch error
        # Never print arbitrary HTTP/Cmd objects containing private headers.
        message = error isa ArgumentError ? error.msg : "Runtime action failed; inspect owned runtime status."
        println(stderr, "lcm runtime: ", message)
        return 2
    end
end
