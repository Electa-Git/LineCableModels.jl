# Real registered UI drivers, isolated database, no broker or scientific work.
using LineCableModelsRuntime, HTTP
port = parse(Int, ARGS[1])
directory = abspath(ARGS[2])
isdir(directory) || error("Expected an existing owned temporary directory")
store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
owner = Principal("ui-parity-fixture")
supervisor = UIHostSupervisor(store, default_applications(), joinpath(directory, "hosts");
    limits=RunLimits(max_runs=4, max_runs_per_owner=4, startup_seconds=180,
        disconnect_grace_seconds=600, shutdown_seconds=2))
site = PublishedSite(normpath(joinpath(@__DIR__, "..", "..", "_site")))
server = start_gateway(supervisor, LocalIdentity("http://127.0.0.1:$port", owner); port, site)
shutdown = @async begin
    readline(stdin)
    close(server)
end
try
    for application in ("ichqp-showcase", "cable-study", "template-workbench", "toolkit-gallery")
        # Stagger compiler startup; the comparison needs four concurrent hosts,
        # not four simultaneous cold compilations on a developer laptop.
        run = start_ui!(supervisor, owner, application)
        while get_run(store, owner, run.id).state in (:reserved, :starting)
            sleep(.1)
        end
        get_run(store, owner, run.id).state == :running || error("UI fixture failed: $application")
    end
    println("UI parity fixture listening on $port")
    flush(stdout)
    wait(server)
finally
    close(server)
    close(supervisor)
    close(store)
    wait(shutdown)
end
