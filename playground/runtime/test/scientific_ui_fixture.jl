# Actual registered deck and workbench drivers; no broker or numerical fixture.
using LineCableModelsRuntime, HTTP
length(ARGS) == 2 || error("Expected port and owned temporary directory")
port = parse(Int, ARGS[1])
directory = abspath(ARGS[2])
isdir(directory) || error("Fixture directory must already exist")
store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
owner = Principal("scientific-ui-fixture")
supervisor = UIHostSupervisor(store, default_applications(; xray=true), joinpath(directory,"hosts");
    limits=RunLimits(max_runs=2, max_runs_per_owner=2, startup_seconds=120, shutdown_seconds=1))
site = PublishedSite(normpath(joinpath(@__DIR__,"..","..","_site")))
server = start_gateway(supervisor, LocalIdentity("http://127.0.0.1:$port",owner); port, site)
shutdown = @async begin
    readline(stdin)
    close(server)
end
try
    start_ui!(supervisor,owner,"ichqp-showcase")
    start_ui!(supervisor,owner,"cable-study")
    println("Registered scientific UI fixture listening on $port")
    flush(stdout)
    wait(server)
finally
    close(server)
    close(supervisor)
    close(store)
    wait(shutdown)
end
