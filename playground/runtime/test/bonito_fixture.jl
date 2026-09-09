using LineCableModelsRuntime, HTTP
length(ARGS) == 2 || error("Expected port and an owned temporary directory")
port = parse(Int, ARGS[1])
directory = abspath(ARGS[2])
isdir(directory) || error("Fixture directory must already exist")
store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
registry = ApplicationRegistry()
register!(registry, LocalApplication(
    ApplicationDefinition("bonito-fixture", "Bonito boundary fixture", :workbench, "/counter"),
    dirname(dirname(@__DIR__)), joinpath(@__DIR__, "bonito_child.jl")))
owner = Principal("fixture-owner")
supervisor = UIHostSupervisor(store, registry, joinpath(directory, "hosts");
    limits=RunLimits(max_runs=2, max_runs_per_owner=2, startup_seconds=120, shutdown_seconds=1))
server = start_gateway(supervisor, LocalIdentity("http://127.0.0.1:$port", owner); port)
shutdown = @async begin
    readline(stdin)
    close(server)
end
try
    first = start_ui!(supervisor, owner, "bonito-fixture")
    second = start_ui!(supervisor, owner, "bonito-fixture")
    println("Bonito boundary fixture listening on $port")
    flush(stdout)
    wait(server)
finally
    close(server)
    close(supervisor)
    close(store)
    wait(shutdown)
end
