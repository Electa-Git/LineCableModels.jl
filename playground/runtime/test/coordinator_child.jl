using LineCableModelsRuntime, JSON3
directory = only(ARGS)
store = RuntimeStore(joinpath(directory, "recovery.sqlite"))
registry = ApplicationRegistry()
register!(registry, LocalApplication(
    ApplicationDefinition("mock", "Mock", :workbench, "/mock"),
    dirname(@__DIR__), joinpath(@__DIR__, "ui_child.jl")))
supervisor = UIHostSupervisor(store, registry, joinpath(directory, "hosts");
    limits=RunLimits(startup_seconds=30))
owner = Principal("alice")
try
    run = start_ui!(supervisor, owner, "mock")
    timedwait(() -> get_run(store, owner, run.id).state == :running, 40) == :ok || error("Child not ready")
    receipt = joinpath(directory, "coordinator-ready.json")
    write(receipt * ".pending", JSON3.write((run_id=string(run.id), pid=getpid(supervisor.handles[run.id].process))))
    mv(receipt * ".pending", receipt)
    sleep(600)
finally
    close(supervisor)
    close(store)
end
