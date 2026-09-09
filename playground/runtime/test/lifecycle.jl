struct RejectedApplication <: AbstractApplication end
RT.describe(::RejectedApplication) = ApplicationDefinition("rejected", "Rejected launch", :workbench, "/mock")
RT.ui_command(::RejectedApplication, ::HostContext) = error("This trusted hook rejected launch")

@testset "one invalid receipt cannot prevent sibling teardown" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        registry = ApplicationRegistry()
        register!(registry, LocalApplication(
            ApplicationDefinition("mock", "Mock", :workbench, "/mock"),
            dirname(@__DIR__), joinpath(@__DIR__, "ui_child.jl")))
        supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts");
            limits=RunLimits(startup_seconds=30, shutdown_seconds=0.1))
        alice = Principal("alice")
        try
            runs = [start_ui!(supervisor, alice, "mock") for _ in 1:2]
            @test timedwait(() -> all(r -> get_run(store, alice, r.id).state == :running, runs), 40) == :ok
            handles = [supervisor.handles[r.id] for r in runs]
            marker = joinpath(first(handles).directory, "ownership.json")
            receipt = read(marker)
            write(marker, "{}")
            close(supervisor)
            @test all(h -> !Base.process_running(h.process), handles)
            @test isempty(supervisor.handles)
            @test isdir(first(handles).directory)
            @test !ispath(last(handles).directory)
            @test occursin("operator attention", get_run(store, alice, first(runs).id).reason)
            recovered = UIHostSupervisor(store, registry, joinpath(dir, "hosts"))
            try
                @test isdir(first(handles).directory) # recovery preserves unverified data
                write(marker, receipt)
                RT.remove_ui_directory!(recovered, first(runs).id)
                @test isempty(readdir(recovered.root))
            finally
                close(recovered)
            end
        finally
            close(supervisor); close(store)
        end
    end
end

@testset "bounded failure and disconnected lifetime" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        registry = ApplicationRegistry()
        register!(registry, RejectedApplication())
        register!(registry, LocalApplication(
            ApplicationDefinition("unready", "Unready", :workbench, "/mock"),
            dirname(@__DIR__), joinpath(@__DIR__, "unready_child.jl")))
        supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts");
            limits=RunLimits(max_runs=1, max_runs_per_owner=1, startup_seconds=0.3, shutdown_seconds=0.1))
        alice = Principal("alice")
        try
            for _ in 1:3
                @test_throws ArgumentError start_ui!(supervisor, alice, "rejected")
                @test isempty(supervisor.handles)
                @test isempty(readdir(supervisor.root))
                run = start_ui!(supervisor, alice, "unready")
                @test timedwait(() -> get_run(store, alice, run.id).state == :failed, 10) == :ok
                @test get_run(store, alice, run.id).reason == "UI host startup timed out"
                @test isempty(supervisor.handles)
                @test isempty(readdir(supervisor.root))
            end
        finally
            close(supervisor); close(store)
        end
    end
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        registry = ApplicationRegistry()
        register!(registry, LocalApplication(
            ApplicationDefinition("mock", "Mock", :workbench, "/mock"),
            dirname(@__DIR__), joinpath(@__DIR__, "ui_child.jl")))
        supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts");
            limits=RunLimits(startup_seconds=30, shutdown_seconds=0.1, disconnect_grace_seconds=1))
        alice = Principal("alice")
        try
            run = start_ui!(supervisor, alice, "mock")
            @test timedwait(() -> get_run(store, alice, run.id).state == :running, 40) == :ok
            @test timedwait(() -> get_run(store, alice, run.id).state == :stopped, 10) == :ok
            @test occursin("Disconnected", get_run(store, alice, run.id).reason)
            @test isempty(readdir(supervisor.root))
        finally
            close(supervisor); close(store)
        end
    end
end

@testset "coordinator death, child death and receipt-based recovery" begin
    mktempdir() do dir
        command = `$(Base.julia_cmd()) --startup-file=no --threads=2 --project=$(dirname(@__DIR__)) $(joinpath(@__DIR__, "coordinator_child.jl")) $dir`
        parent = run(pipeline(ignorestatus(command); stdout=devnull, stderr=devnull); wait=false)
        try
            receipt = joinpath(dir, "coordinator-ready.json")
            @test timedwait(() -> isfile(receipt), 60) == :ok
            child = JSON3.read(read(receipt, String))
            kill(parent, Base.SIGKILL)
            wait(parent)
            # Observe only this fixture's child; never reconstruct a kill handle
            # from a persisted PID. Zombies are already dead and need adoption.
            child_dead() = begin
                status = "/proc/$(child.pid)/stat"
                !isfile(status) || occursin(r"\) Z ", read(status, String))
            end
            @test timedwait(child_dead, 10) == :ok
            store = RuntimeStore(joinpath(dir, "recovery.sqlite"))
            registry = ApplicationRegistry()
            supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts"))
            try
                recovered = get_run(store, Principal("alice"), UUID(child.run_id))
                @test recovered.state == :failed
                @test occursin("coordinator restart", recovered.reason)
                @test isempty(supervisor.handles)
                @test isempty(readdir(supervisor.root))
            finally
                close(supervisor); close(store)
            end
        finally
            Base.process_running(parent) && kill(parent, Base.SIGKILL)
            wait(parent)
        end
    end
end
