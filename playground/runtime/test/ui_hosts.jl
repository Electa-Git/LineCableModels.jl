@testset "independent supervised UI processes" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        registry = ApplicationRegistry()
        app = LocalApplication(ApplicationDefinition("mock", "Mock", :workbench, "/mock"),
            dirname(@__DIR__), joinpath(@__DIR__, "ui_child.jl"))
        register!(registry, app)
        supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts");
            limits=RunLimits(max_runs=2, max_runs_per_owner=1, startup_seconds=30, shutdown_seconds=0.5))
        alice, bob = Principal("alice"), Principal("bob")
        try
            @test_throws ArgumentError UIHostSupervisor(store, registry, joinpath(dir, "other"))
            request_id = uuid4()
            first = start_ui!(supervisor, alice, "mock"; request_id)
            second = start_ui!(supervisor, bob, "mock")
            @test first.state == :starting
            @test first.id != second.id
            @test start_ui!(supervisor, alice, "mock"; request_id).id == first.id
            @test_throws CapacityUnavailable start_ui!(supervisor, alice, "mock")
            @test timedwait(() -> get_run(store, alice, first.id).state != :starting, 40) == :ok
            @test timedwait(() -> get_run(store, bob, second.id).state != :starting, 40) == :ok
            @test get_run(store, alice, first.id).state == :running
            @test get_run(store, bob, second.id).state == :running
            @test_throws AccessDenied stop_ui!(supervisor, bob, first.id)
            first_handle, second_handle = supervisor.handles[first.id], supervisor.handles[second.id]
            @test !occursin(first_handle.key, repr(first_handle))
            kill(first_handle.process, Base.SIGKILL)
            @test timedwait(() -> get_run(store, alice, first.id).state == :failed, 10) == :ok
            @test get_run(store, bob, second.id).state == :running
            @test Base.process_running(second_handle.process)
            @test !ispath(joinpath(dir, "hosts", string(first.id)))
            @test stop_ui!(supervisor, bob, second.id).state == :stopped
            @test stop_ui!(supervisor, bob, second.id).state == :stopped
            @test isempty(supervisor.handles)
            @test isempty(readdir(joinpath(dir, "hosts")))
        finally
            close(supervisor)
            close(store)
        end
    end
end
