include("agent_resource_fixture.jl")

@testset "production coordinator and agent schedulers enforce cleanup before replacement reports" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        profiles, applications = ProfileRegistry(), ApplicationRegistry()
        register!(profiles, ProfileDefinition("line-parameters", "/approved", repeat("a",64); operations=("system.echo",)))
        definition = ApplicationDefinition("study", "Study", :workbench, "/study";
            requirements=(RuntimeRequirement("main", ("line-parameters",)),))
        register!(applications, definition)
        operator, alice = Principal("operator"; administrator=true), Principal("alice")
        trusts = [WorkerTrust(id, "credential-$id", ("line-parameters",)) for id in ("worker-a", "worker-b")]
        for trust in trusts
            enroll_worker!(store, operator, trust)
            set_registration_state!(store, operator, trust.worker_id, :approved; expected_revision=1)
        end
        service = ControlService(ControlConfig(endpoint("coordinator"), profiles, trusts), store, applications)
        resources = AgentResourceFixture(profiles)
        agent = AgentService(AgentConfig("worker-a", endpoint("worker-a"), profiles, joinpath(dir, "a")), resources)
        # B is a healthy supervisor with no locally usable profile. Discovery
        # must remain available while matching execution is strictly unavailable.
        b = AgentService(AgentConfig("worker-b", endpoint("worker-b"), profiles, joinpath(dir, "b")),
            AgentResourceFixture(ProfileRegistry()))
        try
            start_agent!(agent); start_agent!(b); start_control!(service)
            @test start_agent!(agent) === agent
            @test resources.recovered == 1
            @test timedwait(() -> length(service.inventory.presence) == 2 &&
                all(p -> p.report.sequence >= 2, values(service.inventory.presence)), 20) == :ok
            @test service.state == :online && agent.state == :online
            @test isempty(service.inventory.presence["worker-b"].report.profiles)
            run = reserve_run!(store, alice, definition)
            @test_throws AccessDenied reserve_assignment!(service.coordinator.assignments, alice,
                run.id, "main", "line-parameters"; placement=PinnedPlacement("worker-b"))
            lease = reserve_assignment!(service.coordinator.assignments, alice, run.id, "main",
                "line-parameters"; placement=PinnedPlacement("worker-a"))
            id = UUID(lease.fence.lease_id)
            grant_assignment!(service.coordinator, alice, id)
            @test timedwait(() -> assignment_usable(service.coordinator, alice, id), 5) == :ok
            resources.release_allowed[] = false
            release_assignment!(service.coordinator, alice, id)
            @test timedwait(() -> !isempty(agent.cleanups), 5) == :ok
            sequence = service.inventory.presence["worker-a"].report.sequence
            @test timedwait(() -> service.inventory.presence["worker-a"].report.sequence > sequence, 5) == :ok
            @test get_assignment(store, alice, id).state == :releasing
            @test !assignment_usable(service.coordinator, alice, id)
            @test isempty(resources.released)
            # Coordinator replacement cannot trust old SQLite authority, nor can
            # the agent report to the replacement before its held cleanup ends.
            close(service)
            replacement = ControlService(service.config, store, applications)
            service = replacement
            start_control!(service)
            @test timedwait(() -> agent.state == :reconciling, 8) == :ok
            @test !haskey(service.inventory.presence, "worker-a")
            @test get_assignment(store, alice, id).state == :releasing
            resources.release_allowed[] = true
            @test timedwait(() -> haskey(service.inventory.presence, "worker-a"), 10) == :ok
            @test get_assignment(store, alice, id).state == :expired
            @test isempty(agent.cleanups) && !isempty(resources.released)
            @test !agent_lease_usable(agent.ledger, lease.fence)

            # A fresh generation can be acknowledged, and ordinary coordinator
            # shutdown should persist its resource-confirmed release.
            next_lease = reserve_assignment!(service.coordinator.assignments, alice, run.id,
                "main", "line-parameters"; placement=PinnedPlacement("worker-a"))
            next_id = UUID(next_lease.fence.lease_id)
            grant_assignment!(service.coordinator, alice, next_id)
            @test timedwait(() -> assignment_usable(service.coordinator, alice, next_id), 5) == :ok
            close(service)
            @test get_assignment(store, alice, next_id).state == :released
        finally
            resources.release_allowed[] = true
            close(service)
            close(agent); close(b)
            close(store)
        end
        @test istaskdone(agent.task) && agent.state == :stopped
        @test isempty(agent.cleanups) && agent.link.control === nothing
        @test resources.closed
    end
end
