mutable struct LeaseTestTransport
    messages::Vector{RT.Protocol.LeaseControl}
    available::Bool
end
function RT.send_control!(transport::LeaseTestTransport, command::RT.Protocol.LeaseControl)
    transport.available || throw(RT.BrokerUnavailable())
    push!(transport.messages, command)
    return nothing
end

function coordinated_lease_fixture(f; profile_kind=:scientific, duration_ms=5000)
    allocation_fixture(; capacity=1, profile_kind) do fixture
        transport = LeaseTestTransport(RT.Protocol.LeaseControl[], true)
        coordinator = LeaseCoordinator(fixture.manager, transport; duration_ms, ack_seconds=1)
        agent = AgentLeaseLedger("worker-a", fixture.profiles; clock=()->fixture.clock[])
        probe = probe_worker!(fixture.inventory, "worker-a")
        receive_probe!(agent, probe)
        report = RT.Protocol.WorkerAnnouncement("2.0", "worker-a", agent.boot_id,
            fixture.inventory.coordinator_id, probe.challenge, 1, 1,
            [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))])
        reconcile_worker_report!(coordinator, "worker-a", report)
        lease = reserve_assignment!(fixture.manager, fixture.alice, fixture.a_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        f((; fixture..., transport, coordinator, agent, lease, id=UUID(lease.fence.lease_id)))
    end
end

@testset "durable reservation becomes usable only through its exact acknowledged control" begin
    coordinated_lease_fixture() do f
        (; coordinator, agent, lease, id, alice, bob, clock, store, transport) = f
        @test !assignment_usable(coordinator, alice, id)
        @test_throws AccessDenied grant_assignment!(coordinator, bob, id)
        request = uuid4()
        grant = grant_assignment!(coordinator, alice, id; request_id=request)
        @test grant.revision == 1 && grant.fence == lease.fence
        @test !assignment_usable(coordinator, alice, id)
        ack = handle_lease_control!(agent, grant)
        @test !accept_lease_ack!(coordinator, "worker-b", ack)
        @test !accept_lease_ack!(coordinator, "worker-a", change_runtime_record(ack; request_id=string(uuid4())))
        @test accept_lease_ack!(coordinator, "worker-a", ack)
        @test assignment_usable(coordinator, alice, id) && agent_lease_usable(agent, lease.fence)
        @test !accept_lease_ack!(coordinator, "worker-a", ack)
        @test get_assignment(store, alice, id).state == :active
        clock[] = 0.5
        sent = length(transport.messages)
        @test grant_assignment!(coordinator, alice, id; request_id=request) == grant
        @test length(transport.messages) == sent # acknowledged retry is inert
        renewal = renew_assignment!(coordinator, alice, id)
        @test renewal.revision == 2
        @test assignment_usable(coordinator, alice, id) # prior authority survives while renewal is pending
        renewed = handle_lease_control!(agent, renewal)
        @test accept_lease_ack!(coordinator, "worker-a", renewed)
        @test coordinator.flights[id].authority_until == 5.5
        @test !accept_lease_ack!(coordinator, "worker-a", ack)
        release = release_assignment!(coordinator, alice, id)
        @test !assignment_usable(coordinator, alice, id)
        @test handle_lease_control!(agent, release) === nothing
        @test get_assignment(store, alice, id).state == :releasing
        @test_throws AccessDenied reserve_assignment!(f.manager, bob, f.b_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        finished = complete_agent_cleanup!(agent, lease.fence)
        @test accept_lease_ack!(coordinator, "worker-a", finished)
        @test get_assignment(store, alice, id).state == :released
        @test reserve_assignment!(f.manager, bob, f.b_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a")).fence.owner == "bob"
    end
end

@testset "lost or reordered grants cannot free capacity or revive released authority" begin
    coordinated_lease_fixture() do f
        grant = grant_assignment!(f.coordinator, f.alice, f.id)
        f.clock[] = 1
        tick_leases!(f.coordinator)
        release = last(f.transport.messages)
        @test release.action == "release" && release.revision == 2
        @test !assignment_usable(f.coordinator, f.alice, f.id)
        # Release overtakes the undelivered grant. Its tombstone must survive
        # after the coordinator makes the physical slot available again.
        finished = handle_lease_control!(f.agent, release)
        @test finished.accepted
        @test accept_lease_ack!(f.coordinator, "worker-a", finished)
        @test !handle_lease_control!(f.agent, grant).accepted
        @test !agent_lease_usable(f.agent, f.lease.fence)
        tick_leases!(f.coordinator)
        @test isempty(f.coordinator.flights)
    end
    coordinated_lease_fixture() do f
        grant = grant_assignment!(f.coordinator, f.alice, f.id)
        accepted_but_lost = handle_lease_control!(f.agent, grant)
        f.clock[] = 1.1
        @test !accept_lease_ack!(f.coordinator, "worker-a", accepted_but_lost) # late ACK cannot revive authority
        tick_leases!(f.coordinator)
        release = last(f.transport.messages)
        @test release.action == "release"
        @test handle_lease_control!(f.agent, release) === nothing
        finished = complete_agent_cleanup!(f.agent, f.lease.fence)
        # Deliberately lose the cleanup reply. Coordinator still reserves capacity.
        @test get_assignment(f.store, f.alice, f.id).state == :releasing
        @test_throws AccessDenied reserve_assignment!(f.manager, f.bob, f.b_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        f.clock[] = 2.2
        tick_leases!(f.coordinator)
        @test last(f.transport.messages) == release
        @test handle_lease_control!(f.agent, release) == finished
        @test accept_lease_ack!(f.coordinator, "worker-a", finished)
        @test get_assignment(f.store, f.alice, f.id).state == :released
    end
end

@testset "broker outage retains reservations and retries only cleanup" begin
    coordinated_lease_fixture() do f
        f.transport.available = false
        grant_assignment!(f.coordinator, f.alice, f.id)
        @test isempty(f.transport.messages)
        f.clock[] = 1.1
        tick_leases!(f.coordinator)
        @test get_assignment(f.store, f.alice, f.id).state == :releasing
        @test isempty(f.transport.messages)
        f.transport.available = true
        f.clock[] = 2.2
        tick_leases!(f.coordinator)
        command = only(f.transport.messages)
        @test command.action == "release" # no delayed grant replay
        @test accept_lease_ack!(f.coordinator, "worker-a", handle_lease_control!(f.agent, command))
        @test !assignment_usable(f.coordinator, f.alice, f.id)
    end
end

@testset "coordinator restart reconciles cleanup instead of trusting persisted active rows" begin
    coordinated_lease_fixture() do f
        command = grant_assignment!(f.coordinator, f.alice, f.id)
        acknowledgement = handle_lease_control!(f.agent, command)
        accept_lease_ack!(f.coordinator, "worker-a", acknowledgement)
        restarted_inventory = WorkerInventory(f.store, f.profiles; clock=()->f.clock[])
        restarted = LeaseCoordinator(AssignmentManager(restarted_inventory, f.applications),
            LeaseTestTransport(RT.Protocol.LeaseControl[], true))
        @test get_assignment(f.store, f.alice, f.id).state == :active
        @test !assignment_usable(restarted, f.alice, f.id)
        @test !accept_lease_ack!(restarted, "worker-a", acknowledgement)
        probe = probe_worker!(restarted_inventory, "worker-a")
        @test !receive_probe!(f.agent, probe)
        @test !agent_lease_usable(f.agent, f.lease.fence)
        complete_agent_cleanup!(f.agent, f.lease.fence)
        @test receive_probe!(f.agent, probe)
        report = RT.Protocol.WorkerAnnouncement("2.0", "worker-a", f.agent.boot_id,
            probe.coordinator_id, probe.challenge, 2, 1,
            [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))])
        reconcile_worker_report!(restarted, "worker-a", report)
        @test get_assignment(f.store, f.alice, f.id).state == :expired
        next = reserve_assignment!(restarted.assignments, f.alice, f.a_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        @test next.fence.generation == 2 && next.fence.coordinator_id == probe.coordinator_id
        @test !handle_lease_control!(f.agent, command).accepted
    end
end

@testset "abandoned reservation and stopped run both initiate bounded release" begin
    coordinated_lease_fixture() do f
        tick_leases!(f.coordinator)
        @test isempty(f.transport.messages)
        f.clock[] = 1.1
        tick_leases!(f.coordinator)
        release = only(f.transport.messages)
        @test release.action == "release" && release.revision == 1
        @test get_assignment(f.store, f.alice, f.id).state == :releasing
        @test isempty(f.coordinator.unissued)
        @test accept_lease_ack!(f.coordinator, "worker-a", handle_lease_control!(f.agent, release))
        tick_leases!(f.coordinator)
        @test isempty(f.coordinator.flights) && isempty(f.coordinator.unissued)
    end
    coordinated_lease_fixture() do f
        grant = grant_assignment!(f.coordinator, f.alice, f.id)
        accept_lease_ack!(f.coordinator, "worker-a", handle_lease_control!(f.agent, grant))
        transition_run!(f.store, f.alice, f.a_run.id, :stopped)
        @test !assignment_usable(f.coordinator, f.alice, f.id)
        tick_leases!(f.coordinator)
        @test last(f.transport.messages).action == "release"
        @test get_assignment(f.store, f.alice, f.id).state == :releasing
    end
end
