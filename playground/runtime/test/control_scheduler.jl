@testset "scheduled coordinator probes, renews and expires over actual TLS" begin
    mktempdir() do state_directory
        store = RuntimeStore(joinpath(state_directory, "runtime.sqlite"))
        profiles, applications = ProfileRegistry(), ApplicationRegistry()
        register!(profiles, ProfileDefinition("line-parameters", "/approved/environment", repeat("a",64);
            operations=("system.echo",)))
        definition = ApplicationDefinition("study", "Study", :workbench, "/study";
            requirements=(RuntimeRequirement("main", ("line-parameters",)),))
        register!(applications, definition)
        trust = WorkerTrust("worker-a", "credential-a", ("line-parameters",))
        operator, alice = Principal("operator"; administrator=true), Principal("alice")
        enroll_worker!(store, operator, trust)
        set_registration_state!(store, operator, "worker-a", :approved; expected_revision=1)
        service = ControlService(ControlConfig(endpoint("coordinator"), profiles, [trust]), store, applications)
        agent_wire = BrokerControl(endpoint("worker-a"), WorkerIdentity("worker-a"))
        ledger = AgentLeaseLedger("worker-a", profiles)
        running, paused = Ref(true), Ref(false)
        sequence = Ref(0)
        # This fixture exercises the production coordinator scheduler against an
        # actual authenticated agent connection. It owns no scientific executor.
        task = @async begin
            while running[]
                if !paused[]
                    for envelope in poll_control!(agent_wire)
                        record = envelope.record
                        if record isa P.WorkerProbe
                            if receive_probe!(ledger, record)
                                sequence[] += 1
                                send_control!(agent_wire, P.WorkerAnnouncement("2.0", "worker-a", ledger.boot_id,
                                    record.coordinator_id, record.challenge, sequence[], 1,
                                    [P.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a",64))]))
                            end
                        else
                            ack = handle_lease_control!(ledger, record)
                            ack === nothing || send_control!(agent_wire, ack)
                        end
                    end
                end
                for fence in expire_agent_leases!(ledger)
                    ack = complete_agent_cleanup!(ledger, fence) # no resources in this fixture
                    ack === nothing || send_control!(agent_wire, ack)
                end
                sleep(0.01)
            end
        end
        try
            @test start_control!(service) === service
            original_task = service.task
            @test start_control!(service).task === original_task
            @test timedwait(() -> !isempty(service.inventory.presence) &&
                first(values(service.inventory.presence)).report.sequence >= 2, 20) == :ok
            @test service.state == :online
            @test only(RT.control_snapshot(service, alice).workers).liveness == "online"
            @test only(RT.control_snapshot(service, alice).workers).preparation == "unknown"
            run = reserve_run!(store, alice, definition)
            lease = reserve_assignment!(service.coordinator.assignments, alice, run.id, "main", "line-parameters")
            id = UUID(lease.fence.lease_id)
            grant_assignment!(service.coordinator, alice, id)
            @test timedwait(() -> assignment_usable(service.coordinator, alice, id), 5) == :ok
            @test timedwait(() -> get_assignment(store, alice, id).revision >= 2 &&
                !service.coordinator.flights[id].pending, 8) == :ok
            @test agent_lease_usable(ledger, lease.fence)
            @test assignment_usable(service.coordinator, alice, id)

            # Pausing control reception prevents new probes/renewals from being
            # acknowledged. Local deadlines still revoke authority and capacity
            # is not free while the remote release has not been acknowledged.
            paused[] = true
            @test timedwait(() -> !assignment_usable(service.coordinator, alice, id), 12) == :ok
            @test get_assignment(store, alice, id).state in (:active, :releasing, :reconciling)
            @test timedwait(() -> !agent_lease_usable(ledger, lease.fence), 12) == :ok
            @test_throws AccessDenied reserve_assignment!(service.coordinator.assignments, alice,
                run.id, "main", "line-parameters"; placement=PinnedPlacement("worker-a"))
            paused[] = false
            @test timedwait(() -> get_assignment(store, alice, id).state == :released, 10) == :ok
            @test timedwait(() -> only(RT.control_snapshot(service, alice).workers).liveness == "online", 10) == :ok
            @test length(control_events(service.events, alice).records) <= service.events.capacity
            @test !occursin("fixture-password", JSON3.write(control_events(service.events, alice)))
            transition_run!(store, alice, run.id, :stopped; reason="fixture completed")
        finally
            close(service)
            running[] = false
            wait(task)
            close(agent_wire)
            close(store)
        end
        @test istaskdone(service.task) && service.state == :stopped
        @test service.link.control === nothing
    end
end
