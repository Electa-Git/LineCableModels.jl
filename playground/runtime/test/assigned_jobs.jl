# Executable transport tests use a synthetic scientific result, not an engine.
@testset "targeted durable jobs retain exact authority and result-before-ack recovery" begin
    mktempdir() do directory
        store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
        connections = Any[]
        try
            cw = BrokerControl(endpoint("coordinator"), CoordinatorIdentity(); worker_ids=("worker-a",))
            push!(connections, cw)
            aw = BrokerControl(endpoint("worker-a"), WorkerIdentity("worker-a"))
            push!(connections, aw)
            cd = BrokerJobs(endpoint("coordinator"), CoordinatorIdentity())
            push!(connections, cd)
            ad = BrokerJobs(endpoint("worker-a"), WorkerIdentity("worker-a"))
            push!(connections, ad)
            bd = BrokerJobs(endpoint("worker-b"), WorkerIdentity("worker-b"))
            push!(connections, bd)
            profiles, applications = ProfileRegistry(), ApplicationRegistry()
            register!(profiles, ProfileDefinition("line-parameters", "/approved/project", repeat("a", 64);
                operations=("system.echo",)))
            definition = ApplicationDefinition("study", "Study", :workbench, "/study";
                requirements=(RuntimeRequirement("main", ("line-parameters",)),))
            register!(applications, definition)
            operator, alice, bob = Principal("operator"; administrator=true), Principal("alice"), Principal("bob")
            trust = WorkerTrust("worker-a", "credential-a", ("line-parameters",))
            enroll_worker!(store, operator, trust)
            set_registration_state!(store, operator, "worker-a", :approved; expected_revision=1)
            @test ensure_worker_streams!(cd, trust) === nothing
            @test_throws AccessDenied ensure_worker_streams!(cd,
                WorkerTrust("worker-a", "credential-a", ("line-parameters",); capacity=2))
            clock = Ref(0.0)
            inventory = WorkerInventory(store, profiles; clock=()->clock[])
            coordinator = LeaseCoordinator(AssignmentManager(inventory, applications), cw)
            agent = AgentLeaseLedger("worker-a", profiles; clock=()->clock[])
            probe = probe_worker!(inventory, "worker-a")
            send_control!(cw, probe)
            @test receive_probe!(agent, only(records(aw)).record)
            send_control!(aw, P.WorkerAnnouncement("2.0", "worker-a", agent.boot_id,
                probe.coordinator_id, probe.challenge, 1, 1,
                [P.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))]))
            envelope = only(records(cw))
            reconcile_worker_report!(coordinator, envelope.worker_id, envelope.record)
            app_run = reserve_run!(store, alice, definition)
            lease = reserve_assignment!(coordinator.assignments, alice, app_run.id, "main", "line-parameters";
                placement=PinnedPlacement("worker-a"))
            id = UUID(lease.fence.lease_id)
            request = P.new_job_request("system.echo", Dict("value"=>1); session_id=string(app_run.id))
            execution=P.PreparedExecution(string(uuid4()),1,repeat("b",64))
            job = P.AssignedJob("2.0", lease.fence, request,execution)
            @test_throws AccessDenied publish_assigned_job!(cd, coordinator, alice, job)
            grant_assignment!(coordinator, alice, id)
            send_control!(aw, handle_lease_control!(agent, only(records(aw)).record))
            envelope = only(records(cw))
            @test accept_lease_ack!(coordinator, envelope.worker_id, envelope.record)
            @test_throws AccessDenied publish_assigned_job!(cd, coordinator, bob, job)
            forbidden = P.AssignedJob("2.0", lease.fence, P.new_job_request("system.delay", Dict("value"=>1);
                session_id=string(app_run.id)),execution)
            @test_throws AccessDenied publish_assigned_job!(cd, coordinator, alice, forbidden)
            acknowledgement = publish_assigned_job!(cd, coordinator, alice, job)
            @test acknowledgement.stream == RT.job_stream("worker-a")
            @test publish_assigned_job!(cd, coordinator, alice, job).duplicate === true
            b_agent = AgentLeaseLedger("worker-b", profiles)
            @test poll_assigned_job!(bd, b_agent) === nothing
            delivery = poll_assigned_job!(ad, agent)
            @test delivery isa AssignedDelivery && delivery.job == job
            @test assigned_result(ad, lease.fence, request.job_id) === nothing
            function outcome_for(job, value)
                result = P.JobResult("1.0", job.request.job_id, job.request.operation, "1.0",
                    job.request.input_hash, "transport-fixture", job.fence.fingerprint, job.fence.worker_id,
                    "miss", P.utc_timestamp(), P.utc_timestamp(), Dict{String,Any}("value"=>value),
                    nothing, nothing, String[])
                P.AssignedResult("2.0", job.fence, result,job.execution)
            end
            outcome = outcome_for(job, 1)
            @test persist_assigned_result!(ad, agent, delivery, outcome) == outcome
            @test assigned_result(ad, lease.fence, request.job_id) == outcome
            @test assigned_result(cd, lease.fence, request.job_id) == outcome
            @test_throws AccessDenied assigned_result(bd, lease.fence, request.job_id)
            @test_throws NATS.NATSError NATS.request(bd.connection,
                "\$JS.API.DIRECT.GET.$(RT.result_stream("worker-a"))",
                JSON3.write((last_by_subj=P.assigned_result_subject(lease.fence, request.job_id),)); timeout=0.2)
            @test poll_assigned_job!(ad, agent) === nothing

            # Persist a second result but deliberately lose the input ACK.
            next_job = P.AssignedJob("2.0", lease.fence,
                P.new_job_request("system.echo", Dict("value"=>2); session_id=string(app_run.id)),execution)
            publish_assigned_job!(cd, coordinator, alice, next_job)
            first_delivery = poll_assigned_job!(ad, agent)
            saved = outcome_for(next_job, 2)
            subject = P.assigned_result_subject(lease.fence, next_job.request.job_id)
            RT.job_request(ad, NATS.JetStream.PubAck, subject,
                (P.encode_message(saved), ["Nats-Msg-Id"=>subject]))
            NATS.JetStream.consumer_ack(ad.connection, first_delivery.message, "-NAK")
            redelivery = poll_assigned_job!(ad, agent)
            @test redelivery.job == next_job
            @test assigned_result(ad, lease.fence, next_job.request.job_id) == saved
            @test persist_assigned_result!(ad, agent, redelivery, outcome_for(next_job, 999)) == saved
            @test assigned_result(cd, lease.fence, next_job.request.job_id) == saved
            @test poll_assigned_job!(ad, agent) === nothing

            # Malformed requests are terminated, never handed to an executor.
            NATS.JetStream.stream_publish(cd.connection, P.assigned_job_subject(lease.fence), "{")
            @test poll_assigned_job!(ad, agent) === nothing
            stale = P.AssignedJob("2.0", lease.fence,
                P.new_job_request("system.echo", Dict("value"=>3); session_id=string(app_run.id)),execution)
            publish_assigned_job!(cd, coordinator, alice, stale)
            clock[] = 11
            @test poll_assigned_job!(ad, agent) === nothing
            @test poll_assigned_job!(ad, agent) === nothing
            @test_throws AccessDenied publish_assigned_job!(cd, coordinator, alice, stale)
            @test_throws AccessDenied persist_assigned_result!(ad, agent, redelivery, saved)
            @test assigned_result(cd, lease.fence, request.job_id) == outcome # history, not fresh authority
            @test get_assignment(store, alice, id).state == :active # durable row alone proves nothing
            @test !assignment_usable(coordinator, alice, id)
        finally
            foreach(close, reverse(connections))
            close(store)
        end
    end
end
