using Test, UUIDs, Dates, LineCableModelsRuntime
import NATS, JSON3, Logging
const RT = LineCableModelsRuntime
const P = RT.Protocol
directory, port, artifact_port = ARGS
certs = joinpath(directory, "certs")
endpoint(id) = BrokerEndpoint("tls://127.0.0.1:$port", joinpath(directory, id * ".password");
    ca_file=joinpath(certs, "ca.pem"), certificate_file=joinpath(certs, "worker-cert.pem"),
    key_file=joinpath(certs, "worker-key.pem"), server_name="localhost")
artifact_location(id)=S3RuntimeArtifacts("https://127.0.0.1:$artifact_port","lcm-runtime-private","runtime-v1",
    joinpath(directory,"artifact-$id.toml");ca_file=joinpath(certs,"ca.pem"))
function records(control; timeout=5)
    found = RT.ControlEnvelope[]
    timedwait(timeout; pollint=0.01) do
        append!(found, poll_control!(control))
        !isempty(found)
    end
    return found
end
function pull(connection, id; timeout=1)
    NATS.request(connection, "\$JS.API.CONSUMER.MSG.NEXT.$(RT.job_stream(id)).agent",
        JSON3.write((batch=1, no_wait=true)); timeout)
end
function raw_connection(user, password)
    Logging.with_logger(Logging.NullLogger()) do
        NATS.connect("tls://127.0.0.1:$port"; user, pass=password,
            auth_token=nothing, jwt=nothing, nkey=nothing, nkey_seed=nothing,
            tls_ca_path=joinpath(certs, "ca.pem"), tls_cert_path=joinpath(certs, "worker-cert.pem"),
            tls_key_path=joinpath(certs, "worker-key.pem"), tls_server_name="localhost",
            retry_on_init_fail=false, connect_timeout=3.0, send_enqueue_when_disconnected=false,
            drain_timeout=0.5, drain_poll=0.01)
    end
end

@testset "actual TLS control subjects bind distinct worker credentials" begin
    coordinator = BrokerControl(endpoint("coordinator"), CoordinatorIdentity(); worker_ids=("worker-a", "worker-b"))
    a = b = legacy = data = nothing
    try
        a = BrokerControl(endpoint("worker-a"), WorkerIdentity("worker-a"))
        b = BrokerControl(endpoint("worker-b"), WorkerIdentity("worker-b"))
        @test_throws RT.BrokerUnavailable RT.connect_broker(endpoint("worker-a"), WorkerIdentity("worker-b"))
        probe = P.WorkerProbe("2.0", "worker-a", string(uuid4()), string(uuid4()))
        send_control!(coordinator, probe)
        @test only(records(a)).record == probe
        @test isempty(poll_control!(b))
        boot = string(uuid4())
        installed = [P.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))]
        report = P.WorkerAnnouncement("2.0", "worker-a", boot, probe.coordinator_id,
            probe.challenge, 1, 1, installed)
        send_control!(a, report)
        envelope = only(records(coordinator))
        @test envelope.worker_id == "worker-a" && envelope.record == report

        # Bypass our typed sender deliberately: the broker must reject A's
        # publication on B's authenticated identity subject.
        NATS.publish(a.connection, "lcm.report.v2.worker-b", P.encode_message(report))
        NATS.ping(a.connection; measure=false)
        @test isempty(records(coordinator; timeout=0.2))
        # The owned receiver also rejects a foreign payload on A's legal subject.
        foreign = P.WorkerAnnouncement("2.0", "worker-b", boot, probe.coordinator_id,
            probe.challenge, 2, 1, installed)
        NATS.publish(a.connection, "lcm.report.v2.worker-a", P.encode_message(foreign))
        @test isempty(records(coordinator; timeout=0.2))
        @test coordinator.rejected == 1

        fence = P.AssignmentFence(string(uuid4()), string(uuid4()), "alice", "parameters",
            "worker-a", boot, probe.coordinator_id, "line-parameters", "1.0.0", repeat("a", 64), 1)
        grant = P.LeaseControl("2.0", string(uuid4()), "grant", fence, 1, 10_000)
        send_control!(coordinator, grant)
        @test only(records(a)).record == grant
        @test isempty(poll_control!(b))
        ack = P.LeaseAcknowledgement("2.0", grant.request_id, fence, 1, true, "accepted")
        send_control!(a, ack)
        @test only(records(coordinator)).record == ack
        @test_throws AccessDenied send_control!(b, ack)
        NATS.publish(a.connection, "lcm.control.v2.worker-a.lease", P.encode_message(grant))
        NATS.ping(a.connection; measure=false)
        @test isempty(records(a; timeout=0.2)) # agents cannot grant themselves authority

        # Separate streams prevent a consumer API from bypassing subject ACLs.
        data = BrokerJobs(endpoint("coordinator"), CoordinatorIdentity())
        for id in ("worker-a", "worker-b")
            trust = WorkerTrust(id, "credential-$id", ("line-parameters",))
            ensure_worker_streams!(data, trust)
            @test ensure_worker_streams!(data, trust) === nothing
        end
        NATS.JetStream.stream_publish(coordinator.connection, P.assigned_job_subject(fence), "target-a")
        @test_throws NATS.NATSError pull(b.connection, "worker-a"; timeout=0.2)
        first_delivery = pull(a.connection, "worker-a")
        @test NATS.payload(first_delivery) == "target-a"
        # An explicit negative acknowledgement must re-deliver immediately,
        # not silently become an empty positive ACK.
        NATS.JetStream.consumer_ack(a.connection, first_delivery, "-NAK")
        redelivery = pull(a.connection, "worker-a")
        @test NATS.payload(redelivery) == "target-a"
        @test redelivery.reply_to != first_delivery.reply_to
        NATS.JetStream.consumer_ack(a.connection, redelivery, "+TERM")
        @test_throws NATS.NATSError pull(a.connection, "worker-a")

        # A worker cannot manufacture a consumer with a broader filter.
        @test_throws NATS.NATSError NATS.request(a.connection,
            "\$JS.API.CONSUMER.CREATE.$(RT.job_stream("worker-b")).stolen",
            JSON3.write((stream_name=RT.job_stream("worker-b"),
                config=(name="stolen", ack_policy="explicit", filter_subjects=["lcm.jobs.v2.worker-b.>"]))); timeout=0.2)
        legacy = raw_connection("worker", "legacy-worker-fixture-password")
        @test_throws NATS.NATSError pull(legacy, "worker-a"; timeout=0.2)
        @test_throws NATS.NATSError NATS.request(legacy,
            "\$JS.API.CONSUMER.CREATE.$(RT.job_stream("worker-a")).legacy",
            "{}"; timeout=0.2)

        # Invalid frames consume only the bounded poll budget.
        for _ in 1:10
            NATS.publish(a.connection, "lcm.report.v2.worker-a", "{}")
        end
        send_control!(b, P.WorkerAnnouncement("2.0", "worker-b", string(uuid4()),
            probe.coordinator_id, string(uuid4()), 1, 1, installed))
        observed = records(coordinator)
        @test any(item -> item.worker_id == "worker-b", observed)
        @test coordinator.rejected >= 2
        @test isempty(poll_control!(coordinator))
    finally
        data === nothing || close(data)
        legacy === nothing || NATS.drain(legacy)
        b === nothing || close(b)
        a === nothing || close(a)
        close(coordinator)
        close(coordinator)
    end
    @test coordinator.closed && isempty(coordinator.subscriptions)
    @test NATS.status(coordinator.connection) == NATS.DRAINED
    @test isempty(coordinator.connection.sub_data)
end

@testset "SQLite reservation and agent lease authority round-trip through actual TLS" begin
    mktempdir() do state_directory
        store = RuntimeStore(joinpath(state_directory, "runtime.sqlite"))
        coordinator_wire = worker_wire = nothing
        try
            coordinator_wire = BrokerControl(endpoint("coordinator"), CoordinatorIdentity(); worker_ids=("worker-a",))
            worker_wire = BrokerControl(endpoint("worker-a"), WorkerIdentity("worker-a"))
            profiles, applications = ProfileRegistry(), ApplicationRegistry()
            register!(profiles, ProfileDefinition("line-parameters", "/approved/project", repeat("a", 64);
                operations=("system.echo",)))
            definition = ApplicationDefinition("study", "Study", :workbench, "/study";
                requirements=(RuntimeRequirement("main", ("line-parameters",)),))
            register!(applications, definition)
            operator, alice = Principal("operator"; administrator=true), Principal("alice")
            enroll_worker!(store, operator, WorkerTrust("worker-a", "credential-a", ("line-parameters",)))
            set_registration_state!(store, operator, "worker-a", :approved; expected_revision=1)
            run = reserve_run!(store, alice, definition)
            inventory = WorkerInventory(store, profiles)
            coordinator = LeaseCoordinator(AssignmentManager(inventory, applications), coordinator_wire)
            agent = AgentLeaseLedger("worker-a", profiles)
            sequence = Ref(0)
            function heartbeat()
                probe = probe_worker!(inventory, "worker-a")
                send_control!(coordinator_wire, probe)
                received = only(records(worker_wire)).record
                @test receive_probe!(agent, received)
                sequence[] += 1
                report = P.WorkerAnnouncement("2.0", "worker-a", agent.boot_id, probe.coordinator_id,
                    received.challenge, sequence[], 1,
                    [P.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))])
                send_control!(worker_wire, report)
                envelope = only(records(coordinator_wire))
                reconcile_worker_report!(coordinator, envelope.worker_id, envelope.record)
            end
            heartbeat()
            lease = reserve_assignment!(coordinator.assignments, alice, run.id, "main", "line-parameters";
                placement=PinnedPlacement("worker-a"))
            id = UUID(lease.fence.lease_id)
            @test !assignment_usable(coordinator, alice, id)
            # Refresh after cold reservation compilation; preparation/readiness
            # must not be inferred from the elapsed startup work in this test.
            heartbeat()
            grant_assignment!(coordinator, alice, id)
            grant = only(records(worker_wire)).record
            send_control!(worker_wire, handle_lease_control!(agent, grant))
            acknowledgement = only(records(coordinator_wire))
            @test accept_lease_ack!(coordinator, acknowledgement.worker_id, acknowledgement.record)
            @test assignment_usable(coordinator, alice, id)
            @test agent_lease_usable(agent, lease.fence)

            renew_assignment!(coordinator, alice, id)
            renewal = only(records(worker_wire)).record
            @test renewal.action == "renew"
            send_control!(worker_wire, handle_lease_control!(agent, renewal))
            envelope = only(records(coordinator_wire))
            @test accept_lease_ack!(coordinator, envelope.worker_id, envelope.record)
            @test assignment_usable(coordinator, alice, id)

            release_assignment!(coordinator, alice, id)
            release = only(records(worker_wire)).record
            @test handle_lease_control!(agent, release) === nothing
            @test !assignment_usable(coordinator, alice, id)
            @test !agent_lease_usable(agent, lease.fence)
            @test get_assignment(store, alice, id).state == :releasing
            # This transport fixture owns no executor. The real supervisor must
            # perform teardown before this same completion call in the agent.
            send_control!(worker_wire, complete_agent_cleanup!(agent, lease.fence))
            envelope = only(records(coordinator_wire))
            @test accept_lease_ack!(coordinator, envelope.worker_id, envelope.record)
            @test get_assignment(store, alice, id).state == :released
            tick_leases!(coordinator)
            @test isempty(coordinator.flights)
        finally
            worker_wire === nothing || close(worker_wire)
            coordinator_wire === nothing || close(coordinator_wire)
            close(store)
        end
    end
end

include("assigned_jobs.jl")
include("control_scheduler.jl")
include("agent_scheduler_tls.jl")
include("artifacts_tls.jl")
include("science_tls.jl")
include("jobs_tls.jl")
