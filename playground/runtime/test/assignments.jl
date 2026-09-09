function allocation_fixture(f; capacity=2, limits=AssignmentLimits(), profile_kind=:scientific)
    mktempdir() do directory
        store = RuntimeStore(joinpath(directory, "runtime.sqlite"))
        profiles = ProfileRegistry()
        profile = profile_kind == :terminal ?
            ProfileDefinition("line-parameters", "fixture@sha256:" * repeat("a",64), repeat("a",64);
                kind=:terminal, isolation=:container) :
            ProfileDefinition("line-parameters", "/approved/project", repeat("a",64); operations=("system.echo",))
        register!(profiles, profile)
        applications = ApplicationRegistry()
        definition = ApplicationDefinition("study", "Study", :workbench, "/study";
            requirements=Tuple(RuntimeRequirement(role, ("line-parameters",)) for role in ("main", "aux", "third")))
        register!(applications, definition)
        admin, alice, bob = Principal("operator"; administrator=true), Principal("alice"), Principal("bob")
        a_run = reserve_run!(store, alice, definition)
        b_run = reserve_run!(store, bob, definition)
        clock = Ref(0.0)
        inventory = WorkerInventory(store, profiles; clock=()->clock[])
        function announce(id; sequence=1, boot=string(uuid4()))
            probe = probe_worker!(inventory, id)
            accept_report!(inventory, id, RT.Protocol.WorkerAnnouncement("2.0", id, boot,
                inventory.coordinator_id, probe.challenge, sequence, capacity,
                [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))]))
        end
        for id in ("worker-a", "worker-b")
            enroll_worker!(store, admin, WorkerTrust(id, "credential-$id", ("line-parameters",); capacity))
            set_registration_state!(store, admin, id, :approved; expected_revision=1)
            announce(id)
        end
        manager = AssignmentManager(inventory, applications; limits)
        try
            f((; store, profiles, applications, definition, admin, alice, bob, a_run, b_run,
                clock, inventory, announce, manager))
        finally
            close(store)
        end
    end
end

@testset "assignments authorize and reserve without pretending to be ready" begin
    allocation_fixture() do fixture
        (; store, alice, bob, admin, a_run, b_run, manager, inventory) = fixture
        request = uuid4()
        @test_throws AccessDenied reserve_assignment!(manager, bob, a_run.id, "main", "line-parameters")
        @test_throws AccessDenied reserve_assignment!(manager, alice, a_run.id, "missing", "line-parameters")
        @test_throws AccessDenied reserve_assignment!(manager, alice, a_run.id, "main", "not-permitted")
        assigned = reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters"; request_id=request)
        @test assigned.state == :reserving && assigned.revision == 0
        @test assigned.fence.owner == "alice" && assigned.fence.run_id == string(a_run.id)
        @test assigned.fence.worker_id == "worker-a" && assigned.fence.generation == 1
        @test assigned.fence.worker_boot == inventory.presence["worker-a"].report.boot_id
        @test reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters";
            request_id=request).fence == assigned.fence
        @test_throws AccessDenied reserve_assignment!(manager, alice, a_run.id, "aux", "line-parameters";
            request_id=request)
        @test_throws AccessDenied reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters")
        @test_throws AccessDenied get_assignment(store, bob, UUID(assigned.fence.lease_id))
        @test isempty(list_assignments(store, bob))
        @test_throws AccessDenied list_assignments(store, bob; run_id=a_run.id)
        @test only(list_assignments(store, admin)).fence == assigned.fence
        # Operator action is still owned/charged to Bob, with stable retries.
        operator_request = uuid4()
        other = reserve_assignment!(manager, admin, b_run.id, "main", "line-parameters"; request_id=operator_request)
        @test other.fence.owner == "bob" && other.fence.worker_id == "worker-b"
        @test reserve_assignment!(manager, admin, b_run.id, "main", "line-parameters";
            request_id=operator_request).fence == other.fence
        @test length(list_assignments(store, admin)) == 2
    end
end

@testset "pinned loss and dedicated allocation cannot silently share or fall back" begin
    allocation_fixture(; capacity=3) do fixture
        (; manager, alice, bob, a_run, b_run, store, clock, announce) = fixture
        first = reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters";
            placement=DedicatedPlacement("worker-a"))
        @test first.fence.worker_id == "worker-a" && first.placement == :dedicated
        @test_throws AccessDenied reserve_assignment!(manager, bob, b_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        sibling = reserve_assignment!(manager, alice, a_run.id, "aux", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        @test sibling.fence.worker_id == "worker-a"
        other = reserve_assignment!(manager, bob, b_run.id, "main", "line-parameters")
        @test other.fence.worker_id == "worker-b"
        @test_throws AccessDenied reserve_assignment!(manager, alice, a_run.id, "third", "line-parameters";
            placement=DedicatedPlacement("worker-b"))
        # Release only fixture bookkeeping; executable cleanup is tested by the
        # lease lifecycle, not asserted by this allocation-only test.
        RT.sql_rows(store.db, "UPDATE leases SET state='released'")
        clock[] = 5
        announce("worker-b"; sequence=2, boot=fixture.inventory.presence["worker-b"].report.boot_id)
        @test_throws AccessDenied reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        replacement = reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters")
        @test replacement.fence.worker_id == "worker-b" && replacement.fence.generation == 2
        @test replacement.fence.lease_id != first.fence.lease_id
    end
end

@testset "reserved and restart-unreconciled capacity is still occupied" begin
    allocation_fixture(; capacity=1) do fixture
        (; manager, store, profiles, applications, alice, bob, a_run, b_run, clock) = fixture
        reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters"; placement=PinnedPlacement("worker-a"))
        @test_throws AccessDenied reserve_assignment!(manager, bob, b_run.id, "main", "line-parameters";
            placement=PinnedPlacement("worker-a"))
        another_store = RuntimeStore(store.path)
        try
            restarted = WorkerInventory(another_store, profiles; clock=()->clock[])
            second = AssignmentManager(restarted, applications)
            @test_throws AccessDenied reserve_assignment!(second, bob, b_run.id, "main", "line-parameters")
            probe = probe_worker!(restarted, "worker-a")
            accept_report!(restarted, "worker-a", RT.Protocol.WorkerAnnouncement("2.0", "worker-a", string(uuid4()),
                restarted.coordinator_id, probe.challenge, 1, 1,
                [RT.Protocol.ProfileAdvertisement("line-parameters", "1.0.0", repeat("a", 64))]))
            # A fresh report alone does not erase the previous reservation.
            @test_throws AccessDenied reserve_assignment!(second, bob, b_run.id, "main", "line-parameters";
                placement=PinnedPlacement("worker-a"))
            @test length(list_assignments(another_store, fixture.admin)) == 1
        finally
            close(another_store)
        end
        set_registration_state!(store, fixture.admin, "worker-b", :draining; expected_revision=2)
        @test_throws AccessDenied reserve_assignment!(manager, bob, b_run.id, "main", "line-parameters")
    end
    allocation_fixture(; limits=AssignmentLimits(total=2, per_owner=1, per_run=1)) do fixture
        (; manager, alice, a_run) = fixture
        reserve_assignment!(manager, alice, a_run.id, "main", "line-parameters")
        @test_throws CapacityUnavailable reserve_assignment!(manager, alice, a_run.id, "aux", "line-parameters")
    end
    @test_throws ArgumentError AssignmentLimits(total=1, per_owner=2)
    @test_throws ArgumentError AssignmentLimits(per_run=true)
    @test_throws ArgumentError PinnedPlacement("worker.*")
end

@testset "schema two migration retains approval and idempotent enrollment" begin
    mktempdir() do directory
        path = joinpath(directory, "runtime.sqlite")
        store = RuntimeStore(path)
        admin = Principal("operator"; administrator=true)
        request = uuid4()
        trust = WorkerTrust("worker-a", "credential-a", ("line-parameters",))
        enroll_worker!(store, admin, trust; request_id=request)
        set_registration_state!(store, admin, "worker-a", :approved; expected_revision=1)
        RT.sql_rows(store.db, "DROP TABLE job_cancellations")
        RT.sql_rows(store.db, "DROP TABLE jobs")
        RT.sql_rows(store.db, "DROP TABLE leases")
        RT.sql_rows(store.db, "PRAGMA user_version=2")
        close(store)
        @test_throws ArgumentError RuntimeStore(path)
        backup = migrate_runtime!(path)
        @test occursin(".schema-2-backup-", backup)
        old = SQLite.DB(backup)
        try
            @test only(RT.sql_rows(old, "PRAGMA user_version")).user_version == 2
            @test only(RT.sql_rows(old, "SELECT state FROM workers")).state == "approved"
        finally
            close(old)
        end
        current = RuntimeStore(path)
        try
            @test only(list_registrations(current, admin)).state == :approved
            @test enroll_worker!(current, admin, trust; request_id=request).state == :pending
            @test only(list_registrations(current, admin)).state == :approved
            @test isempty(list_assignments(current, admin))
        finally
            close(current)
        end
        @test migrate_runtime!(path) === nothing
    end
end

@testset "independent processes cannot oversubscribe one worker slot" begin
    allocation_fixture(; capacity=1) do fixture
        barrier = mktempdir(dirname(fixture.store.path); prefix="allocation-race-")
        processes = Base.Process[]
        logs = IO[]
        try
            for (principal, app_run) in ((fixture.alice, fixture.a_run), (fixture.bob, fixture.b_run))
                logfile = open(joinpath(barrier, principal.id * ".log"), "w")
                push!(logs, logfile)
                project = normpath(joinpath(@__DIR__, ".."))
                script = joinpath(@__DIR__, "allocation_race.jl")
                boot = fixture.inventory.presence["worker-a"].report.boot_id
                command = `$(Base.julia_cmd()) --startup-file=no --compiled-modules=existing --project=$project $script $(fixture.store.path) $(principal.id) $(string(app_run.id)) $boot $barrier`
                push!(processes, run(pipeline(command; stdin=devnull, stdout=logfile, stderr=logfile); wait=false))
            end
            @test timedwait(() -> all(isfile(joinpath(barrier, id * ".ready")) for id in ("alice", "bob")),
                45; pollint=0.025) == :ok
            write(joinpath(barrier, "go"), "go")
            @test timedwait(() -> all(process_exited, processes), 30; pollint=0.025) == :ok
            @test all(success, processes)
            outcomes = sort([read(joinpath(barrier, id * ".result"), String) for id in ("alice", "bob")])
            @test outcomes == ["reserved", "unavailable"]
            @test length(list_assignments(fixture.store, fixture.admin)) == 1
        finally
            for process in processes
                if process_running(process)
                    kill(process, Base.SIGTERM)
                    timedwait(() -> process_exited(process), 3; pollint=0.025)
                    process_running(process) && kill(process, Base.SIGKILL)
                end
                wait(process)
            end
            foreach(close, logs)
        end
    end
end
