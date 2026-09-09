@testset "operator enrollment and durable idempotent decisions" begin
    mktempdir() do directory
        path = joinpath(directory, "runtime.sqlite")
        store = RuntimeStore(path)
        admin = Principal("operator"; administrator=true)
        alice = Principal("alice")
        trust = WorkerTrust("worker-a", "credential-a", ("line-parameters",); capacity=2)
        request = uuid4()
        try
            @test_throws AccessDenied enroll_worker!(store, alice, trust)
            pending = enroll_worker!(store, admin, trust; request_id=request)
            @test pending.state == :pending && pending.revision == 1
            @test enroll_worker!(store, admin, trust; request_id=request).created_at == pending.created_at
            @test_throws AccessDenied enroll_worker!(store, admin, trust)
            @test_throws AccessDenied enroll_worker!(store, admin,
                WorkerTrust("worker-b", "credential-a", ("line-parameters",)))
            @test_throws AccessDenied enroll_worker!(store, admin,
                WorkerTrust("worker-b", "credential-b", ("line-parameters",)); request_id=request)
            @test length(list_registrations(store, alice)) == 1
            @test_throws AccessDenied set_registration_state!(store, alice, "worker-a", :approved; expected_revision=1)
            @test_throws AccessDenied set_registration_state!(store, admin, "worker-a", :draining; expected_revision=1)
            approval = uuid4()
            active = set_registration_state!(store, admin, "worker-a", :approved;
                expected_revision=1, request_id=approval)
            @test active.state == :approved && active.revision == 2
            drained = set_registration_state!(store, admin, "worker-a", :draining; expected_revision=2)
            @test drained.state == :draining && drained.revision == 3
            @test_throws AccessDenied set_registration_state!(store, admin, "worker-a", :approved; expected_revision=2)
            replay = set_registration_state!(store, admin, "worker-a", :approved;
                expected_revision=1, request_id=approval)
            @test replay.state == :approved && replay.revision == 2
            @test only(list_registrations(store, alice)).state == :draining # replay has no side effect
            @test_throws AccessDenied set_registration_state!(store, admin, "worker-a", :disabled;
                expected_revision=3, request_id=approval)
            @test set_registration_state!(store, admin, "worker-a", :disabled; expected_revision=3).revision == 4
        finally
            close(store)
        end
        reopened = RuntimeStore(path)
        try
            @test only(list_registrations(reopened, alice)).state == :disabled
            @test only(list_registrations(reopened, alice)).profiles == ("line-parameters",)
        finally
            close(reopened)
        end
    end
    @test_throws ArgumentError WorkerTrust("worker.*", "identity", ("line-parameters",))
    @test_throws ArgumentError WorkerTrust("worker-a", "identity", ("line-parameters",); capacity=true)
    @test_throws ArgumentError WorkerTrust("worker-a", "identity", ("line-parameters", "line-parameters"))
end

@testset "explicit backed-up schema migration excludes live supervisors" begin
    mktempdir() do directory
        path = joinpath(directory, "runtime.sqlite")
        store = RuntimeStore(path)
        principal = Principal("alice")
        saved = reserve_run!(store, principal, ApplicationDefinition("mock", "Mock", :workbench, "/mock"))
        supervisor = UIHostSupervisor(store, ApplicationRegistry(), joinpath(directory, "hosts"))
        try
            @test_throws ArgumentError migrate_runtime!(path)
        finally
            close(supervisor)
        end
        # Reproduce the exact earlier schema: only run bookkeeping existed.
        RT.sql_rows(store.db, "DROP TABLE job_cancellations")
        RT.sql_rows(store.db, "DROP TABLE jobs")
        RT.sql_rows(store.db, "DROP TABLE leases")
        RT.sql_rows(store.db, "DROP TABLE worker_actions")
        RT.sql_rows(store.db, "DROP TABLE workers")
        RT.sql_rows(store.db, "PRAGMA user_version = 1")
        close(store)
        @test_throws ArgumentError RuntimeStore(path)
        backup = migrate_runtime!(path)
        @test isfile(backup) && stat(backup).mode & 0o777 == 0o600
        old = SQLite.DB(backup)
        try
            @test only(RT.sql_rows(old, "PRAGMA user_version")).user_version == 1
            @test only(RT.sql_rows(old, "SELECT run_id FROM runs")).run_id == string(saved.id)
            @test isempty(RT.sql_rows(old, "SELECT name FROM sqlite_master WHERE name='workers'"))
        finally
            close(old)
        end
        migrated = RuntimeStore(path)
        try
            @test get_run(migrated, principal, saved.id).application == "mock"
            @test isempty(list_registrations(migrated, principal))
            @test only(RT.sql_rows(migrated.db, "PRAGMA user_version")).user_version == RT.RUN_SCHEMA_VERSION
        finally
            close(migrated)
        end
        @test migrate_runtime!(path) === nothing
        @test length(filter(name -> occursin(".schema-1-backup-", name), readdir(directory))) == 1
    end
end
