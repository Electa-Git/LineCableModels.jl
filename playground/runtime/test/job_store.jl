function receipt_fixture_job(f;value=1)
    P=RT.Protocol
    P.AssignedJob("2.0",f.lease.fence,
        P.new_job_request("system.echo",Dict("value"=>value);session_id=f.lease.fence.run_id),
        P.PreparedExecution(string(uuid4()),1,repeat("b",64)))
end

@testset "job receipts preserve request identity and bound admission across connections" begin
    coordinated_lease_fixture() do f
        (;store,alice,bob,coordinator,agent,id)=f
        job=receipt_fixture_job(f)
        request=uuid4()
        @test_throws AccessDenied reserve_job!(store,alice,job,request)
        grant=grant_assignment!(coordinator,alice,id)
        @test accept_lease_ack!(coordinator,"worker-a",handle_lease_control!(agent,grant))
        first=reserve_job!(store,alice,job,request)
        @test first.state==:queued && first.job==job && first.request_id==request
        # A lost HTTP reply must not replace deadline, job ID or preparation.
        retry=reserve_job!(store,alice,receipt_fixture_job(f),request)
        @test retry.job==first.job && retry.created_at==first.created_at
        @test prior_job_submission(store,alice,id,"system.echo",Dict("value"=>1),request).job==job
        @test prior_job_submission(store,alice,id,"system.echo",Dict("value"=>1),uuid4())===nothing
        @test_throws AccessDenied reserve_job!(store,alice,receipt_fixture_job(f;value=2),request)
        @test_throws AccessDenied prior_job_submission(store,alice,id,"system.echo",Dict("value"=>2),request)
        @test_throws AccessDenied prior_job_submission(store,bob,id,"system.echo",Dict("value"=>1),request)
        @test_throws AccessDenied get_job(store,bob,UUID(job.request.job_id))
        @test_throws AccessDenied get_job(store,alice,uuid4())
        @test_throws AccessDenied list_jobs(store,bob,UUID(job.fence.run_id))
        @test_throws AccessDenied reserve_job!(store,bob,job,request)
        @test_throws AccessDenied reserve_job!(store,alice,receipt_fixture_job(f),uuid4())
        other=RuntimeStore(store.path)
        try
            @test get_job(other,alice,UUID(job.request.job_id)).job==job
            @test_throws AccessDenied reserve_job!(other,alice,receipt_fixture_job(f),uuid4())
            @test transition_job!(other,alice,UUID(job.request.job_id),:submitted).state==:submitted
            @test transition_job!(store,alice,UUID(job.request.job_id),:submitted).state==:submitted
            @test transition_job!(store,alice,UUID(job.request.job_id),:succeeded).state==:succeeded
            @test_throws AccessDenied transition_job!(store,alice,UUID(job.request.job_id),:submitted)
            @test_throws AccessDenied transition_job!(store,bob,UUID(job.request.job_id),:failed)
            @test_throws ArgumentError transition_job!(store,alice,UUID(job.request.job_id),:arbitrary)
            @test reserve_job!(other,alice,receipt_fixture_job(f),request).job==job
            second=reserve_job!(other,alice,receipt_fixture_job(f;value=2),uuid4())
            @test second.state==:queued && second.job.request.job_id!=job.request.job_id
            @test length(list_jobs(store,alice,UUID(job.fence.run_id)))==2
            # Losing authority cannot revive a saved submission or admit a new
            # one. Historical receipts remain readable without broker access.
            release_assignment!(coordinator,alice,id)
            @test reserve_job!(store,alice,receipt_fixture_job(f),request).job==job
            transition_job!(store,alice,UUID(second.job.request.job_id),:revoked)
            @test_throws AccessDenied reserve_job!(store,alice,receipt_fixture_job(f;value=3),uuid4())
            @test get_job(store,alice,UUID(second.job.request.job_id)).state==:revoked
        finally
            close(other)
        end
        @test_throws AccessDenied RT.job_submission_hash(id,"system.echo",Dict("large"=>repeat("x",65536)))
        @test_throws ArgumentError RT.job_submission_hash(id,"system.echo",Dict("callable"=>identity))
    end
end

@testset "job receipt history has a finite per-assignment lifetime" begin
    coordinated_lease_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        for index in 1:256
            receipt=reserve_job!(f.store,f.alice,receipt_fixture_job(f;value=index),uuid4())
            transition_job!(f.store,f.alice,UUID(receipt.job.request.job_id),:failed)
        end
        @test length(list_jobs(f.store,f.alice,f.a_run.id))==256
        @test_throws AccessDenied reserve_job!(f.store,f.alice,receipt_fixture_job(f),uuid4())
    end
end

@testset "schema three migration retains existing run and lease bookkeeping" begin
    coordinated_lease_fixture() do f
        path=f.store.path
        RT.sql_rows(f.store.db,"DROP TABLE job_cancellations")
        RT.sql_rows(f.store.db,"DROP TABLE jobs")
        RT.sql_rows(f.store.db,"PRAGMA user_version=3")
        @test_throws ArgumentError RuntimeStore(path)
        backup=migrate_runtime!(path)
        @test occursin(".schema-3-backup-",backup) && stat(backup).mode & 0o777==0o600
        old=SQLite.DB(backup)
        try
            @test only(RT.sql_rows(old,"PRAGMA user_version")).user_version==3
            @test only(RT.sql_rows(old,"SELECT lease_id FROM leases")).lease_id==string(f.id)
            @test isempty(RT.sql_rows(old,"SELECT name FROM sqlite_master WHERE name='jobs'"))
        finally
            close(old)
        end
        current=RuntimeStore(path)
        try
            @test get_assignment(current,f.alice,f.id).fence==f.lease.fence
            @test get_run(current,f.alice,f.a_run.id).id==f.a_run.id
            @test isempty(list_jobs(current,f.alice,f.a_run.id))
            @test only(RT.sql_rows(current.db,"PRAGMA user_version")).user_version==RT.RUN_SCHEMA_VERSION
        finally
            close(current)
        end
        @test migrate_runtime!(path)===nothing
    end
end
