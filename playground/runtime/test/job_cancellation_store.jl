@testset "durable cancellation intent is idempotent and not a terminal job result" begin
    coordinated_lease_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        receipt=reserve_job!(f.store,f.alice,receipt_fixture_job(f),uuid4())
        id=UUID(receipt.job.request.job_id)
        request=uuid4()
        @test job_cancellation(f.store,f.alice,id)===nothing
        @test_throws AccessDenied request_job_cancellation!(f.store,f.bob,id,request)
        intent=request_job_cancellation!(f.store,f.alice,id,request)
        @test intent.job_id==id && intent.request_id==request && !intent.acknowledged
        @test request_job_cancellation!(f.store,f.alice,id,uuid4()).request_id==request
        other=RuntimeStore(f.store.path)
        try
            @test job_cancellation(other,f.alice,id).request_id==request
            @test_throws AccessDenied job_cancellation(other,f.bob,id)
            @test_throws AccessDenied RT.acknowledge_job_cancellation!(other,f.alice,id,uuid4())
            @test RT.acknowledge_job_cancellation!(other,f.alice,id,request).acknowledged
            @test get_job(other,f.alice,id).state==:queued # an ACK is not a result
            transition_job!(other,f.alice,id,:canceled)
            next=reserve_job!(other,f.alice,receipt_fixture_job(f),uuid4())
            @test_throws AccessDenied request_job_cancellation!(other,f.alice,UUID(next.job.request.job_id),request)
            @test request_job_cancellation!(other,f.alice,id,request).acknowledged
            transition_job!(other,f.alice,UUID(next.job.request.job_id),:succeeded)
            @test request_job_cancellation!(other,f.alice,UUID(next.job.request.job_id),uuid4())===nothing
        finally
            close(other)
        end
    end
end

@testset "schema four migration preserves exact submission receipts" begin
    coordinated_lease_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        receipt=reserve_job!(f.store,f.alice,receipt_fixture_job(f),uuid4())
        RT.sql_rows(f.store.db,"DROP TABLE job_cancellations")
        RT.sql_rows(f.store.db,"PRAGMA user_version=4")
        @test_throws ArgumentError RuntimeStore(f.store.path)
        backup=migrate_runtime!(f.store.path)
        @test occursin(".schema-4-backup-",backup) && stat(backup).mode & 0o777==0o600
        old=SQLite.DB(backup)
        try
            @test only(RT.sql_rows(old,"PRAGMA user_version")).user_version==4
            @test only(RT.sql_rows(old,"SELECT job_json FROM jobs")).job_json==RT.Protocol.encode_message(receipt.job)
        finally
            close(old)
        end
        current=RuntimeStore(f.store.path)
        try
            @test get_job(current,f.alice,UUID(receipt.job.request.job_id)).job==receipt.job
            @test job_cancellation(current,f.alice,UUID(receipt.job.request.job_id))===nothing
            @test only(RT.sql_rows(current.db,"PRAGMA user_version")).user_version==5
        finally
            close(current)
        end
    end
end
