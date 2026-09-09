function job_coordinator_fixture(f)
    coordinated_lease_fixture() do fixture
        endpoint=BrokerEndpoint("tls://broker.invalid","/unused-password")
        config=ControlConfig(endpoint,fixture.profiles,[WorkerTrust("worker-a","credential-worker-a",("line-parameters",))])
        science=ScientificCoordinator(endpoint,fixture.coordinator,["worker-a"])
        jobs=JobCoordinator(config,fixture.coordinator,science,ControlEvents())
        try
            f((;fixture...,science,jobs))
        finally
            close(jobs)
            close(science)
        end
    end
end

@testset "job coordinator is inert, owner-fenced and reconciles authority without a broker" begin
    job_coordinator_fixture() do f
        (;jobs,science,store,alice,bob,id)=f
        @test jobs.connection===nothing && jobs.task===nothing && jobs.connector===nothing
        @test isempty(jobs.flights) && jobs.state==:offline && jobs.result_readers==0 && jobs.artifact_readers==0
        @test_throws AccessDenied submit_job!(jobs,bob,id,"system.echo",Dict("value"=>1))
        @test_throws RT.BrokerUnavailable submit_job!(jobs,alice,id,"system.echo",Dict("value"=>1))
        grant=grant_assignment!(f.coordinator,alice,id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        receipt=reserve_job!(store,alice,receipt_fixture_job(f),uuid4())
        jobid=UUID(receipt.job.request.job_id)
        @test_throws AccessDenied RT.owned_job_artifact(jobs,bob,jobid)
        @test_throws ArtifactUnavailable RT.owned_job_artifact(jobs,alice,jobid)
        @test submit_job!(jobs,alice,id,"system.echo",Dict("value"=>1);request_id=receipt.request_id).job==receipt.job
        @test_throws AccessDenied submit_job!(jobs,alice,id,"system.echo",Dict("value"=>2);request_id=receipt.request_id)
        @test_throws AccessDenied owned_job_result(jobs,bob,jobid)
        @test_throws RT.BrokerUnavailable owned_job_result(jobs,alice,jobid)
        @test jobs.result_readers==0
        @test_throws AccessDenied cancel_job!(jobs,bob,jobid)
        @test cancel_job!(jobs,alice,jobid).state==:queued
        @test !job_cancellation(store,alice,jobid).acknowledged
        @test length(control_events(jobs.events,alice).records)==1
        @test isempty(control_events(jobs.events,bob).records)
        event=only(control_events(jobs.events,alice).records)
        @test event.job_id==string(jobid) && event.executor_id==receipt.job.execution.executor_id
        @test event.stage=="cancel_requested" && !haskey(event,:owner) && !haskey(event,:parameters)
        @test_throws ArgumentError RT.record_event!(jobs.events,:job_canceled)
        @test_throws ArgumentError RT.record_event!(jobs.events,:control_connected;job=receipt.job)
        @test_throws RT.BrokerUnavailable RT.reconcile_job!(jobs,jobid)
        @test get_job(store,alice,jobid).state==:queued
        f.clock[]=6 # lease loss must progress even when no job connection exists
        RT.tick_jobs!(jobs)
        @test timedwait(()->get_job(store,alice,jobid).state==:revoked,2)==:ok
        @test get_job(store,alice,jobid).job==receipt.job
        @test !job_cancellation(store,alice,jobid).acknowledged
        @test !RT.current_job_assignment(jobs,receipt)
        @test submit_job!(jobs,alice,id,"system.echo",Dict("value"=>1);request_id=receipt.request_id).state==:revoked
        @test cancel_job!(jobs,alice,jobid).state==:revoked
        @test science.task===nothing && isempty(science.lanes)
    end
end

@testset "job owner drains active authorized readers before closing" begin
    job_coordinator_fixture() do f
        f.jobs.result_readers=1
        f.jobs.artifact_readers=1
        closing=@async close(f.jobs)
        @test timedwait(()->f.jobs.closed,1)==:ok
        @test !istaskdone(closing)
        lock(f.jobs.lock) do; f.jobs.result_readers=0; end
        yield()
        @test !istaskdone(closing)
        lock(f.jobs.lock) do; f.jobs.artifact_readers=0; end
        @test timedwait(()->istaskdone(closing),2)==:ok
        @test fetch(closing)===nothing && f.jobs.state==:stopped
        @test close(f.jobs)===nothing
    end
end

@testset "uncertain job receipt accepts later evidence but cannot be replayed" begin
    job_coordinator_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        receipt=reserve_job!(f.store,f.alice,receipt_fixture_job(f),uuid4())
        id=UUID(receipt.job.request.job_id)
        @test transition_job!(f.store,f.alice,id,:uncertain).state==:uncertain
        @test_throws AccessDenied transition_job!(f.store,f.alice,id,:submitted)
        @test submit_job!(f.jobs,f.alice,f.id,"system.echo",Dict("value"=>1);request_id=receipt.request_id).state==:uncertain
        @test RT.reconcile_job!(f.jobs,id)===nothing
        @test transition_job!(f.store,f.alice,id,:succeeded).state==:succeeded
        @test_throws AccessDenied transition_job!(f.store,f.alice,id,:uncertain)
        @test_throws AccessDenied transition_job!(f.store,f.alice,id,:failed)
    end
end
