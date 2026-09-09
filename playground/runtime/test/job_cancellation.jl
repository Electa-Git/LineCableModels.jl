@testset "queued cancellation fences admission across preparation and child replacement" begin
    mktempdir() do dir
        (;driver,resources,agent,clock,fences)=scientific_fixture(dir;capacity=2)
        first,second=fences
        try
            recover_owned!(resources)
            canceled=assigned_fixture_job(resources,first)
            @test !cancel_job!(resources,first,canceled.request.job_id)
            @test !cancel_job!(resources,first,canceled.request.job_id)
            @test isempty(driver.processes) && isempty(resources.handles)
            @test length(resources.canceled_jobs[first.lease_id].ids)==1
            @test_throws RT.ExecutionCore.OperationCanceled execute_assigned!(resources,canceled)
            @test_throws AccessDenied cancel_job!(resources,change_runtime_record(first;owner="foreign"),string(uuid4()))
            @test_throws ArgumentError cancel_job!(resources,first,"not-a-uuid")
            fetch(prepare_assigned!(resources,first,Dict{String,Any}()))
            fetch(prepare_assigned!(resources,second,Dict{String,Any}()))
            canceled=change_runtime_record(canceled;execution=prepared_execution(resources,first))
            @test_throws RT.ExecutionCore.OperationCanceled execute_assigned!(resources,canceled)
            @test fetch(execute_assigned!(resources,assigned_fixture_job(resources,second))).value["value"]==3
            @test release_owned!(resources,first) # replacing an executor is not release of its lease
            @test haskey(resources.canceled_jobs,first.lease_id)
            fetch(prepare_assigned!(resources,first,Dict{String,Any}()))
            canceled=change_runtime_record(canceled;execution=prepared_execution(resources,first))
            @test_throws RT.ExecutionCore.OperationCanceled execute_assigned!(resources,canceled)
            slow=assigned_fixture_job(resources,first,"fixture.delay",Dict("seconds"=>10))
            running=execute_assigned!(resources,slow)
            @test timedwait(()->scientific_status(resources,first).phase==:executing,5)==:ok
            @test cancel_job!(resources,first,slow.request.job_id)
            @test_throws TaskFailedException fetch(running)
            @test scientific_status(resources,first).failure=="canceled"
            @test scientific_status(resources,second).preparation==:ready
            @test !cancel_job!(resources,first,slow.request.job_id)
            for _ in 1:254
                cancel_job!(resources,first,string(uuid4()))
            end
            @test length(resources.canceled_jobs[first.lease_id].ids)==256
            @test_throws AccessDenied cancel_job!(resources,first,string(uuid4()))
            @test !cancel_job!(resources,first,slow.request.job_id) # bounded history retains existing IDs
            clock[]=6
            @test_throws AccessDenied cancel_job!(resources,first,string(uuid4()))
            release_owned!(resources,first)
            @test isempty(resources.canceled_jobs)
        finally
            close(agent)
        end
        @test isempty(resources.canceled_jobs) && isempty(driver.processes)
    end
end

@testset "job cancellation travels through the existing ordered scientific channel" begin
    mktempdir() do dir
        (;resources,agent,fences)=scientific_fixture(dir;capacity=1)
        fence=only(fences)
        reports=RT.Protocol.ScientificReport[]
        emit=report->push!(reports,report)
        try
            recover_owned!(resources)
            job=assigned_fixture_job(resources,fence)
            command=scientific_command(fence,1,"cancel_job";target_id=job.request.job_id)
            receive_scientific!(agent.science,command;emit)
            @test last(reports).accepted && last(reports).reason=="cancel_recorded"
            @test last(reports).preparation=="cold"
            @test_throws RT.ExecutionCore.OperationCanceled execute_assigned!(resources,job)
            receive_scientific!(agent.science,command;emit)
            @test last(reports).accepted && length(resources.canceled_jobs[fence.lease_id].ids)==1
            successor=change_runtime_record(command;revision=2)
            receive_scientific!(agent.science,successor;emit)
            @test last(reports).accepted && length(resources.canceled_jobs[fence.lease_id].ids)==1
            receive_scientific!(agent.science,command;emit)
            @test !last(reports).accepted && last(reports).reason=="stale_command"
        finally
            close(agent)
        end
    end
end
