@testset "durable execution owner is inert and binds exact result provenance" begin
    mktempdir() do directory
        (;driver,resources,agent,fences)=scientific_fixture(directory;capacity=1)
        job=assigned_fixture_job(resources,only(fences))
        try
            @test agent.jobs isa AgentJobService
            @test agent.jobs.resources===resources && agent.jobs.ledger===agent.ledger
            @test agent.jobs.connection===nothing && agent.jobs.task===nothing
            @test isempty(driver.processes) && isempty(agent.jobs.flights)
            @test_throws ArgumentError start_jobs!(agent.jobs)
            @test_throws ArgumentError AgentJobService(agent.config.endpoint,resources,
                AgentLeaseLedger("worker-a",driver.profiles))
            started=RT.Protocol.utc_timestamp()
            output=ScientificOutput(Dict{String,Any}("value"=>7),"1.2",["fixture warning"])
            result=RT.job_outcome(job,started,output)
            @test result.execution==job.execution && result.fence==job.fence
            @test result.result.schema_version=="1.2" && result.result.inline_result==output.value
            @test result.result.warnings==output.warnings
            @test result.result.engine_version=="environment-sha256:"*job.fence.fingerprint
            @test result.result.cache_status=="bypass"
            @test RT.matching_result(job,result)
            @test !RT.matching_result(change_runtime_record(job;
                execution=change_runtime_record(job.execution;executor_generation=2)),result)
            for code in ("not_admitted","operation_rejected","canceled","deadline","executor_failed",
                    "execution_uncertain","result_payload_limit","artifact_unavailable")
                failure=RT.job_failure(job,started,code)
                @test RT.matching_result(job,failure)
                @test failure.result.failure.category==code && failure.result.inline_result===nothing
                @test !failure.result.failure.retryable
            end
            @test_throws ArgumentError RT.job_failure(job,started,"/private/exception details")
        finally
            close(agent)
        end
        @test agent.jobs.closed && agent.jobs.state==:stopped
        @test isempty(driver.processes) && agent.cleanup_complete
    end
end
