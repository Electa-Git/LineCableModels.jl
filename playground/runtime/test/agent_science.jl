function scientific_command(fence,revision,action;parameters=Dict{String,Any}(),target_id=nothing)
    RT.Protocol.ScientificCommand("2.0",string(uuid4()),fence,revision,action,parameters,target_id)
end

@testset "agent preparation channel is explicit, ordered and non-blocking" begin
    mktempdir() do directory
        (;driver,resources,agent,clock,fences)=scientific_fixture(directory;capacity=1)
        fence=only(fences)
        science=agent.science
        reports=RT.Protocol.ScientificReport[]
        emit=report->push!(reports,report)
        try
            @test science isa AgentScientificService
            @test science.ledger===agent.ledger && science.resources===resources
            @test_throws ArgumentError start_science!(science)
            recover_owned!(resources)
            status=scientific_command(fence,1,"status")
            receive_scientific!(science,status;emit)
            @test last(reports).preparation=="cold" && last(reports).valid_for_ms==0
            @test isempty(driver.processes) # status does not create or prepare an executor
            prepare=scientific_command(fence,2,"prepare";parameters=Dict{String,Any}("seconds"=>0.15))
            receive_scientific!(science,prepare;emit)
            @test last(reports).accepted && last(reports).preparation=="preparing"
            pending=resources.handles[fence.lease_id].task
            @test !istaskdone(pending)
            receive_scientific!(science,prepare;emit)
            @test resources.handles[fence.lease_id].task===pending
            @test last(reports).reason=="duplicate"
            @test_throws AccessDenied receive_scientific!(science,change_runtime_record(prepare;parameters=Dict{String,Any}());emit)
            receive_scientific!(science,status;emit)
            @test !last(reports).accepted && last(reports).reason=="stale_command"
            @test fetch(pending)["cache_status"]=="miss"
            # A duplicated acceptance cannot advertise readiness from a local cache.
            receive_scientific!(science,prepare;emit)
            @test last(reports).preparation=="unknown"
            query=scientific_command(fence,3,"status")
            receive_scientific!(science,query;emit)
            wait(science.lanes[fence.lease_id].task)
            ready=last(reports)
            @test ready.request_id==query.request_id && ready.preparation=="ready"
            @test 0<ready.valid_for_ms<=5000 && ready.preparation_key!==nothing
            @test ready.executor_id==scientific_status(resources,fence).executor_id
            receive_scientific!(science,query;emit)
            @test last(reports).preparation=="unknown" && last(reports).valid_for_ms==0
            slow=scientific_command(fence,4,"prepare";parameters=Dict{String,Any}("seconds"=>10.0))
            receive_scientific!(science,slow;emit)
            pending=resources.handles[fence.lease_id].task
            rejected=scientific_command(fence,5,"prepare";parameters=Dict{String,Any}("seconds"=>11.0))
            receive_scientific!(science,rejected;emit)
            @test !last(reports).accepted && last(reports).reason=="preparation_not_admitted"
            receive_scientific!(science,rejected;emit)
            @test !last(reports).accepted # duplicate rejection must not turn into acceptance
            query=scientific_command(fence,6,"status")
            receive_scientific!(science,query;emit)
            @test last(reports).preparation=="preparing"
            @test last(reports).current_request_id==slow.request_id
            tick_agent!(agent)
            @test (@elapsed tick_agent!(agent))<0.1
            receive_scientific!(science,scientific_command(fence,7,"cancel";target_id=string(uuid4()));emit)
            @test last(reports).reason=="not_pending" && !istaskdone(pending)
            cancel=scientific_command(fence,8,"cancel";target_id=slow.request_id)
            receive_scientific!(science,cancel;emit)
            @test last(reports).reason=="cancel_requested"
            @test_throws TaskFailedException fetch(pending)
            receive_scientific!(science,scientific_command(fence,9,"status");emit)
            @test last(reports).preparation!="ready" && last(reports).failure=="canceled"
            receive_scientific!(science,cancel;emit)
            @test last(reports).reason=="stale_command"
            @test_throws AccessDenied receive_scientific!(science,
                scientific_command(change_runtime_record(fence;owner="foreign"),10,"status");emit)
            clock[]=6
            @test_throws AccessDenied receive_scientific!(science,scientific_command(fence,10,"status");emit)
            tick_science!(science)
            @test isempty(science.lanes)
        finally
            close(agent)
        end
        @test science.closed && science.state==:stopped && isempty(science.lanes)
        @test isempty(driver.processes) && agent.cleanup_complete
    end
end
