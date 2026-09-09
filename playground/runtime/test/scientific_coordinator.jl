function ready_scientific_report(command; valid_for_ms=1000)
    RT.Protocol.ScientificReport("2.0",command.request_id,command.fence,command.revision,true,
        "accepted","idle","ready",string(uuid4()),1,nothing,repeat("b",64),1000,100,0,nothing,valid_for_ms)
end

@testset "queued cancellation safely retries lost replies without repeating preparation" begin
    coordinated_lease_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        science=ScientificCoordinator(BrokerEndpoint("tls://broker.invalid","/unused-password"),f.coordinator,["worker-a"])
        sent=RT.Protocol.ScientificCommand[]
        send=command->push!(sent,command)
        request,target=uuid4(),string(uuid4())
        try
            first=request_scientific!(science,f.alice,f.id,"cancel_job";target_id=target,request_id=request,send)
            @test request_scientific!(science,f.alice,f.id,"cancel_job";target_id=target,request_id=request,send)===first
            @test length(sent)==1
            # A newer status request can overtake a reply; a safe cancellation
            # retry still carries the same target under a newer command revision.
            rejected=change_runtime_record(RT.unavailable_scientific_report(first.command,"stale_command");accepted=false)
            @test accept_scientific_report!(science,"worker-a",rejected)
            second=request_scientific!(science,f.alice,f.id,"cancel_job";target_id=target,request_id=request,send)
            @test second.command.revision==2 && second.command.request_id==first.command.request_id
            @test length(science.lanes[string(f.id)].mutations)==1 && length(sent)==2
            accepted=change_runtime_record(RT.unavailable_scientific_report(second.command,"cancel_recorded");accepted=true)
            @test accept_scientific_report!(science,"worker-a",accepted)
            query=request_scientific!(science,f.alice,f.id,"status";send)
            @test request_scientific!(science,f.alice,f.id,"cancel_job";target_id=target,request_id=request,send)===second
            @test length(sent)==3 && query.command.action=="status"
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"cancel_job";
                target_id=target,request_id=UUID(query.command.request_id),send)
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"cancel_job";
                target_id=string(uuid4()),request_id=request,send)
        finally
            close(science)
        end
    end
end

@testset "remote preparation correlates exact reports without replaying mutations" begin
    coordinated_lease_fixture() do f
        science=ScientificCoordinator(BrokerEndpoint("tls://broker.invalid","/unused-password"),f.coordinator,["worker-a"])
        sent=RT.Protocol.ScientificCommand[]
        send=command->push!(sent,command)
        try
            @test isempty(science.lanes) && science.task===nothing && science.connection===nothing
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"prepare";send)
            grant=grant_assignment!(f.coordinator,f.alice,f.id)
            @test accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
            science.state=:online
            @test remote_scientific_status(science,f.alice,f.id).preparation=="unknown"
            @test_throws AccessDenied request_scientific!(science,f.bob,f.id,"status";send)
            @test_throws AccessDenied remote_scientific_status(science,f.bob,f.id)
            @test isempty(sent)
            request=uuid4()
            prepare=request_scientific!(science,f.alice,f.id,"prepare";parameters=Dict("value"=>2),request_id=request,send)
            @test prepare.command.action=="prepare" && prepare.command.revision==1
            @test request_scientific!(science,f.alice,f.id,"prepare";parameters=Dict("value"=>2),request_id=request,send)===prepare
            @test length(sent)==1
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"prepare";parameters=Dict("value"=>3),request_id=request,send)
            @test request_scientific!(science,f.alice,f.id,"status";send)===prepare
            @test_throws ArgumentError request_scientific!(science,f.alice,f.id,"status";parameters=Dict("value"=>2),send)
            @test !accept_scientific_report!(science,"worker-a",ready_scientific_report(prepare.command))
            acknowledgement=RT.unavailable_scientific_report(prepare.command,"preparation_not_admitted")
            @test accept_scientific_report!(science,"worker-a",acknowledgement)
            @test remote_scientific_status(science,f.alice,f.id).accepted===false
            query=request_scientific!(science,f.alice,f.id,"status";send)
            @test query.command.revision==2 && length(sent)==2
            report=ready_scientific_report(query.command)
            @test_throws AccessDenied accept_scientific_report!(science,"worker-b",report)
            @test !accept_scientific_report!(science,"worker-a",change_runtime_record(report;revision=1))
            @test !accept_scientific_report!(science,"worker-a",change_runtime_record(report;request_id=string(uuid4())))
            @test !accept_scientific_report!(science,"worker-a",change_runtime_record(report;
                fence=change_runtime_record(report.fence;generation=2)))
            f.clock[]=0.4
            @test accept_scientific_report!(science,"worker-a",report)
            @test query.ready_until==1.0 # full round trip is subtracted, not added
            @test remote_scientific_status(science,f.alice,f.id).preparation=="ready"
            next_query=request_scientific!(science,f.alice,f.id,"status";send)
            @test remote_scientific_status(science,f.alice,f.id).preparation=="ready"
            @test remote_scientific_status(science,f.alice,f.id).valid_for_ms<=600
            @test remote_scientific_status(science,f.alice,f.id).pending
            f.clock[]=0.9
            @test !accept_scientific_report!(science,"worker-a",report)
            @test query.ready_until==1.0
            f.clock[]=0.9999
            @test remote_scientific_status(science,f.alice,f.id).preparation=="unknown"
            @test remote_scientific_status(science,f.alice,f.id).valid_for_ms==0
            f.clock[]=1.01
            @test remote_scientific_status(science,f.alice,f.id).preparation=="unknown"
            @test remote_scientific_status(science,f.alice,f.id).preparation_key===nothing
            @test request_scientific!(science,f.alice,f.id,"prepare";parameters=Dict("value"=>2),request_id=request,send)===next_query
            @test length(sent)==3 # retry an old mutation after polling must remain inert
            query2=request_scientific!(science,f.alice,f.id,"status";send)
            @test accept_scientific_report!(science,"worker-a",ready_scientific_report(query2.command))
            science.state=:offline
            @test remote_scientific_status(science,f.alice,f.id).preparation=="unknown"
            f.clock[]=6
            @test_throws AccessDenied remote_scientific_status(science,f.alice,f.id)
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"prepare";send)
            tick_science!(science)
            @test isempty(science.lanes)
        finally
            close(science)
        end
        @test science.closed && science.state==:stopped
        @test close(science)===nothing
    end
end

@testset "uncertain scientific sends retain bounded retry identity and no ready report" begin
    coordinated_lease_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        science=ScientificCoordinator(BrokerEndpoint("tls://broker.invalid","/unused-password"),f.coordinator,["worker-a"])
        sent=RT.Protocol.ScientificCommand[]
        send=command->push!(sent,command)
        request=uuid4()
        try
            @test_throws RT.BrokerUnavailable request_scientific!(science,f.alice,f.id,"prepare";request_id=request,
                send=command->throw(RT.BrokerUnavailable()))
            @test !remote_scientific_status(science,f.alice,f.id).pending
            uncertain=science.lanes[string(f.id)].flight
            @test !accept_scientific_report!(science,"worker-a",RT.unavailable_scientific_report(uncertain.command,"late"))
            request_scientific!(science,f.alice,f.id,"prepare";request_id=request,send)
            @test isempty(sent)
            query=request_scientific!(science,f.alice,f.id,"status";send)
            @test query.command.revision==2 && length(sent)==1
            for _ in 1:255
                request_scientific!(science,f.alice,f.id,"cancel";target_id=string(request),send)
            end
            @test length(science.lanes[string(f.id)].mutations)==256
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"prepare";send)
            @test request_scientific!(science,f.alice,f.id,"prepare";request_id=request,send)!==nothing
            close(science)
            @test_throws AccessDenied request_scientific!(science,f.alice,f.id,"status";send)
        finally
            close(science)
        end
    end
end
