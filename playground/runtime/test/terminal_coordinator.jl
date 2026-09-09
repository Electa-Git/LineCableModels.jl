function coordinated_terminal_fixture(action)
    coordinated_lease_fixture(;profile_kind=:terminal,duration_ms=60000) do f
        command=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,command))
        terminals=TerminalCoordinator(BrokerEndpoint("tls://broker.invalid","/unused-password"),f.coordinator,["worker-a"])
        sent=RT.Protocol.TerminalCommand[]
        try
            action((;f...,terminals,sent,send=command->push!(sent,command)))
        finally
            close(terminals)
        end
    end
end

function terminal_coordinator_report(command;session_id=command.session_id,options...)
    base=RT.Protocol.TerminalReport("2.0",command.request_id,command.fence,command.revision,true,
        "accepted",session_id,"ready",true,command.input_sequence,0,0,false,UInt8[],nothing,false)
    change_runtime_record(base;options...)
end
terminal_sent(flight)=(wait(flight.sender);flight)

@testset "terminal handoff waits only for exact closing ownership" begin
    coordinated_terminal_fixture() do f
        t=f.terminals
        prior=RT.claim_terminal_attachment!(t,f.alice,f.id)
        @test_throws AccessDenied RT.claim_terminal_attachment!(t,f.alice,f.id)
        prior.closed=true
        @test_throws AccessDenied RT.claim_terminal_attachment!(t,f.bob,f.id)
        successor=@async RT.claim_terminal_attachment!(t,f.alice,f.id)
        sleep(0.1)
        @test !istaskdone(successor)
        @test t.attachments[prior.fence.lease_id]===prior
        lock(()->delete!(t.attachments,prior.fence.lease_id),t.lock)
        @test timedwait(()->istaskdone(successor),1)==:ok
        next=fetch(successor)
        @test next!==prior && next.fence==prior.fence
        lock(()->delete!(t.attachments,next.fence.lease_id),t.lock)
    end
    coordinated_terminal_fixture() do f
        t=f.terminals
        prior=RT.claim_terminal_attachment!(t,f.alice,f.id);prior.closed=true
        successor=@async RT.claim_terminal_attachment!(t,f.alice,f.id)
        sleep(0.1)
        transition_run!(f.store,f.alice,f.a_run.id,:stopped)
        @test timedwait(()->istaskdone(successor),1)==:ok
        @test_throws TaskFailedException fetch(successor)
        @test t.attachments[prior.fence.lease_id]===prior
        lock(()->delete!(t.attachments,prior.fence.lease_id),t.lock)
    end
    coordinated_terminal_fixture() do f
        t=f.terminals
        prior=RT.claim_terminal_attachment!(t,f.alice,f.id);prior.closed=true
        successor=@async RT.claim_terminal_attachment!(t,f.alice,f.id)
        @test timedwait(()->istaskdone(successor),7)==:ok
        @test_throws TaskFailedException fetch(successor)
        @test t.attachments[prior.fence.lease_id]===prior
        lock(()->delete!(t.attachments,prior.fence.lease_id),t.lock)
    end
end

@testset "terminal coordinator binds private access to the actual owner, not inventory administrators" begin
    coordinated_terminal_fixture() do f
        t=f.terminals
        @test t.connection===nothing && t.task===nothing && isempty(t.lanes)
        writer=string(uuid4())
        for principal in (f.bob,f.admin)
            @test_throws AccessDenied request_terminal!(t,principal,f.id,"open";writer_id=writer,columns=100,rows=30,send=f.send)
        end
        @test isempty(f.sent) && isempty(t.lanes)
        request=uuid4()
        flight=terminal_sent(request_terminal!(t,f.alice,f.id,"open";writer_id=writer,columns=100,rows=30,request_id=request,send=f.send))
        command=only(f.sent)
        @test command.revision==1 && command.fence==f.lease.fence
        @test request_terminal!(t,f.alice,f.id,"open";writer_id=writer,columns=100,rows=30,request_id=request,send=f.send)===flight
        @test length(f.sent)==1
        @test_throws AccessDenied request_terminal!(t,f.alice,f.id,"open";writer_id=writer,columns=80,rows=30,request_id=request,send=f.send)
        @test_throws AccessDenied request_terminal!(t,f.alice,f.id,"status";session_id=string(uuid4()),send=f.send)
        session=string(uuid4());report=terminal_coordinator_report(command;session_id=session)
        @test_throws AccessDenied accept_terminal_report!(t,"worker-b",report)
        @test !accept_terminal_report!(t,"worker-a",change_runtime_record(report;request_id=string(uuid4())))
        @test !accept_terminal_report!(t,"worker-a",change_runtime_record(report;revision=2))
        @test !accept_terminal_report!(t,"worker-a",change_runtime_record(report;fence=change_runtime_record(report.fence;owner="bob")))
        @test accept_terminal_report!(t,"worker-a",report)
        @test !accept_terminal_report!(t,"worker-a",report)
        @test flight.command===nothing && flight.report===report
        @test terminal_flight(t,f.alice,f.id,request)===flight
        @test_throws AccessDenied terminal_flight(t,f.bob,f.id,request)
        @test_throws AccessDenied terminal_flight(t,f.admin,f.id,request)
        @test !occursin(writer,repr(t)) && !occursin(writer,repr(flight))
        transition_run!(f.store,f.alice,f.a_run.id,:stopped)
        @test_throws AccessDenied terminal_flight(t,f.alice,f.id,request)
        RT.tick_terminals!(t)
        @test isempty(t.lanes)
        close(t)
        @test t.closed && t.state==:stopped && close(t)===nothing
    end
    coordinated_lease_fixture() do f
        grant=grant_assignment!(f.coordinator,f.alice,f.id)
        accept_lease_ack!(f.coordinator,"worker-a",handle_lease_control!(f.agent,grant))
        t=TerminalCoordinator(BrokerEndpoint("tls://broker.invalid","/unused-password"),f.coordinator,["worker-a"])
        try
            @test_throws AccessDenied request_terminal!(t,f.alice,f.id,"open";writer_id=string(uuid4()),columns=80,rows=24)
            @test isempty(t.lanes)
        finally;close(t);end
    end
end

@testset "terminal requests discard private bytes and require exact explicit mutation retries" begin
    coordinated_terminal_fixture() do f
        t=f.terminals;writer=string(uuid4());session=string(uuid4());request=uuid4()
        input=Vector{UInt8}(codeunits("private_value = 41\r"))
        failed=command->throw(RT.BrokerUnavailable())
        first=terminal_sent(request_terminal!(t,f.alice,f.id,"input";session_id=session,writer_id=writer,
            input_sequence=1,bytes=input,request_id=request,send=failed))
        @test first.deadline==0 && first.report===nothing
        command=first.command
        RT.tick_terminals!(t)
        @test first.command===nothing && first.sender===nothing
        @test !accept_terminal_report!(t,"worker-a",terminal_coordinator_report(command))
        @test request_terminal!(t,f.alice,f.id,"input";session_id=session,writer_id=writer,
            input_sequence=1,bytes=input,request_id=request,send=f.send)===first
        @test isempty(f.sent) # no implicit mutation replay on reconnect or ordinary retry
        @test_throws AccessDenied request_terminal!(t,f.alice,f.id,"input";session_id=session,writer_id=writer,
            input_sequence=1,bytes=UInt8[1],request_id=request,retry=true,send=f.send)
        retried=terminal_sent(request_terminal!(t,f.alice,f.id,"input";session_id=session,writer_id=writer,
            input_sequence=1,bytes=input,request_id=request,retry=true,send=f.send))
        @test retried===first && only(f.sent)==command
        bad=terminal_coordinator_report(command;input_sequence=2)
        @test !accept_terminal_report!(t,"worker-a",bad)
        @test !accept_terminal_report!(t,"worker-a",terminal_coordinator_report(command;session_id=string(uuid4())))
        @test accept_terminal_report!(t,"worker-a",terminal_coordinator_report(command))
        @test first.command===nothing
        @test !occursin("private_value",repr(first))
        query=terminal_sent(request_terminal!(t,f.alice,f.id,"read";session_id=session,after=12,send=f.send))
        read_command=last(f.sent)
        @test query.revision==2
        @test_throws AccessDenied request_terminal!(t,f.alice,f.id,"input";session_id=session,writer_id=writer,
            input_sequence=1,bytes=input,request_id=request,retry=true,send=f.send)
        @test !accept_terminal_report!(t,"worker-a",terminal_coordinator_report(read_command;bytes=UInt8[65],cursor=14,output_sequence=14))
        @test accept_terminal_report!(t,"worker-a",terminal_coordinator_report(read_command;bytes=UInt8[65],cursor=13,output_sequence=14))
        @test query.report.bytes==UInt8[65]
        @test_throws AccessDenied terminal_flight(t,f.alice,f.id,request)
        next=terminal_sent(request_terminal!(t,f.alice,f.id,"read";session_id=session,after=13,send=f.send))
        @test accept_terminal_report!(t,"worker-a",terminal_coordinator_report(last(f.sent);bytes=UInt8[66],gap=true,cursor=100,output_sequence=100))
        @test next.report.gap
        restart=terminal_sent(request_terminal!(t,f.alice,f.id,"restart";session_id=session,writer_id=writer,columns=80,rows=24,send=f.send))
        @test !accept_terminal_report!(t,"worker-a",terminal_coordinator_report(last(f.sent)))
        @test accept_terminal_report!(t,"worker-a",terminal_coordinator_report(last(f.sent);session_id=string(uuid4())))
        @test restart.report.session_id!=session
    end
end

@testset "terminal deadlines fence late reports and shutdown rejects scheduled sends" begin
    coordinated_terminal_fixture() do f
        t=f.terminals;session=string(uuid4())
        query=terminal_sent(request_terminal!(t,f.alice,f.id,"status";session_id=session,send=f.send))
        command=only(f.sent)
        f.clock[]=6
        f.announce("worker-a";sequence=2,boot=f.agent.boot_id)
        RT.tick_terminals!(t)
        @test query.command===nothing
        @test !accept_terminal_report!(t,"worker-a",terminal_coordinator_report(command))
        next=terminal_sent(request_terminal!(t,f.alice,f.id,"status";session_id=session,send=f.send))
        @test next.revision==2
        @test !accept_terminal_report!(t,"worker-a",terminal_coordinator_report(command))
        @test accept_terminal_report!(t,"worker-a",terminal_coordinator_report(last(f.sent)))
        @test_throws ArgumentError request_terminal!(t,f.alice,f.id,"status";session_id=session,columns=20,send=f.send)
        @test length(f.sent)==2
        pending=request_terminal!(t,f.alice,f.id,"status";session_id=session,send=f.send)
        close(t) # queued @async sender must see closure before publishing
        @test length(f.sent)==2 && istaskdone(pending.sender)
        @test_throws AccessDenied request_terminal!(t,f.alice,f.id,"status";session_id=session,send=f.send)
        @test isempty(t.lanes)
    end
end
