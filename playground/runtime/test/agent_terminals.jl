function terminal_fixture_command(fence,revision,action;session=nothing,writer=nothing,sequence=0,after=0,
        columns=0,rows=0,bytes=UInt8[],request_id=string(uuid4()))
    RT.Protocol.TerminalCommand("2.0",request_id,fence,revision,action,session,writer,sequence,after,columns,rows,collect(bytes))
end
function terminal_fixture_reply(service,command)
    report=receive_terminal_command!(service,command)
    if report === nothing
        lane=service.lanes[command.fence.lease_id]
        @test timedwait(()->istaskdone(lane.task),30)==:ok
        RT.completed_terminal_reports!(service)
        report=lane.report
    end
    @test report isa RT.Protocol.TerminalReport
    @test RT.Protocol.validate(report)===report
    return report
end

@testset "writer presence expires independently of broker and read traffic" begin
    with_terminal_sessions(;capacity=1,limits=TerminalSessionLimits(presence_seconds=1,disconnect_seconds=1)) do f
        recover_owned!(f.resources)
        fence=only(f.fences);writer=string(uuid4())
        id=open_terminal!(f.resources,fence,writer)
        terminal_wait_ready(f.resources,fence,id)
        session=only(values(f.resources.handles))
        f.clock[]=0.8
        command=terminal_fixture_command(fence,1,"keepalive";session=id,writer)
        @test terminal_fixture_reply(f.agent.terminals,command).accepted
        @test session.writer_seen_at==0.8 && session.activity_at==0
        f.clock[]=1.7
        @test terminal_fixture_reply(f.agent.terminals,command).accepted # cached ACK, not fresh presence
        @test session.writer_seen_at==0.8
        f.clock[]=1.9
        @test timedwait(()->session.disconnected_at !== nothing,2)==:ok
        @test isapprox(session.disconnected_at,1.8)
        @test terminal_status(f.resources,fence,id).writer_connected==false
        @test read_terminal(f.resources,fence,id,0) isa RT.TerminalRead
        rejected=terminal_fixture_reply(f.agent.terminals,terminal_fixture_command(fence,2,"keepalive";session=id,writer))
        @test !rejected.accepted && session.disconnected_at==1.8
        @test open_terminal!(f.resources,fence,writer)==id
        @test session.disconnected_at===nothing && session.writer_seen_at==1.9
        f.clock[]=4.1
        @test timedwait(()->istaskdone(session.task),3)==:ok
        @test session.failure==:disconnect_timeout && session.cleanup_complete
    end
end

@testset "terminal agent orders private actions and restart does not replay" begin
    with_terminal_sessions() do f
        r=f.resources;service=f.agent.terminals;a,b=f.fences
        @test service.resources === r && service.ledger === f.agent.ledger
        @test service.connection === nothing && service.task === nothing && isempty(service.lanes)
        recover_owned!(r)
        writer=string(uuid4())
        opening=terminal_fixture_command(a,1,"open";writer,columns=100,rows=30)
        first=terminal_fixture_reply(service,opening)
        @test first.accepted && first.session_id !== nothing
        id=first.session_id;terminal_wait_ready(r,a,id)
        @test terminal_fixture_reply(service,opening) === first
        @test length(r.handles)==1
        altered=change_runtime_record(opening;columns=80)
        @test !terminal_fixture_reply(service,altered).accepted
        input=terminal_fixture_command(a,2,"input";session=id,writer,sequence=1,bytes=codeunits("counter=41\r"))
        accepted=terminal_fixture_reply(service,input)
        @test accepted.accepted && accepted.input_sequence==1
        terminal_expect(r,a,id,"counter=41")
        @test terminal_fixture_reply(service,input) === accepted
        @test fieldnames(RT.AgentTerminalLane)==(:fence,:revision,:request_id,:digest,:task,:report)
        @test service.lanes[a.lease_id].task === nothing
        @test !occursin("counter=",repr(service.lanes[a.lease_id]))
        @test terminal_fixture_reply(service,opening).reason=="stale_revision"
        other=terminal_fixture_reply(service,terminal_fixture_command(b,1,"open";writer=string(uuid4()),columns=80,rows=24))
        terminal_wait_ready(r,b,other.session_id)
        old_process=r.handles[a.lease_id].process.process
        restart=terminal_fixture_command(a,3,"restart";session=id,writer,columns=120,rows=40)
        fresh=terminal_fixture_reply(service,restart)
        @test fresh.accepted && fresh.session_id!=id
        @test !process_running(old_process)
        terminal_wait_ready(r,a,fresh.session_id)
        @test terminal_fixture_reply(service,restart) === fresh
        @test r.handles[a.lease_id].id==fresh.session_id
        check=terminal_fixture_command(a,4,"input";session=fresh.session_id,writer,sequence=1,
            bytes=codeunits("println(\"FRESH:\",isdefined(Main,:counter))\r"))
        @test terminal_fixture_reply(service,check).accepted
        terminal_expect(r,a,fresh.session_id,"FRESH:false")
        @test !terminal_fixture_reply(service,terminal_fixture_command(a,5,"input";session=id,writer,sequence=2,bytes=UInt8[3])).accepted
        stopped=terminal_fixture_reply(service,terminal_fixture_command(a,6,"stop";session=fresh.session_id,writer))
        @test stopped.accepted && stopped.phase=="closing"
        @test terminal_status(r,b,other.session_id).phase==:ready
        @test timedwait(()->r.handles[a.lease_id].cleanup_complete,10)==:ok
        @test length(r.handles)==2 # explicit stop retains its ended identity
        @test_throws AccessDenied receive_terminal_command!(service,change_runtime_record(opening;
            fence=change_runtime_record(a;owner="foreign")))
    end
end

@testset "terminal channel loss starts writer grace independently of control" begin
    with_terminal_sessions(;capacity=1,limits=TerminalSessionLimits(disconnect_seconds=1)) do f
        r=f.resources;service=f.agent.terminals;fence=only(f.fences);writer=string(uuid4())
        recover_owned!(r)
        id=open_terminal!(r,fence,writer);terminal_wait_ready(r,fence,id)
        RT.terminal_service_state!(service,:online)
        RT.terminal_service_state!(service,:offline)
        @test !terminal_status(r,fence,id).writer_connected
        f.clock[]=0.5;RT.terminal_service_state!(service,:offline)
        @test r.handles[fence.lease_id].disconnected_at==0
        @test agent_lease_usable(f.agent.ledger,fence)
        f.clock[]=1.1
        @test timedwait(()->r.handles[fence.lease_id].cleanup_complete,10)==:ok
        @test r.handles[fence.lease_id].failure==:disconnect_timeout
    end
end
