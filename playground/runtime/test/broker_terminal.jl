using Test, UUIDs, LineCableModelsRuntime
import NATS, JSON3
const RT=LineCableModelsRuntime
const P=RT.Protocol
directory,port=ARGS
certs=joinpath(directory,"certs")
endpoint(id)=BrokerEndpoint("tls://127.0.0.1:$port",joinpath(directory,id*".password");
    ca_file=joinpath(certs,"ca.pem"),certificate_file=joinpath(certs,"worker-cert.pem"),
    key_file=joinpath(certs,"worker-key.pem"),server_name="localhost")
function terminal_records(channel;timeout=5)
    result=RT.ControlEnvelope[]
    timedwait(timeout;pollint=0.01) do
        append!(result,poll_terminal!(channel));!isempty(result)
    end
    return result
end

# Initial fixture connection may spend a cold compilation interval inside a
# finite TLS handshake. Production services already retry connection attempts;
# mirror that bounded setup behavior without replaying any terminal operation.
function terminal_wire(args...;kwargs...)
    deadline=time()+15
    while true
        try
            return BrokerTerminal(args...;kwargs...)
        catch error
            error isa RT.BrokerUnavailable && time()<deadline || rethrow()
            sleep(0.1)
        end
    end
end

@testset "TLS terminal subjects enforce roles and exact assignment identities" begin
    coordinator=terminal_wire(endpoint("coordinator"),CoordinatorIdentity();worker_ids=("worker-a","worker-b"),capacity=2)
    a=b=control=worker_control=nothing
    try
        a=terminal_wire(endpoint("worker-a"),WorkerIdentity("worker-a"))
        b=terminal_wire(endpoint("worker-b"),WorkerIdentity("worker-b"))
        fence=P.AssignmentFence(string(uuid4()),string(uuid4()),"alice","terminal","worker-a",string(uuid4()),
            string(uuid4()),"julia-terminal","1.0.0",repeat("a",64),1)
        other=P.AssignmentFence(string(uuid4()),string(uuid4()),"bob","terminal","worker-b",string(uuid4()),
            fence.coordinator_id,"julia-terminal","1.0.0",repeat("a",64),1)
        command=P.TerminalCommand("2.0",string(uuid4()),fence,1,"open",nothing,string(uuid4()),0,0,100,30,UInt8[])
        @test isempty(coordinator.subscriptions) && isempty(a.subscriptions)
        @test_throws AccessDenied send_terminal!(coordinator,command)
        for (channel,target) in ((coordinator,fence),(coordinator,other),(a,fence),(b,other))
            @test watch_terminal!(channel,target) === nothing
            @test watch_terminal!(channel,target) === nothing
        end
        @test length(coordinator.subscriptions)==2 && length(a.subscriptions)==1
        @test_throws AccessDenied watch_terminal!(a,other)
        send_terminal!(coordinator,command)
        @test only(terminal_records(a)).record==command
        @test isempty(poll_terminal!(b))
        report=RT.terminal_rejection(command,"fixture_only")
        send_terminal!(a,report)
        @test only(terminal_records(coordinator)).record==report
        @test_throws AccessDenied send_terminal!(b,report)
        @test_throws AccessDenied send_terminal!(a,command)
        @test_throws AccessDenied send_terminal!(coordinator,report)
        NATS.publish(a.connection,P.terminal_subject(fence,:command),P.encode_message(command))
        NATS.ping(a.connection;measure=false)
        @test isempty(terminal_records(a;timeout=0.2)) # worker cannot self-open
        NATS.publish(a.connection,P.terminal_subject(other,:report),P.encode_message(report))
        NATS.ping(a.connection;measure=false)
        @test isempty(terminal_records(coordinator;timeout=0.2)) # cannot impersonate B
        foreign=P.AssignmentFence(fence.lease_id,fence.run_id,"bob",fence.role,fence.worker_id,fence.worker_boot,
            fence.coordinator_id,fence.profile_id,fence.profile_version,fence.fingerprint,fence.generation)
        forged=P.TerminalReport(report.protocol_version,report.request_id,foreign,report.revision,report.accepted,
            report.reason,report.session_id,report.phase,report.writer_connected,report.input_sequence,report.cursor,
            report.output_sequence,report.gap,report.bytes,report.failure,report.cleanup_pending)
        NATS.publish(a.connection,P.terminal_subject(fence,:report),P.encode_message(forged))
        @test isempty(terminal_records(coordinator;timeout=0.2))
        @test coordinator.rejected==1
        @test_throws AccessDenied watch_terminal!(coordinator,foreign)
        @test unwatch_terminal!(coordinator,foreign) === nothing && length(coordinator.subscriptions)==2
        for _ in 1:20;NATS.publish(a.connection,P.terminal_subject(fence,:report),"{}");end
        other_command=P.TerminalCommand("2.0",string(uuid4()),other,1,"open",nothing,string(uuid4()),0,0,80,24,UInt8[])
        send_terminal!(b,RT.terminal_rejection(other_command,"fixture_only"))
        observed=terminal_records(coordinator)
        @test any(record->record.worker_id=="worker-b",observed)
        @test coordinator.rejected>1
        control=BrokerControl(endpoint("coordinator"),CoordinatorIdentity();worker_ids=("worker-a",))
        worker_control=BrokerControl(endpoint("worker-a"),WorkerIdentity("worker-a"))
        probe=P.WorkerProbe("2.0","worker-a",fence.coordinator_id,string(uuid4()))
        send_control!(control,probe)
        @test timedwait(()->!isempty(poll_control!(worker_control)),5)==:ok
        @test control.connection !== coordinator.connection && worker_control.connection !== a.connection
        @test unwatch_terminal!(coordinator,fence) === nothing
        @test length(coordinator.subscriptions)==1
        @test_throws AccessDenied send_terminal!(coordinator,command)
        @test close(a) === nothing
        @test_throws RT.BrokerUnavailable send_terminal!(a,report)
    finally
        worker_control === nothing || close(worker_control)
        control === nothing || close(control)
        a === nothing || close(a);b === nothing || close(b);close(coordinator)
    end
end

include("agent_leases.jl")
include("terminal_resources.jl")
include("agent_terminals.jl")
include("terminal_sockets.jl")

@testset "TLS terminal channel carries the real private REPL without durable jobs" begin
    with_terminal_sessions(;capacity=1) do f
        resources=f.resources;service=f.agent.terminals;fence=only(f.fences)
        service.endpoint=endpoint("worker-a")
        wire=terminal_wire(endpoint("coordinator"),CoordinatorIdentity();worker_ids=("worker-a",))
        writer=string(uuid4());revision=Ref(0)
        function exchange(action;options...)
            revision[]+=1
            command=terminal_fixture_command(fence,revision[],action;options...)
            send_terminal!(wire,command)
            records=terminal_records(wire)
            @test length(records)==1
            report=only(records).record
            @test report.request_id==command.request_id && report.revision==command.revision && report.fence==fence
            return report,command
        end
        try
            recover_owned!(resources);start_terminals!(service);watch_terminal!(wire,fence)
            @test timedwait(()->service.connection !== nothing && length(service.connection.subscriptions)==1,10)==:ok
            opened,_=exchange("open";writer,columns=100,rows=30)
            @test opened.accepted
            id=opened.session_id
            terminal_wait_ready(resources,fence,id)
            input,command=exchange("input";session=id,writer,sequence=1,bytes=codeunits("println(\"TLS-REPL-λ\")\r"))
            @test input.accepted && input.input_sequence==1
            send_terminal!(wire,command)
            @test only(terminal_records(wire)).record==input
            terminal_expect(resources,fence,id,"TLS-REPL-λ")
            output,_=exchange("read";session=id)
            @test output.accepted && length(output.bytes)<=P.MAX_TERMINAL_CHUNK_BYTES
            @test output.cursor>0 && output.output_sequence>=output.cursor
            @test !occursin("lcm-terminal-ready",String(copy(output.bytes)))
            denied,_=exchange("input";session=id,writer=string(uuid4()),sequence=2,bytes=UInt8[3])
            @test !denied.accepted && denied.reason=="authority_rejected"
            @test terminal_status(resources,fence,id).input_sequence==1
            fresh,_=exchange("restart";session=id,writer,columns=80,rows=24)
            @test fresh.accepted && fresh.session_id!=id
            terminal_wait_ready(resources,fence,fresh.session_id)
            previous=service.connection
            close(previous) # exact test-owned channel loss, not broker-wide shutdown
            @test timedwait(()->!terminal_status(resources,fence,fresh.session_id).writer_connected,5)==:ok
            f.clock[]=3 # advance the fixture clock past the bounded reconnect delay
            @test timedwait(()->service.connection !== nothing && service.connection !== previous &&
                length(service.connection.subscriptions)==1,10)==:ok
            resumed,_=exchange("open";writer,columns=80,rows=24)
            @test resumed.accepted && resumed.session_id==fresh.session_id
            @test resumed.input_sequence==0
            @test !occursin(writer,repr(service))
        finally
            close(wire)
        end
    end
end

include("terminal_gateway_tls.jl")
