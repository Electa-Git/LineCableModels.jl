using Test, UUIDs, LineCableModelsRuntime
const SessionRT=LineCableModelsRuntime

# This required-interface fixture owns finite native REPL children exclusively
# for byte/session tests. Production terminal admission remains container-only.
mutable struct TerminalDriverFixture <: AbstractTerminalDriver
    profiles::ProfileRegistry
    handles::Dict{String,Tuple{SessionRT.Protocol.AssignmentFence,SessionRT.TerminalProcess}}
    recovered::Int
    verified::Int
    closed::Bool
    allow_release::Bool
    emit_marker::Bool
    gate::Union{Nothing,Channel{Nothing}}
end

@testset "startup markers are bounded, exact and omit only the startup prefix" begin
    process=SessionRT.TerminalProcess(;limits=SessionRT.TerminalIOLimits(output_bytes=8192))
    marker=collect(codeunits("\x1elcm-terminal-ready:"*string(uuid4())*"\x1f"))
    try
        @test SessionRT.terminal_startup_offset(process,marker) === nothing
        prefix=collect(codeunits("Julia startup\r\n"))
        SessionRT.append_terminal!(process.output,[prefix;marker[1:end-1]])
        @test SessionRT.terminal_startup_offset(process,marker) === nothing
        SessionRT.append_terminal!(process.output,[last(marker);collect(codeunits("julia> "))])
        @test SessionRT.terminal_startup_offset(process,marker)==length(prefix)+length(marker)
        @test SessionRT.terminal_startup_offset(process,[marker;0x01]) === nothing
        @test_throws SessionRT.TerminalFailure SessionRT.terminal_startup_offset(process,zeros(UInt8,129))
        SessionRT.append_terminal!(process.output,zeros(UInt8,8193))
        @test_throws SessionRT.TerminalFailure SessionRT.terminal_startup_offset(process,marker)
    finally
        close(process)
    end
end
SessionRT.installed_profiles(d::TerminalDriverFixture)=d.profiles
SessionRT.recover_owned!(d::TerminalDriverFixture)=(d.recovered+=1;nothing)
SessionRT.verify_terminal!(d::TerminalDriverFixture,::ProfileDefinition,::SessionRT.Protocol.AssignmentFence)=(d.verified+=1;nothing)
function SessionRT.terminal_for!(d::TerminalDriverFixture,::ProfileDefinition,fence::SessionRT.Protocol.AssignmentFence)
    found=get!(d.handles,fence.lease_id) do
        (fence,SessionRT.TerminalProcess(;limits=SessionRT.TerminalIOLimits(
            chunk_bytes=128,input_bytes=1024,output_bytes=2048,lifetime_seconds=45)))
    end
    found[1]==fence || error("fixture fence differs")
    return found[2]
end
function SessionRT.terminal_ready_marker(d::TerminalDriverFixture,p::ProfileDefinition,f::SessionRT.Protocol.AssignmentFence)
    process=SessionRT.terminal_for!(d,p,f)
    return collect(codeunits("\x1elcm-terminal-ready:"*string(process.id)*"\x1f"))
end
function SessionRT.start_owned_terminal!(d::TerminalDriverFixture,p::ProfileDefinition,f::SessionRT.Protocol.AssignmentFence,cols::Int,rows::Int)
    process=SessionRT.terminal_for!(d,p,f)
    process.started_at==0 || return process
    d.gate === nothing || take!(d.gate)
    julia=joinpath(Sys.BINDIR,Base.julia_exename())
    child=joinpath(@__DIR__,"terminal_session_child.jl")
    project=dirname(@__DIR__)
    command=setenv(`setpriv --pdeathsig KILL $julia --startup-file=no --history-file=no --compiled-modules=existing --project=$project --load=$child --interactive --color=yes`,
        ["PATH"=>ENV["PATH"],"JULIA_DEPOT_PATH"=>join(DEPOT_PATH,':'),"JULIA_LOAD_PATH"=>"@:@stdlib",
         "JULIA_NUM_THREADS"=>"1","OPENBLAS_NUM_THREADS"=>"1","TERM"=>"xterm-256color",
         "LCM_TERMINAL_READY"=>string(process.id),"LCM_TEST_READY"=>d.emit_marker ? "yes" : "no"])
    SessionRT.start_terminal!(process,command;columns=cols,rows)
    return process
end
function SessionRT.release_owned!(d::TerminalDriverFixture,fence::SessionRT.Protocol.AssignmentFence)
    d.allow_release || return false
    found=get(d.handles,fence.lease_id,nothing)
    if found !== nothing
        found[1]==fence || error("fixture cleanup fence differs")
        close(found[2]);delete!(d.handles,fence.lease_id)
    end
    return true
end
function Base.close(d::TerminalDriverFixture)
    for (fence,_) in collect(values(d.handles))
        release_owned!(d,fence) || error("fixture cleanup refused")
    end
    d.closed=true
    return nothing
end

function with_terminal_sessions(action;capacity=2,limits=TerminalSessionLimits(),emit_marker=true)
    mktempdir() do directory
        profile=ProfileDefinition("terminal","registry.invalid/lcm@sha256:"*repeat("a",64),repeat("a",64);
            kind=:terminal,isolation=:container)
        profiles=ProfileRegistry();register!(profiles,profile)
        driver=TerminalDriverFixture(profiles,Dict{String,Tuple{SessionRT.Protocol.AssignmentFence,SessionRT.TerminalProcess}}(),
            0,0,false,true,emit_marker,nothing)
        resources=TerminalResources(driver;limits)
        config=AgentConfig("worker-a",BrokerEndpoint("tls://broker.invalid","/unused-password"),profiles,
            joinpath(directory,"owned");capacity)
        clock=Ref(0.0)
        agent=AgentService(config,resources;clock=()->clock[])
        coordinator=string(uuid4())
        receive_probe!(agent.ledger,SessionRT.Protocol.WorkerProbe("2.0","worker-a",coordinator,string(uuid4())))
        fences=[SessionRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"owner-$i","terminal","worker-a",
            agent.ledger.boot_id,coordinator,profile.id,"1.0.0",profile.fingerprint,1) for i in 1:capacity]
        for fence in fences;handle_lease_control!(agent.ledger,lease_command(fence;duration_ms=60_000));end
        try
            action((;driver,resources,agent,clock,fences))
        finally
            driver.allow_release=true
            driver.gate === nothing || isready(driver.gate) || put!(driver.gate,nothing)
            close(agent)
        end
        @test driver.closed && isempty(driver.handles) && isempty(resources.handles)
    end
end

function terminal_wait_ready(r,f,id)
    @test timedwait(()->terminal_status(r,f,id).phase!=:starting,30;pollint=0.02)==:ok
    @test terminal_status(r,f,id).phase==:ready
end
function terminal_contains(r,f,id,text)
    cursor=0;bytes=UInt8[]
    while true
        result=read_terminal(r,f,id,cursor);append!(bytes,result.bytes);cursor=result.cursor
        cursor==result.sequence && break
    end
    return occursin(text,String(bytes))
end
function terminal_expect(r,f,id,text)
    @test timedwait(()->terminal_contains(r,f,id,text),10;pollint=0.02)==:ok
end

struct MissingTerminalDriver <: AbstractTerminalDriver end
@testset "terminal session admission, sole writer and actual REPL are private" begin
    @test_throws ArgumentError TerminalResources(MissingTerminalDriver())
    @test_throws ArgumentError TerminalSessionLimits(idle_seconds=0)
    with_terminal_sessions() do f
        r=f.resources;a,b=f.fences;wa,wb=string(uuid4()),string(uuid4())
        @test isempty(f.driver.handles) && f.driver.verified==0
        @test_throws AccessDenied open_terminal!(r,a,wa)
        @test recover_owned!(r) === nothing
        @test recover_owned!(r) === nothing && f.driver.recovered==1
        @test_throws ArgumentError bind_agent!(r,AgentLeaseLedger("worker-a",f.driver.profiles))
        f.driver.gate=Channel{Nothing}(1)
        id=open_terminal!(r,a,wa)
        @test open_terminal!(r,a,wa)==id
        @test_throws AccessDenied open_terminal!(r,a,wb)
        @test_throws SessionRT.TerminalFailure write_terminal!(r,a,id,wa,1,UInt8[0x03])
        @test terminal_status(r,a,id).phase==:starting
        put!(f.driver.gate,nothing)
        terminal_wait_ready(r,a,id)
        f.driver.gate=nothing
        @test f.driver.verified==1
        @test !terminal_contains(r,a,id,"lcm-terminal-ready")
        @test_throws AccessDenied terminal_status(r,b,id)
        @test_throws AccessDenied terminal_status(r,change_runtime_record(a;owner="foreign"),id)
        @test_throws AccessDenied write_terminal!(r,a,id,wb,1,UInt8[0x03])
        @test write_terminal!(r,a,id,wa,1,codeunits("counter=1\r"))==1
        terminal_expect(r,a,id,"counter=1")
        @test write_terminal!(r,a,id,wa,2,codeunits("counter+=1\r"))==2
        @test write_terminal!(r,a,id,wa,2,codeunits("counter+=1\r"))==2
        @test_throws AccessDenied write_terminal!(r,a,id,wa,2,codeunits("counter+=2\r"))
        @test_throws AccessDenied write_terminal!(r,a,id,wa,4,UInt8[0x03])
        @test_throws AccessDenied write_terminal!(r,a,id,wa,1,codeunits("counter=1\r"))
        @test write_terminal!(r,a,id,wa,3,codeunits("println(\"COUNT:\",counter)\r"))==3
        terminal_expect(r,a,id,"COUNT:2")
        @test resize_terminal!(r,a,id,wa,120,40) === nothing
        @test_throws ArgumentError resize_terminal!(r,a,id,wa,0,40)
        @test_throws ArgumentError read_terminal(r,a,id,-1)
        @test_throws ArgumentError write_terminal!(r,a,id,wa,4,zeros(UInt8,129))
        @test terminal_status(r,a,id).input_sequence==3
        other=open_terminal!(r,b,wb);terminal_wait_ready(r,b,other)
        @test other!=id
        @test write_terminal!(r,b,other,wb,1,codeunits("println(\"SEPARATE:\",isdefined(Main,:counter))\r"))==1
        terminal_expect(r,b,other,"SEPARATE:false")
        @test disconnect_terminal!(r,a,id,wa) === nothing
        f.clock[]=1
        @test disconnect_terminal!(r,a,id,wa) === nothing
        @test r.handles[a.lease_id].disconnected_at==0
        @test_throws SessionRT.TerminalFailure write_terminal!(r,a,id,wa,4,UInt8[0x03])
        @test_throws AccessDenied open_terminal!(r,a,wb)
        @test open_terminal!(r,a,wa)==id
        @test terminal_status(r,a,id).writer_connected
        @test write_terminal!(r,a,id,wa,4,codeunits("sleep(10)\r"))==4
        sleep(0.2)
        @test write_terminal!(r,a,id,wa,5,UInt8[0x03])==5
        terminal_expect(r,a,id,"InterruptException")
        @test terminal_status(r,b,other).phase==:ready
        @test !occursin(wa,repr(r.handles[a.lease_id]))
        @test !occursin("counter",repr(r))
        @test write_terminal!(r,a,id,wa,6,codeunits("exit()\r"))==6
        @test timedwait(()->terminal_status(r,a,id).phase==:exited,10)==:ok
        @test_throws SessionRT.TerminalFailure open_terminal!(r,a,wa)
        @test terminal_status(r,b,other).phase==:ready
        @test release_owned!(r,a)
    end
end

@testset "terminal deadlines and lease loss retire the exact process" begin
    for reason in (:disconnect_timeout,:idle_timeout,:lease_lost)
        with_terminal_sessions(;capacity=1,limits=TerminalSessionLimits(disconnect_seconds=1,idle_seconds=2)) do f
            r=f.resources;fence=only(f.fences);writer=string(uuid4())
            recover_owned!(r);id=open_terminal!(r,fence,writer);terminal_wait_ready(r,fence,id)
            handle=r.handles[fence.lease_id];process=handle.process.process
            reason==:disconnect_timeout && disconnect_terminal!(r,fence,id,writer)
            f.clock[]=reason==:disconnect_timeout ? 1.1 : reason==:idle_timeout ? 2.1 : 11.0
            @test timedwait(()->istaskdone(handle.task),10)==:ok
            @test handle.failure==reason && handle.closing && handle.cleanup_complete
            @test !process_running(process) && isempty(f.driver.handles)
            @test haskey(r.handles,fence.lease_id) # no implicit replacement in this lease
            @test release_owned!(r,fence)
        end
    end
end

@testset "startup evidence is required and unresolved cleanup keeps ownership" begin
    with_terminal_sessions(;capacity=1,limits=TerminalSessionLimits(startup_seconds=2),emit_marker=false) do f
        r=f.resources;fence=only(f.fences);writer=string(uuid4())
        recover_owned!(r);id=open_terminal!(r,fence,writer)
        handle=r.handles[fence.lease_id]
        @test timedwait(()->istaskdone(handle.task),15)==:ok
        @test handle.phase==:failed && handle.failure==:startup_timeout && handle.cleanup_complete
        @test_throws SessionRT.TerminalFailure read_terminal(r,fence,id,0)
        @test_throws SessionRT.TerminalFailure open_terminal!(r,fence,writer)
    end
    with_terminal_sessions(;capacity=1,emit_marker=false,limits=TerminalSessionLimits(disconnect_seconds=1)) do f
        r=f.resources;fence=only(f.fences);writer=string(uuid4())
        recover_owned!(r);id=open_terminal!(r,fence,writer)
        handle=r.handles[fence.lease_id]
        @test timedwait(()->handle.process !== nothing && handle.process.started_at>0,10)==:ok
        disconnect_terminal!(r,fence,id,writer);f.clock[]=1.1
        @test timedwait(()->istaskdone(handle.task),10)==:ok
        @test handle.failure==:disconnect_timeout && handle.cleanup_complete
        @test !process_running(handle.process.process)
    end
    with_terminal_sessions(;capacity=1) do f
        r=f.resources;fence=only(f.fences)
        recover_owned!(r);id=open_terminal!(r,fence,string(uuid4()));terminal_wait_ready(r,fence,id)
        handle=r.handles[fence.lease_id]
        f.driver.allow_release=false
        @test !release_owned!(r,fence)
        @test !process_running(handle.process.process)
        @test haskey(r.handles,fence.lease_id) && haskey(f.driver.handles,fence.lease_id)
        @test terminal_status(r,fence,id).cleanup_pending
        @test_throws ArgumentError close(r)
        f.driver.allow_release=true
        @test close(r) === nothing
        @test handle.cleanup_complete && isempty(r.handles) && isempty(f.driver.handles)
    end
end
