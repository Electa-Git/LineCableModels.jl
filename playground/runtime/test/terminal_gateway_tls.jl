import HTTP, JSON3
include("terminal_live_browser.jl")

# Real gateway, control/terminal coordinators, TLS broker and host agent. Only
# the physical driver is an explicit finite native REPL fixture, not a sandbox.
@testset "protected browser terminal relays a live leased REPL over TLS" begin
    mktempdir() do dir
        profiles=ProfileRegistry()
        register!(profiles,ProfileDefinition("terminal","fixture@sha256:"*repeat("a",64),repeat("a",64);
            kind=:terminal,isolation=:container))
        applications=ApplicationRegistry()
        definition=ApplicationDefinition("terminal-fixture","Terminal fixture",:workbench,"/terminal-fixture";
            requirements=(RuntimeRequirement("terminal",("terminal",)),))
        register!(applications,definition)
        driver=TerminalDriverFixture(profiles,Dict{String,Tuple{P.AssignmentFence,RT.TerminalProcess}}(),
            0,0,false,true,true,nothing)
        resources=TerminalResources(driver)
        agent=AgentService(AgentConfig("worker-a",endpoint("worker-a"),profiles,joinpath(dir,"owned");capacity=2),resources)
        store=RuntimeStore(joinpath(dir,"runtime.sqlite"))
        trust=WorkerTrust("worker-a","worker-a",("terminal",);capacity=2)
        service=ControlService(ControlConfig(endpoint("coordinator"),profiles,[trust]),store,applications)
        supervisor=UIHostSupervisor(store,applications,joinpath(dir,"hosts"))
        key=repeat("terminal-fixture-proxy-",3)
        policy=ProxyIdentity("https://lcm.test",["127.0.0.1"],key;administrators=["operator"])
        server=start_gateway(supervisor,policy;control=service)
        base="http://127.0.0.1:$(HTTP.port(server))"
        headers(user;origin="https://lcm.test")=["X-LCM-Proxy-Key"=>key,"X-LCM-Principal"=>user,"Origin"=>origin]
        alice,bob,operator=Principal("alice"),Principal("bob"),Principal("operator";administrator=true)
        sockets=Any[]
        socket_locks=IdDict{Any,ReentrantLock}()
        heartbeat=nothing;heartbeat_stop=Ref(false)
        function attach(id,user,writer)
            socket=HTTP.WebSockets.open(replace(base,"http:"=>"ws:")*"/runtime/api/assignments/$id/terminal";
                headers=headers(user),proxy=nothing,cookies=false,request_timeout=5,maxframesize=65536,read_idle_timeout=10)
            push!(sockets,socket)
            socket_locks[socket]=ReentrantLock()
            @test JSON3.read(HTTP.WebSockets.receive(socket)).kind=="hello"
            HTTP.WebSockets.send(socket,JSON3.write((action="attach",writer_id=writer)))
            @test JSON3.read(HTTP.WebSockets.receive(socket)).kind=="attached"
            return socket
        end
        function action(socket,name;session=nothing,sequence=0,after=0,columns=0,rows=0,bytes=UInt8[],request_id=string(uuid4()),retry=false)
            lock(socket_locks[socket]) do
                HTTP.WebSockets.send(socket,JSON3.write((action=name,request_id,session_id=session,
                    input_sequence=sequence,after,columns,rows,bytes=collect(bytes),retry)))
                reply=JSON3.read(HTTP.WebSockets.receive(socket))
                @test reply.kind=="report" && reply.request_id==request_id
                return reply
            end
        end
        function ready(socket,id)
            latest=Ref{Any}(nothing)
            @test timedwait(15;pollint=0.05) do
                latest[]=action(socket,"status";session=id)
                latest[].phase!="starting"
            end==:ok
            @test latest[].accepted && latest[].phase=="ready"
        end
        function expect_output(socket,id,needle;after=0)
            bytes=UInt8[];cursor=after
            @test timedwait(10;pollint=0.05) do
                report=action(socket,"read";session=id,after=cursor)
                @test report.accepted
                cursor=report.cursor;append!(bytes,UInt8.(report.bytes))
                occursin(needle,String(copy(bytes)))
            end==:ok
            return cursor,String(bytes)
        end
        function admit(principal,run)
            # Cold fixture-only HTTP setup can age the last presence report.
            # Wait for a genuinely fresh challenged report, never extend its
            # validity or admit against the stale report still in the table.
            @test timedwait(15;pollint=0.05) do
                lock(service.inventory.lock) do
                    RT.presence_state(service.inventory,
                        get(service.inventory.presence,"worker-a",nothing))==:online
                end
            end==:ok
            # Cold compilation can consume this combined-process fixture's
            # two-second ACK budget. Exercise confirmed cleanup and a fresh
            # reservation, never extend or revive the expired grant.
            for attempt in 1:3
                lease=reserve_assignment!(service.coordinator.assignments,principal,run.id,"terminal","terminal";
                    placement=PinnedPlacement("worker-a"))
                id=UUID(lease.fence.lease_id)
                grant_assignment!(service.coordinator,principal,id)
                if timedwait(()->assignment_usable(service.coordinator,principal,id),5)==:ok
                    return lease
                end
                println("Terminal fixture setup: grant attempt ",attempt," expired; awaiting exact cleanup before fresh reservation.")
                release_assignment!(service.coordinator,principal,id)
                @test timedwait(()->get_assignment(store,principal,id).state in (:released,:expired),10)==:ok
                @test !assignment_usable(service.coordinator,principal,id)
            end
            error("terminal fixture could not establish live lease authority")
        end
        try
            # This combined-process fixture shares Julia's compiler with both
            # peers. Warm gateway dispatch before starting wall-clock leases;
            # real deployments have separately initialized gateway/agent hosts.
            @test HTTP.get(base*"/health";proxy=nothing).status==200
            warm=base*"/runtime/api/assignments/$(uuid4())/terminal"
            @test HTTP.get(warm;headers=headers("alice"),proxy=nothing,status_exception=false).status==404
            for (fn,types) in ((RT.relay_terminal_socket!,(typeof(service.terminals),RT.TerminalAttachment)),
                    (RT.terminal_browser_request,(typeof(service.terminals),RT.TerminalAttachment,Dict{String,Any})))
                Base.precompile(Tuple{typeof(fn),types...})
            end
            enroll_worker!(store,operator,trust)
            set_registration_state!(store,operator,"worker-a",:approved;expected_revision=1)
            start_agent!(agent);start_control!(service)
            @test timedwait(()->haskey(service.inventory.presence,"worker-a") && service.terminals.state==:online,15)==:ok
            a_run=reserve_run!(store,alice,definition);b_run=reserve_run!(store,bob,definition)
            a=admit(alice,a_run);b=admit(bob,b_run)
            aid,bid=UUID(a.fence.lease_id),UUID(b.fence.lease_id)
            @test timedwait(()->assignment_usable(service.coordinator,alice,aid) && assignment_usable(service.coordinator,bob,bid),5)==:ok
            url=replace(base,"http:"=>"ws:")*"/runtime/api/assignments/$aid/terminal"
            for hs in (Pair{String,String}[],headers("bob"),headers("operator"),headers("alice";origin="https://evil.test"))
                @test_throws Exception HTTP.WebSockets.open(url;headers=hs,proxy=nothing,cookies=false,request_timeout=5)
            end
            @test isempty(service.terminals.attachments) && isempty(resources.handles)
            non_upgrade=HTTP.get(base*"/runtime/api/assignments/$aid/terminal";headers=headers("alice"),proxy=nothing,status_exception=false)
            if non_upgrade.status!=405
                @info "terminal fixture authority" response=String(non_upgrade.body) state=get_assignment(store,alice,aid).state control=service.state agent=agent.state
            end
            @test non_upgrade.status==405
            wa,wb=string(uuid4()),string(uuid4())
            sa=attach(aid,"alice",wa);sb=attach(bid,"bob",wb)
            @test_throws Exception HTTP.WebSockets.open(url;headers=headers("alice"),proxy=nothing,cookies=false,request_timeout=5)
            first=action(sa,"open";columns=100,rows=30);second=action(sb,"open";columns=80,rows=24)
            @test first.accepted && second.accepted && first.session_id!=second.session_id
            ready(sa,first.session_id);ready(sb,second.session_id)
            heartbeat=@async begin
                next=time()+5
                while !heartbeat_stop[]
                    if time()>=next
                        @test action(sb,"keepalive";session=second.session_id).accepted
                        next=time()+5
                    end
                    sleep(0.05)
                end
            end
            input_id=string(uuid4());input=codeunits("private_counter=41\r")
            @test action(sa,"input";session=first.session_id,sequence=1,bytes=input,request_id=input_id).input_sequence==1
            @test action(sa,"input";session=first.session_id,sequence=1,bytes=input,request_id=input_id,retry=true).input_sequence==1
            @test action(sb,"input";session=second.session_id,sequence=1,
                bytes=codeunits("println(\"OTHER:\",isdefined(Main,:private_counter))\r")).accepted
            _,output=expect_output(sb,second.session_id,"OTHER:false")
            @test !occursin("private_counter=41",output)
            @test action(sa,"keepalive";session=first.session_id).accepted
            @test action(sa,"resize";session=first.session_id,columns=120,rows=40).accepted
            size=zeros(UInt16,4)
            @test ccall(:ioctl,Cint,(Cint,Culong,Ptr{UInt16}),resources.handles[a.fence.lease_id].process.master,0x5413,size)==0
            @test size[1:2]==[40,120]
            close(sa)
            @test timedwait(()->!haskey(service.terminals.attachments,string(aid)),7)==:ok
            @test !terminal_status(resources,a.fence,first.session_id).writer_connected
            resumed=attach(aid,"alice",wa)
            reopened=action(resumed,"open";columns=120,rows=40)
            @test reopened.accepted && reopened.session_id==first.session_id
            @test action(resumed,"input";session=first.session_id,sequence=2,
                bytes=codeunits("println(\"RETAINED:\",private_counter)\r")).accepted
            expect_output(resumed,first.session_id,"RETAINED:41")
            fresh=action(resumed,"restart";session=first.session_id,columns=100,rows=30)
            @test fresh.accepted && fresh.session_id!=first.session_id
            ready(resumed,fresh.session_id)
            @test action(resumed,"input";session=fresh.session_id,sequence=1,
                bytes=codeunits("println(\"NEW:\",isdefined(Main,:private_counter))\r")).accepted
            expect_output(resumed,fresh.session_id,"NEW:false")
            @test HTTP.get(base*"/health";proxy=nothing).status==200
            @test assignment_usable(service.coordinator,bob,bid)
            @test isempty(list_jobs(store,alice,a_run.id)) && isempty(list_jobs(store,bob,b_run.id))
            events=HTTP.get(base*"/runtime/api/control/events";headers=headers("alice"),proxy=nothing)
            @test !occursin("private_counter",String(events.body))
            transition_run!(store,alice,a_run.id,:stopped)
            @test timedwait(()->!isopen(resumed.readchannel),5)==:ok
            @test_throws HTTP.WebSockets.WebSocketError HTTP.WebSockets.receive(resumed)
            @test timedwait(()->!haskey(service.terminals.attachments,string(aid)),5)==:ok
            @test action(sb,"keepalive";session=second.session_id).accepted
            @test isempty(supervisor.handles) # gateway did not spawn a second UI/evaluator
            heartbeat_stop[]=true;wait(heartbeat)
            sb.close_transport!() # abrupt transport loss, no WebSocket close frame
            @test timedwait(()->!haskey(service.terminals.attachments,string(bid)),7)==:ok
            @test !terminal_status(resources,b.fence,second.session_id).writer_connected
            @test HTTP.get(base*"/health";proxy=nothing).status==200
            transition_run!(store,bob,b_run.id,:stopped)
            @test timedwait(()->isempty(resources.handles),10)==:ok
            browser_run=reserve_run!(store,alice,definition)
            try
                terminal_live_browser(service,supervisor,alice,browser_run,dir) do
                    admit(alice,browser_run)
                end
                @test isempty(list_jobs(store,alice,browser_run.id))
            finally
                transition_run!(store,alice,browser_run.id,:stopped)
            end
            @test timedwait(()->isempty(resources.handles) && isempty(service.terminals.attachments),10)==:ok
        finally
            heartbeat_stop[]=true
            try
                heartbeat === nothing || wait(heartbeat)
            finally
                foreach(socket->try close(socket) catch end,sockets)
                try close(service) finally
                    try close(agent) finally
                        try close(server) finally
                            try close(supervisor) finally;close(store);end
                        end
                    end
                end
            end
        end
        @test driver.closed && isempty(driver.handles) && isempty(resources.handles)
        @test isempty(service.terminals.attachments) && isempty(service.terminals.lanes)
    end
end
