# A separate browser-protocol client process: its TLS/JIT work must not pause
# the gateway's lease scheduler. No coordinator, agent or engine runs here.
using Test, UUIDs, HTTP, JSON3
include("physical_science.jl")

function physical_terminal_client(directory)
    config = JSON3.read(read(joinpath(directory,"client.json"),String),Dict{String,Any})
    base = config["base"]
    client = config["ca_file"] === nothing ? nothing : HTTP.Client(transport=HTTP.Transport(
        tls_config=HTTP.TLS.Config(ca_file=config["ca_file"]),proxy=nothing))
    options = client === nothing ? (;proxy=nothing) : (;client)
    headers(user) = [String(p[1])=>String(p[2]) for p in config["headers"][user]]
    sockets = Any[]; locks = IdDict{Any,ReentrantLock}(); sessions = IdDict{Any,String}()
    heartbeat = nothing; stopping = Ref(false)
    function action(socket,name;session=nothing,sequence=0,after=0,columns=0,rows=0,bytes=UInt8[])
        lock(locks[socket]) do
            request_id = string(uuid4())
            HTTP.WebSockets.send(socket,JSON3.write((action=name,request_id,session_id=session,
                input_sequence=sequence,after,columns,rows,bytes=collect(bytes),retry=false)))
            reply = JSON3.read(HTTP.WebSockets.receive(socket))
            reply.kind=="report" && reply.request_id==request_id || error("Unexpected terminal reply")
            @test reply.accepted
            reply.accepted || error("Terminal action rejected: $(reply.reason)")
            name in ("open","restart") && (sessions[socket]=reply.session_id)
            reply
        end
    end
    function attach(lease,user)
        socket = HTTP.WebSockets.open(replace(base,r"^http"=>"ws") * "/runtime/api/assignments/$lease/terminal";
            options...,headers=headers(user),cookies=false,request_timeout=5,maxframesize=65536,read_idle_timeout=10)
        push!(sockets,socket); locks[socket]=ReentrantLock()
        @test JSON3.read(HTTP.WebSockets.receive(socket)).kind=="hello"
        HTTP.WebSockets.send(socket,JSON3.write((action="attach",writer_id=string(uuid4()))))
        @test JSON3.read(HTTP.WebSockets.receive(socket)).kind=="attached"
        socket
    end
    function ready(socket,id)
        state = Ref{Any}(nothing)
        @test timedwait(60;pollint=0.1) do
            state[]=action(socket,"status";session=id)
            state[].phase!="starting"
        end == :ok
        state[]!==nothing && state[].phase=="ready" || error("Physical REPL did not become ready")
    end
    function expect(socket,id,needle)
        output=UInt8[]; cursor=0
        @test timedwait(15;pollint=0.1) do
            report=action(socket,"read";session=id,after=cursor)
            cursor=report.cursor;append!(output,UInt8.(report.bytes))
            occursin(needle,String(copy(output)))
        end == :ok
    end
    try
        @test HTTP.get(base*"/health";options...).status==200
        url=replace(base,r"^http"=>"ws") * "/runtime/api/assignments/$(config["lease_a"])/terminal"
        @test_throws Exception HTTP.WebSockets.open(url;options...,headers=headers("bob"),request_timeout=5)
        sa,sb=attach(config["lease_a"],"alice"),attach(config["lease_b"],"bob")
        opened_a=action(sa,"open";columns=100,rows=30)
        opened_b=action(sb,"open";columns=80,rows=24)
        @test opened_a.session_id!=opened_b.session_id
        heartbeat=@async begin
            next=time()+2
            while !stopping[]
                if time()>=next
                    for socket in sockets
                        lock(locks[socket]) do; action(socket,"keepalive";session=sessions[socket]); end
                    end
                    next=time()+2
                end
                sleep(0.05)
            end
        end
        ready(sa,opened_a.session_id);ready(sb,opened_b.session_id)
        action(sa,"input";session=opened_a.session_id,sequence=1,
            bytes=codeunits("private_counter=41; println(\"PRIVATE_\",\"VALUE:\",private_counter+1)\r"))
        expect(sa,opened_a.session_id,"PRIVATE_VALUE:42")
        action(sb,"input";session=opened_b.session_id,sequence=1,
            bytes=codeunits("println(\"SEPARATE_\",\"STATE:\",isdefined(Main,:private_counter))\r"))
        expect(sb,opened_b.session_id,"SEPARATE_STATE:false")
        action(sa,"resize";session=opened_a.session_id,columns=120,rows=40)
        fresh=action(sa,"restart";session=opened_a.session_id,columns=100,rows=30)
        @test fresh.session_id!=opened_a.session_id
        ready(sa,fresh.session_id)
        action(sa,"input";session=fresh.session_id,sequence=1,
            bytes=codeunits("println(\"FRESH_\",\"STATE:\",isdefined(Main,:private_counter))\r"))
        expect(sa,fresh.session_id,"FRESH_STATE:false")
        write(joinpath(directory,"terminal-ready"),"ready\n")
        function admit_science(id)
            response=HTTP.post(base*"/runtime/api/runs/$(config["run_a"])/assignments";
                options...,headers=[headers("alice");"X-LCM-Request"=>"1"],status_exception=false,
                body=JSON3.write((role=id,profile=id,placement=(mode="pinned",worker_id="worker-a"),request_id=string(uuid4()))))
            response.status==202 || error("Physical scientific assignment was not accepted: HTTP $(response.status)")
            @test response.status==202
            lease=JSON3.read(response.body).id
            usable=timedwait(5;pollint=0.1) do
                response=HTTP.get(base*"/runtime/api/assignments/$lease";options...,headers=headers("alice"))
                JSON3.read(response.body).usable
            end
            @test usable==:ok
            usable==:ok || error("Physical scientific assignment was not acknowledged")
            return lease
        end
        profiles=get(config,"scientific_profiles",Any[])
        if !isempty(profiles)
            results=physical_scientific_checks(base,options,headers,admit_science,profiles)
            write(joinpath(directory,"scientific-results.json"),JSON3.write(results))
        end
        write(joinpath(directory,"science-ready"),"ready\n")
        @test timedwait(()->isfile(joinpath(directory,"retire-agent")),1500;pollint=0.1)==:ok
        action(sb,"input";session=opened_b.session_id,sequence=2,
            bytes=codeunits("println(\"SURVIVED_\",\"SCIENCE\")\r"))
        expect(sb,opened_b.session_id,"SURVIVED_SCIENCE")
        stopping[]=true;wait(heartbeat)
        write(joinpath(directory,"client-ready-to-retire"),"ready\n")
        @test timedwait(()->all(s->!isopen(s.readchannel),sockets),60;pollint=0.1)==:ok
        @test HTTP.get(base*"/health";options...).status==200
    finally
        stopping[]=true
        heartbeat===nothing || wait(heartbeat;throw=false)
        for socket in sockets; try close(socket) catch end; end
        client===nothing || close(client)
    end
end

@testset "separate private terminal client" begin
    physical_terminal_client(only(ARGS))
end
