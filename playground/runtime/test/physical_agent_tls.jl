# Opt-in real managed agent, TLS broker, gateway, leases and container REPL.
# No test driver, in-process agent or fabricated isolation evidence is used.
using Test, UUIDs, TOML, LineCableModelsRuntime
import HTTP, JSON3
const RT = LineCableModelsRuntime
include("physical_proxy.jl")
include("physical_science.jl")

function physical_agent_gate(directory, port)
    engine = get(ENV,"LCM_TEST_CONTAINER_RUNTIME","")
    engine in ("docker","podman") || error("Select the physical agent engine explicitly.")
    image = ENV["LCM_TEST_TERMINAL_IMAGE"]
    certs = joinpath(directory,"certs")
    endpoint(id) = BrokerEndpoint("tls://127.0.0.1:$port",joinpath(directory,id * ".password");
        ca_file=joinpath(certs,"ca.pem"),certificate_file=joinpath(certs,"worker-cert.pem"),
        key_file=joinpath(certs,"worker-key.pem"),server_name="localhost")
    root = joinpath(directory,"physical-agent")
    mkpath(root)
    scientific_profiles = physical_scientific_profiles()
    capacity = 2 + length(scientific_profiles)
    path = joinpath(root,"agent.toml")
    data = Dict("schema_version"=>1,
        "agent"=>Dict("worker_id"=>"worker-a","scratch_root"=>joinpath(root,"owned"),
            "capacity"=>capacity,"container_runtime"=>engine),
        "broker"=>Dict("url"=>"tls://127.0.0.1:$port","password_file"=>joinpath(directory,"worker-a.password"),
            "ca_file"=>joinpath(certs,"ca.pem"),"certificate_file"=>joinpath(certs,"worker-cert.pem"),
            "key_file"=>joinpath(certs,"worker-key.pem"),"server_name"=>"localhost"),
        "profiles"=>[Dict("id"=>"julia-terminal","kind"=>"terminal","isolation"=>"container",
            "environment"=>image,"fingerprint"=>last(split(image,"@sha256:")),
            "budget"=>Dict("cpus"=>0.5,"memory_bytes"=>512*1024^2,"pids"=>64,"scratch_bytes"=>8*1024^2))])
    append!(data["profiles"],scientific_profiles)
    open(io->TOML.print(io,data),path,"w"); chmod(path,0o600)
    config = read_agent_config(path)
    runner = CommandRunner(;timeout_seconds=20)
    command(args) = run_owned_command!(runner,setenv(Cmd(args),RT.container_command_environment()))
    manager(args) = command([Sys.which("systemctl"),"--user","--no-pager",args...])
    unit = RT.agent_unit_name(config.worker_id)
    strip(manager(["show",unit,"--property=LoadState","--value"]).output) == "not-found" ||
        error("Refusing an existing worker-a service; use an isolated acceptance account.")
    host = check_container_host(runner;requested=engine)
    function inventory()
        result = RT.scoped_container_command(runner,host,["container","ls","--all","--no-trunc","--format","{{.ID}}"])
        result.exitcode == 0 || error("Container inventory unavailable")
        sort(split(strip(result.output),'\n';keepempty=false))
    end
    baseline = inventory()
    identity = nothing
    launched = false
    retired = false
    function checked_unit()
        record = only(RT.systemd_bus_command(runner,["call","org.freedesktop.systemd1",
            "/org/freedesktop/systemd1","org.freedesktop.systemd1.Manager","GetUnit","s",unit]))
        object = only(record.data)
        properties = RT.systemd_properties(runner,object,"Unit",["Id","InvocationID","ActiveState"])
        service = RT.systemd_properties(runner,object,"Service",["MainPID","ExecStart","ExecStopPost"])
        RT.systemd_property(properties,"Id","s") == unit &&
            RT.systemd_exec_matches(RT.systemd_property(service,"ExecStart","a(sasbttttuii)"),RT.agent_start_command(path)) &&
            RT.systemd_exec_matches(RT.systemd_property(service,"ExecStopPost","a(sasbttttuii)"),RT.agent_recovery_command(config)) ||
            error("Owned service identity changed")
        identity === nothing || bytes2hex(UInt8.(RT.systemd_property(properties,"InvocationID","ay"))) == identity.invocation ||
            error("Owned service incarnation changed")
        return properties,service
    end
    function start_service()
        args = [Sys.which("systemd-run"),"--user","--quiet","--unit=" * unit,
            "--setenv=JULIA_DEPOT_PATH=" * join(DEPOT_PATH,':'),"--setenv=JULIA_LOAD_PATH=@:@stdlib",
            "--setenv=OPENBLAS_NUM_THREADS=1","--property=Type=exec","--property=ExitType=main",
            "--property=RemainAfterExit=no","--property=KillMode=control-group","--property=KillSignal=SIGINT",
            "--property=SendSIGKILL=yes","--property=FinalKillSignal=SIGKILL","--property=TimeoutStartSec=120",
            "--property=TimeoutStopSec=60","--property=TimeoutStopFailureMode=terminate","--property=Restart=no",
            "--property=UMask=0077","--property=StandardInput=null","--property=StandardOutput=journal",
            "--property=StandardError=journal","--property=ExecStopPost=" * RT.systemd_command(RT.agent_recovery_command(config)),
            "--",RT.agent_start_command(path)...]
        command(args).exitcode == 0 || error("Managed agent launch failed")
        launched = true
        _,service = checked_unit()
        pid = Int(RT.systemd_property(service,"MainPID","u"))
        identity = verify_managed_agent(runner,path,config;pid,read_cgroup=()->read("/proc/$pid/cgroup",String))
    end
    applications = ApplicationRegistry()
    definition = ApplicationDefinition("physical-terminal","Physical terminal",:workbench,"/physical-terminal";
        requirements=(RuntimeRequirement("terminal",("julia-terminal",)),
            (RuntimeRequirement(p["id"],(p["id"],)) for p in scientific_profiles)...))
    register!(applications,definition)
    store = RuntimeStore(joinpath(root,"runtime.sqlite"))
    trust = WorkerTrust("worker-a","worker-a",Tuple(sort!(collect(keys(config.profiles.definitions))));capacity)
    control = ControlService(ControlConfig(endpoint("coordinator"),config.profiles,[trust]),store,applications)
    supervisor = UIHostSupervisor(store,applications,joinpath(root,"hosts"))
    proxy_key = string(uuid4()) * string(uuid4())
    with_proxy = haskey(ENV,"LCM_TEST_CADDY")
    proxy_port = with_proxy ? physical_proxy_port() : 0
    origin = with_proxy ? "https://127.0.0.1:$proxy_port" : "https://lcm.test"
    server = start_gateway(supervisor,ProxyIdentity(origin,["127.0.0.1"],proxy_key;
        administrators=["operator"]);control)
    direct = "http://127.0.0.1:$(HTTP.port(server))"
    proxy = nothing
    client_options = (;proxy=nothing)
    base = direct
    alice_id,bob_id = with_proxy ? ("operator","researcher") : ("alice","bob")
    headers(user) = with_proxy ? proxy.headers(user=="alice" ? alice_id : user=="bob" ? bob_id : user) :
        ["X-LCM-Proxy-Key"=>proxy_key,"X-LCM-Principal"=>user,"Origin"=>origin]
    alice,bob,operator = Principal(alice_id),Principal(bob_id),Principal("operator";administrator=true)
    terminal_client = nothing
    client_config = joinpath(root,"client.json")
    function client_signal(name,seconds)
        @test timedwait(()->isfile(joinpath(root,name)) || process_exited(terminal_client),seconds;pollint=0.1)==:ok
        isfile(joinpath(root,name)) || error("Private terminal client exited before $name")
    end
    function online()
        lock(control.inventory.lock) do
            RT.presence_state(control.inventory,get(control.inventory.presence,"worker-a",nothing)) == :online
        end
    end
    function admit(owner,run;role="terminal",profile="julia-terminal")
        @test timedwait(online,120;pollint=0.1) == :ok
        lease = reserve_assignment!(control.coordinator.assignments,owner,run.id,role,profile;
            placement=PinnedPlacement("worker-a"))
        grant_assignment!(control.coordinator,owner,UUID(lease.fence.lease_id))
        @test timedwait(()->assignment_usable(control.coordinator,owner,UUID(lease.fence.lease_id)),5) == :ok
        return lease
    end
    try
        if with_proxy
            proxy = start_physical_proxy(directory,HTTP.port(server),proxy_port,proxy_key)
            base = proxy.origin; client_options = (;client=proxy.client)
            @test HTTP.get(base * "/runtime/api/runs";client_options...,status_exception=false).status == 401
            @test HTTP.get(direct * "/runtime/api/runs";proxy=nothing,status_exception=false).status == 401
            @test HTTP.get(base * "/runtime/api/runs";client_options...,headers=["X-LCM-Principal"=>"operator",
                "X-LCM-Proxy-Key"=>proxy_key],status_exception=false).status == 401
            spoofed = [headers("bob"); "X-LCM-Principal"=>"operator"; "X-LCM-Proxy-Key"=>"forged"]
            response = HTTP.get(base * "/runtime/api/control";client_options...,headers=spoofed)
            @test !JSON3.read(response.body).administrator
            @test HTTP.post(base * "/runtime/api/runs";client_options...,
                headers=[first(headers("alice")),"Origin"=>"https://evil.test"],
                body="{}",status_exception=false).status == 403
        end
        @test HTTP.get(base * "/health";client_options...).status == 200
        @test HTTP.get(base * "/runtime/api/assignments/$(uuid4())/terminal";
            client_options...,headers=headers("alice"),status_exception=false).status == 404
        enroll_worker!(store,operator,trust)
        set_registration_state!(store,operator,"worker-a",:approved;expected_revision=1)
        start_control!(control); start_service()
        @test timedwait(()->online() && control.terminals.state == :online,120;pollint=0.1) == :ok
        @test inventory() == baseline # startup does not prepare user executors
        run_a,run_b = reserve_run!(store,alice,definition),reserve_run!(store,bob,definition)
        a,b = admit(alice,run_a),admit(bob,run_b)
        payload = (;base,lease_a=a.fence.lease_id,lease_b=b.fence.lease_id,run_a=string(run_a.id),scientific_profiles,
            ca_file=with_proxy ? joinpath(certs,"ca.pem") : nothing,
            headers=Dict(user=>[[p.first,p.second] for p in headers(user)] for user in ("alice","bob")))
        open(client_config,"w") do io; chmod(client_config,0o600);write(io,JSON3.write(payload));end
        executable=joinpath(Sys.BINDIR,Base.julia_exename())
        project=normpath(joinpath(@__DIR__,"..")); script=joinpath(@__DIR__,"physical_terminal_client.jl")
        terminal_client=run(pipeline(`$executable --startup-file=no --threads=2 --project=$project $script $root`;
            stdin=devnull,stdout=stdout,stderr=stderr);wait=false)
        client_signal("terminal-ready",180)
        @test length(inventory()) == length(baseline) + 2
        @test HTTP.get(base * "/health";client_options...).status == 200
        @test isempty(list_jobs(store,alice,run_a.id)) && isempty(list_jobs(store,bob,run_b.id))
        client_signal("science-ready",1500)
        write(joinpath(root,"retire-agent"),"retire\n")
        client_signal("client-ready-to-retire",20)
        checked_unit()
        if get(ENV,"LCM_TEST_AGENT_STOP","crash") == "graceful"
            @test manager(["stop","--no-block",unit]).exitcode == 0
            @test timedwait(45;pollint=0.2) do
                state = strip(manager(["show",unit,"--property=ActiveState","--value"]).output)
                state in ("inactive","failed")
            end == :ok
            state = strip(manager(["show",unit,"--property=ActiveState","--value"]).output)
            @test state == "inactive"
            events = command([Sys.which("journalctl"),"--user","_SYSTEMD_INVOCATION_ID=" * identity.invocation,
                "--output=cat","--no-pager"]).output
            @test occursin("LCM agent shutdown complete; owned resources retired.",events)
            @test occursin("Owned agent recovery complete; no lease or preparation restored.",events)
            retired = state == "inactive"
        else
            @test manager(["kill","--kill-whom=main","--signal=KILL",unit]).exitcode == 0
            @test timedwait(90;pollint=0.2) do
                properties,service = checked_unit()
                RT.systemd_property(properties,"ActiveState","s") == "failed" && RT.systemd_property(service,"MainPID","u") == 0
            end == :ok
            _,service = checked_unit()
            recovery = only(RT.systemd_property(service,"ExecStopPost","a(sasbttttuii)"))
            @test recovery[9] == 1 && recovery[10] == 0
        end
        @test inventory() == baseline
        @test timedwait(()->!assignment_usable(control.coordinator,alice,UUID(a.fence.lease_id)),20) == :ok
        @test HTTP.get(base * "/health";client_options...).status == 200
        journal = ResourceJournal(RT.agent_journal_root(config),"worker-a")
        try @test isempty(resource_receipts(journal)) finally close(journal) end
        @test timedwait(()->process_exited(terminal_client),65;pollint=0.1)==:ok
        wait(terminal_client); @test success(terminal_client)
    finally
        if terminal_client!==nothing && !process_exited(terminal_client)
            kill(terminal_client,Base.SIGTERM)
            timedwait(()->process_exited(terminal_client),5)==:ok || kill(terminal_client,Base.SIGKILL)
            wait(terminal_client)
        end
        isfile(client_config) && rm(client_config)
        write(joinpath(root,"control-snapshot.json"),JSON3.write(RT.control_snapshot(control,operator)))
        if launched && !retired
            properties,_ = checked_unit()
            if !(RT.systemd_property(properties,"ActiveState","s") in ("inactive","failed"))
                manager(["stop","--no-block",unit]).exitcode == 0 || error("Owned agent stop failed")
                timedwait(()->strip(manager(["show",unit,"--property=ActiveState","--value"]).output) in ("inactive","failed"),90) == :ok ||
                    error("Owned agent stop unresolved")
            end
            strip(manager(["show",unit,"--property=LoadState","--value"]).output)=="not-found" || begin
                properties,_ = checked_unit()
                RT.systemd_property(properties,"ActiveState","s") == "failed" && manager(["reset-failed",unit])
            end
        end
        if launched
            logs = command([Sys.which("journalctl"),"--user","_SYSTEMD_INVOCATION_ID=" * identity.invocation,
                "--output=cat","--no-pager"])
            write(joinpath(root,"agent.log"),logs.output)
        end
        proxy === nothing || proxy.stop()
        close(server); close(supervisor); close(control); close(store)
        try
            recover_agent_resources!(RT.agent_journal_root(config),"worker-a")
            inventory() == baseline || error("Owned container cleanup unresolved")
        finally
            close(runner)
        end
        println("Physical agent diagnostics: ",root)
    end
end

@testset "physical managed agent and private terminal over TLS" begin
    physical_agent_gate(ARGS...)
end
