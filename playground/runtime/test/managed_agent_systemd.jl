using Sockets
include("managed_agent.jl")

# Opt-in real service-manager gate. It starts the control agent only: the sole
# configured image is deliberately unavailable and is never pulled or executed.
function test_managed_agent_systemd(;native=false)
    directory = mktempdir(;prefix="lcm-managed-agent-",cleanup=false)
    worker = "audit-" * string(uuid4())
    listener = listen(ip"127.0.0.1",0)
    acceptor = @async try
        while isopen(listener)
            close(accept(listener))
        end
    catch error
        isopen(listener) && rethrow()
    end
    path,config,data = managed_config(directory;worker,port=Int(getsockname(listener)[2]))
    if native
        # This opt-in negative admission case requires a host that really lacks
        # native prerequisites. It must remain a control-only agent, not bypass
        # any limit in order to exercise the scientific launch branch.
        probe_runner = CommandRunner()
        try
            @test !isempty(check_native_host(probe_runner).failures)
            isempty(check_native_host(probe_runner).failures) && error("Native negative gate requires unsupported host limits.")
        finally
            close(probe_runner)
        end
        only(data["profiles"])["isolation"] = "trusted_process"
        only(data["profiles"])["environment"] = "/fixture/uninstalled-native-project"
        open(path,"w") do io; TOML.print(io,data); end
        config = read_agent_config(path)
    end
    runner = CommandRunner(;timeout_seconds=20)
    unit = ManagedRT.agent_unit_name(worker)
    identity = nothing
    object_path = nothing
    passed = false
    stopped = false
    seed_host = nothing
    baseline = nothing
    command(arguments) = run_owned_command!(runner,
        setenv(Cmd(arguments),ManagedRT.container_command_environment()))
    manager(arguments) = command([Sys.which("systemctl"),"--user","--no-pager",arguments...])
    function inventory(host)
        result = ManagedRT.scoped_container_command(runner,host,["container","ls","--all","--no-trunc","--format","{{.ID}}"])
        result.exitcode == 0 || error("Container inventory unavailable.")
        return sort(split(strip(result.output),'\n';keepempty=false))
    end
    function capture()
        records = ManagedRT.systemd_bus_command(runner,["call","org.freedesktop.systemd1",
            "/org/freedesktop/systemd1","org.freedesktop.systemd1.Manager","GetUnit","s",unit])
        only(only(records).data)
    end
    function current()
        properties = ManagedRT.systemd_properties(runner,object_path,"Unit",["Id","InvocationID","ActiveState"])
        service = ManagedRT.systemd_properties(runner,object_path,"Service",["MainPID","ExecStart","ExecStopPost"])
        return properties,service
    end
    function same_owned_unit()
        properties,service = current()
        ManagedRT.systemd_property(properties,"Id","s") == unit || error("Test unit identity changed.")
        ManagedRT.systemd_exec_matches(ManagedRT.systemd_property(service,"ExecStart","a(sasbttttuii)"),
            ManagedRT.agent_start_command(path)) || error("Test unit start command changed.")
        ManagedRT.systemd_exec_matches(ManagedRT.systemd_property(service,"ExecStopPost","a(sasbttttuii)"),
            ManagedRT.agent_recovery_command(config)) || error("Test unit recovery command changed.")
        identity === nothing || bytes2hex(UInt8.(ManagedRT.systemd_property(properties,"InvocationID","ay"))) ==
            identity.invocation || error("Test unit incarnation changed.")
        return properties,service
    end
    function wait_until(predicate;seconds=120)
        deadline = time() + seconds
        while time() < deadline
            predicate() && return true
            sleep(0.2)
        end
        return false
    end
    function logs()
        result = command([Sys.which("journalctl"),"--user","--unit",unit,"--output=cat","--no-pager"])
        result.exitcode == 0 || error("Owned test unit diagnostics unavailable.")
        return result.output
    end
    try
        before = manager(["show",unit,"--property=LoadState","--value"])
        strip(before.output) == "not-found" || error("Unique test unit already exists.")
        # A foreground caller cannot open ownership merely by naming a config.
        direct = command(ManagedRT.agent_start_command(path))
        @test direct.exitcode != 0
        @test !ispath(config.scratch_root)
        image = get(ENV,"LCM_TEST_STOPPED_CONTAINER_IMAGE","")
        if !isempty(image)
            occursin(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*@sha256:[a-f0-9]{64}$",image) ||
                error("Stopped recovery seed requires an already cached digest-pinned image.")
            seed_host = check_container_host(runner;requested="podman")
            inspected = ManagedRT.scoped_container_command(runner,seed_host,["image","inspect","--format","{{json .}}",image])
            inspected.exitcode == 0 && image in ManagedRT.JSON3.read(inspected.output).RepoDigests ||
                error("Stopped recovery seed image is not cached; no pull is permitted.")
            baseline = inventory(seed_host)
            journal = ResourceJournal(ManagedRT.agent_journal_root(config),worker)
            try
                profile = ProfileDefinition("old-fixture",image,last(split(image,"@sha256:"));
                    kind=:scientific,isolation=:container,operations=("fixture.echo",))
                fence = ManagedRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"audit","old-run",
                    worker,string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
                receipt = reserve_resource!(journal,fence,seed_host.engine.name,container_scope(runner,seed_host))
                # Deliberately create-only. This tests stale ownership recovery,
                # not image admission, executor launch, or effective limits.
                created = ManagedRT.scoped_container_command(runner,seed_host,container_create_arguments(ContainerPolicy(profile,receipt)))
                created.exitcode == 0 || error("Stopped recovery seed acquisition failed.")
                inspected = ManagedRT.scoped_container_command(runner,seed_host,["container","inspect","--format","{{json .}}",strip(created.output)])
                inspected.exitcode == 0 || error("Stopped recovery seed inspection failed.")
                object = ManagedRT.JSON3.read(inspected.output)
                @test object.State.Running === false && object.State.Pid == 0
                @test length(inventory(seed_host)) == length(baseline) + 1
            finally
                close(journal)
            end
        end
        arguments = [Sys.which("systemd-run"),"--user","--quiet","--unit=" * unit,
            # Tests may use a private provisioned depot. Scope it to this unit;
            # never alter the user's manager environment or unrelated services.
            "--setenv=JULIA_DEPOT_PATH=" * join(DEPOT_PATH,':'),
            "--setenv=JULIA_LOAD_PATH=@:@stdlib","--setenv=OPENBLAS_NUM_THREADS=1",
            "--property=Type=exec","--property=ExitType=main","--property=RemainAfterExit=no","--property=KillMode=control-group",
            "--property=KillSignal=SIGINT","--property=SendSIGKILL=yes","--property=FinalKillSignal=SIGKILL",
            "--property=TimeoutStartSec=120","--property=TimeoutStopSec=60","--property=TimeoutStopFailureMode=terminate","--property=Restart=no",
            "--property=UMask=0077","--property=StandardInput=null","--property=StandardOutput=journal",
            "--property=StandardError=journal","--property=ExecStopPost=" * ManagedRT.systemd_command(ManagedRT.agent_recovery_command(config)),
            "--",ManagedRT.agent_start_command(path)...]
        launched = command(arguments)
        launched.exitcode == 0 || error("Isolated service launch failed.")
        object_path = capture()
        properties,service = same_owned_unit()
        pid = Int(ManagedRT.systemd_property(service,"MainPID","u"))
        identity = verify_managed_agent(runner,path,config;pid,read_cgroup=()->read("/proc/$pid/cgroup",String))
        @test wait_until(()->begin
            properties,_ = same_owned_unit()
            ManagedRT.systemd_property(properties,"ActiveState","s") == "active" || error("Agent exited during startup.")
            occursin("Control scheduling active; no executor prepared on startup.",logs())
        end)
        @test occursin("0 eligible profiles",logs())
        seed_host === nothing || @test inventory(seed_host) == baseline
        @test isfile(joinpath(identity.journal_root,"owner.json"))
        @test_throws ArgumentError ResourceJournal(identity.journal_root,worker)
        @test_throws ArgumentError recover_agent_resources!(identity.journal_root,worker)
        @test verify_managed_agent(runner,path,config;previous=identity,pid,
            read_cgroup=()->read("/proc/$pid/cgroup",String)).invocation == identity.invocation
        same_owned_unit()
        @test manager(["kill","--kill-whom=main","--signal=KILL",unit]).exitcode == 0
        @test wait_until(()->begin
            properties,service = same_owned_unit()
            ManagedRT.systemd_property(properties,"ActiveState","s") == "failed" &&
                ManagedRT.systemd_property(service,"MainPID","u") == 0
        end)
        _,service = same_owned_unit()
        recovery = only(ManagedRT.systemd_property(service,"ExecStopPost","a(sasbttttuii)"))
        @test recovery[4] > 0 && recovery[6] > 0
        @test recovery[9] == 1 && recovery[10] == 0 # CLD_EXITED, success
        @test occursin("Owned agent recovery complete; no lease or preparation restored.",logs())
        journal = ResourceJournal(identity.journal_root,worker)
        try
            @test isempty(resource_receipts(journal))
        finally
            close(journal)
        end
        stopped = true
        @test recover_agent_resources!(identity.journal_root,worker) === nothing
        seed_host === nothing || @test inventory(seed_host) == baseline
        passed = true
    finally
        try
            if object_path !== nothing
                properties,_ = same_owned_unit()
                stopped = ManagedRT.systemd_property(properties,"ActiveState","s") in ("failed","inactive")
                if ManagedRT.systemd_property(properties,"ActiveState","s") in ("active","activating","deactivating")
                    manager(["stop","--no-block",unit]).exitcode == 0 || error("Owned test service stop failed.")
                    stopped = wait_until(()->begin
                        result = manager(["show",unit,"--property=ActiveState","--value"])
                        strip(result.output) in ("failed","inactive")
                    end)
                end
                text = logs()
                log_path = joinpath(directory,"agent.log")
                open(log_path,"w") do io; chmod(log_path,0o600); write(io,text); end
                if stopped
                    properties,_ = same_owned_unit()
                    ManagedRT.systemd_property(properties,"ActiveState","s") == "failed" &&
                        manager(["reset-failed",unit]).exitcode != 0 && error("Owned failed unit reset failed.")
                end
            end
            if seed_host !== nothing && (object_path === nothing || stopped)
                recover_agent_resources!(ManagedRT.agent_journal_root(config),worker)
                inventory(seed_host) == baseline || error("Container inventory differs after owned recovery.")
            end
        finally
            close(listener);wait(acceptor);close(runner)
            println("Managed-agent diagnostics: ",directory)
        end
    end
    passed && stopped || error("Managed-agent lifecycle remains unresolved.")
    println("Real managed-agent crash/recovery probe finished; no numerical or terminal executor started.")
end

@testset "actual managed control-agent lifecycle" begin
    test_managed_agent_systemd()
    get(ENV,"LCM_TEST_NATIVE_UNAVAILABLE","") == "1" && test_managed_agent_systemd(;native=true)
end
