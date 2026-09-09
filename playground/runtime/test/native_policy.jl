using Test, UUIDs, LineCableModelsRuntime
const PolicyRT = LineCableModelsRuntime

function native_policy_fixture(f;bound=false)
    receipt = f.add!(;bound)
    profile = ProfileDefinition("fixture","/fixture/project",repeat("a",64);operations=["fixture.echo"])
    policy = NativePolicy(profile,receipt)
    manager = ManagedAgentIdentity("lcm-agent-worker-a.service",repeat("b",32),"/private/agent",
        "/private/config",f.journal.root,"worker-a")
    host = NativeHostCheck(f.scope,"/usr/bin/systemd-run",1000,1000,())
    environment = EnvironmentFingerprint(profile.fingerprint,"Fixture",uuid4(),1,1)
    options = (;project="/fixture/project",julia="/fixture/julia",depots=["/fixture/depot"],
        env_executable="/usr/bin/env",guard="/fixture/native-guard.jl",entry="/fixture/native-scientific.jl")
    return (;receipt,profile,policy,manager,host,environment,options)
end

@testset "native launch policy is fixed, inert and credential filtered" begin
    with_native_fixture() do f
        p = native_policy_fixture(f)
        command = native_launch_command(p.policy,p.manager,p.host,p.environment;p.options...)
        args = collect(command)
        @test first(args) == "/usr/bin/systemd-run"
        @test all(in(args),("--user","--pipe","--wait","--quiet","--slice=app.slice"))
        @test !("--pty" in args) && !("-e" in args) && !("--shell" in args)
        @test "--property=CPUQuota=100.0%" in args && "--property=MemorySwapMax=0" in args
        @test "--property=TasksMax=128" in args && "--property=MemoryMax=1073741824" in args
        @test "--property=BindsTo=lcm-agent-worker-a.service" in args
        @test "--property=After=lcm-agent-worker-a.service" in args
        @test "--property=KillMode=control-group" in args && "--property=TimeoutStopSec=2" in args
        @test "--property=ProtectSystem=strict" in args && "--property=ReadOnlyPaths=/dev" in args
        @test "--property=TemporaryFileSystem=/tmp:rw,nosuid,nodev,noexec,size=268369920,mode=1777 /dev/shm:rw,nosuid,nodev,noexec,size=65536,mode=1777" in args
        index = findfirst(==("--"),args)
        @test args[index+1:index+2] == ["/usr/bin/env","-i"]
        @test "--load=/fixture/native-guard.jl" in args && last(args)=="Fixture"
        @test "JULIA_DEPOT_PATH=/tmp/depot:/fixture/depot" in args
        @test !any(a->occursin("NATS",a) || occursin("broker",a) || occursin("DBUS",a),args[index+1:end])
        @test "DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus" in command.env
        @test !occursin("/fixture",repr(MIME"text/plain"(),p.policy))
        @test only(resource_receipts(f.journal)).physical_id === nothing && isempty(f.mutations)
        for path in ("/tmp/source","/dev/shm/source","/bad/../source","relative","/bad\$name","/bad%name","/bad:name","/bad\nname")
            @test_throws ArgumentError native_launch_command(p.policy,p.manager,p.host,p.environment;p.options...,project=path)
        end
        @test_throws ArgumentError native_launch_command(p.policy,p.manager,p.host,p.environment;p.options...,depots=String[])
        @test_throws ArgumentError native_launch_command(p.policy,p.manager,p.host,p.environment;p.options...,depots=["/tmp/depot"])
        bad = NativeHostCheck(f.scope,p.host.command,1000,1000,(:cpu_controller_missing,))
        @test_throws CommandFailure native_launch_command(p.policy,p.manager,bad,p.environment;p.options...)
        foreign = ManagedAgentIdentity("lcm-agent-foreign.service",repeat("b",32),"/private/agent","/private/config",f.journal.root,"foreign")
        @test_throws ArgumentError native_launch_command(p.policy,foreign,p.host,p.environment;p.options...)
        changed = EnvironmentFingerprint(repeat("c",64),"Fixture",p.environment.uuid,1,1)
        @test_throws ArgumentError native_launch_command(p.policy,p.manager,p.host,changed;p.options...)
    end
end

@testset "native host prerequisites do not allocate or invent missing quotas" begin
    files = Dict("/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/cgroup.controllers"=>"cpu memory pids",
        "/proc/sys/user/max_user_namespaces"=>"100")
    check() = PolicyRT.native_host_failures(1000,p->files[p],_->(ftype=0x63677270,))
    @test isempty(check())
    for name in ("cpu","memory","pids")
        files["/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/cgroup.controllers"] = join(setdiff(["cpu","memory","pids"],[name])," ")
        @test check() == (Symbol(name*"_controller_missing"),)
    end
    files["/sys/fs/cgroup/user.slice/user-1000.slice/user@1000.service/cgroup.controllers"] = "cpu memory pids"
    files["/proc/sys/user/max_user_namespaces"] = "0"
    @test check() == (:user_namespace_unavailable,)
    @test :cgroup_v2_required in PolicyRT.native_host_failures(1000,p->files[p],_->(ftype=0x1234,))
end

@testset "native service verification binds identity and rejects configured drift" begin
    with_native_fixture() do f
        p = native_policy_fixture(f)
        unit = f.units[resource_name(p.receipt)]
        record(type,data) = Dict{String,Any}("type"=>type,"data"=>data)
        specifications = [("CPUQuotaPerSecUSec","t",1_000_000),("CPUQuotaPeriodUSec","t",100_000),
            ("MemoryMax","t",1024^3),("MemorySwapMax","t",0),("TasksMax","t",128),
            ("NoNewPrivileges","b",true),("PrivateUsers","b",true),("PrivateDevices","b",true),
            ("PrivateNetwork","b",true),("PrivateIPC","b",true),("ProtectSystem","s","strict"),
            ("ProtectControlGroups","b",true),("ProtectKernelTunables","b",true),
            ("WorkingDirectory","s","/tmp"),("CapabilityBoundingSet","t",0),("AmbientCapabilities","t",0)]
        for (key,type,value) in specifications; unit[key] = record(type,value); end
        for key in ("BindsTo","After"); unit[key] = record("as",[p.manager.unit]); end
        command = PolicyRT.native_exec_arguments(p.policy,p.host,p.environment;p.options...)
        unit["ExecStart"] = record("a(sasbttttuii)",[Any[first(command),command,false,0,0,0,0,0,0,0]])
        for key in ("ExecStartPre","ExecStartPost"); unit[key] = record("a(sasbttttuii)",[]); end
        check() = verify_native_service!(f.journal,f.runner,p.policy,p.manager,p.host,p.environment;
            command_options=p.options,bus_options=(;invoke=f.bus_options.invoke,which=f.bus_options.which))
        bound = check()
        @test bound.physical_id == bytes2hex(UInt8.(1:16))
        @test only(resource_receipts(f.journal)).physical_id == bound.physical_id
        @test check().physical_id == bound.physical_id && isempty(f.mutations)
        for (key,_,original) in specifications
            unit[key]["data"] = original isa Bool ? !original : original isa Integer ? original+1 : "other"
            @test_throws CommandFailure check()
            unit[key]["data"] = original
        end
        for key in ("BindsTo","After","ExecStart")
            previous = unit[key]["data"]; unit[key]["data"] = []
            @test_throws CommandFailure check()
            unit[key]["data"] = previous
        end
        for key in ("ExecStartPre","ExecStartPost")
            unit[key]["data"] = unit["ExecStart"]["data"]
            @test_throws CommandFailure check()
            unit[key]["data"] = []
        end
        unit["InvocationID"]["data"] = collect(2:17)
        @test_throws ArgumentError check()
        @test only(resource_receipts(f.journal)).physical_id == bound.physical_id
    end
end
