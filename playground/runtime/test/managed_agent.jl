using Test, TOML, UUIDs, LineCableModelsRuntime
const ManagedRT = LineCableModelsRuntime

function managed_config(directory;worker="worker-a",port=4222)
    secret = joinpath(directory,"broker.password")
    write(secret,"managed-fixture-private-secret");chmod(secret,0o600)
    path = joinpath(directory,"agent.toml")
    data = Dict("schema_version"=>1,
        "agent"=>Dict("worker_id"=>worker,"scratch_root"=>joinpath(directory,"owned-agent"),"container_runtime"=>"podman"),
        "broker"=>Dict("url"=>"nats://127.0.0.1:$port","password_file"=>secret,"allow_loopback_plaintext"=>true),
        "profiles"=>[Dict("id"=>"fixture","kind"=>"scientific","isolation"=>"container",
            "environment"=>"registry.invalid/lcm@sha256:" * repeat("a",64),"fingerprint"=>repeat("a",64),
            "operations"=>["fixture.echo"])])
    open(path,"w") do io
        chmod(path,0o600); TOML.print(io,data)
    end
    return path,read_agent_config(path),data
end

function managed_bus_fixture(path,config)
    record(type,data) = Dict{String,Any}("type"=>type,"data"=>data)
    exec_record(command) = [Any[first(command),command,false,0,0,0,0,0,0,0]]
    unit = ManagedRT.agent_unit_name(config.worker_id)
    cgroup = "/user.slice/user-1000.slice/user@1000.service/app.slice/" * unit
    properties = Dict{String,Any}("Id"=>record("s",unit),"InvocationID"=>record("ay",collect(1:16)),
        "ActiveState"=>record("s","active"),"MainPID"=>record("u",getpid()),"Type"=>record("s","exec"),
        "ExitType"=>record("s","main"),"ControlGroup"=>record("s",cgroup),"KillMode"=>record("s","control-group"),
        "KillSignal"=>record("i",2),"SendSIGKILL"=>record("b",true),"FinalKillSignal"=>record("i",9),
        "RemainAfterExit"=>record("b",false),"TimeoutStopFailureMode"=>record("s","terminate"),
        "ExecStop"=>record("a(sasbttttuii)",[]),
        "TimeoutStopUSec"=>record("t",60_000_000),"Restart"=>record("s","no"),
        "ExecStart"=>record("a(sasbttttuii)",exec_record(ManagedRT.agent_start_command(path))),
        "ExecStopPost"=>record("a(sasbttttuii)",exec_record(ManagedRT.agent_recovery_command(config))))
    object = Ref("/org/freedesktop/systemd1/unit/owned_2eservice")
    history = Vector{String}[]
    function invoke(args)
        push!(history,copy(args))
        output = if "GetUnit" in args
            ManagedRT.JSON3.write(record("o",[object[]]))
        else
            interface = findfirst(a->a in ("org.freedesktop.systemd1.Unit","org.freedesktop.systemd1.Service"),args)
            join((ManagedRT.JSON3.write(properties[key]) for key in args[interface+1:end]),'\n')
        end
        CommandResult(0,output,"")
    end
    options = (;invoke,which=_->"/fixture/busctl",read_cgroup=()->"0::" * cgroup * "\n")
    return (;properties,object,history,options)
end

@testset "managed agent requires actual typed lifetime/cleanup evidence" begin
    mktempdir() do directory
        path,config,_ = managed_config(directory)
        runner = CommandRunner()
        f = managed_bus_fixture(path,config)
        try
            identity = verify_managed_agent(runner,path,config;f.options...)
            @test identity.unit == "lcm-agent-worker-a.service"
            @test identity.journal_root == joinpath(config.scratch_root,"resources")
            @test verify_managed_agent(runner,path,config;f.options...,previous=identity).invocation == identity.invocation
            @test !occursin(directory,repr(MIME"text/plain"(),identity))
            @test !ispath(config.scratch_root)
            @test !ManagedRT.managed_agent_stopping(runner,f.object[],identity;invoke=f.options.invoke,which=f.options.which)
            f.properties["ActiveState"]["data"] = "deactivating"
            @test ManagedRT.managed_agent_stopping(runner,f.object[],identity;invoke=f.options.invoke,which=f.options.which)
            @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...)
            f.properties["ActiveState"]["data"] = "inactive"
            @test_throws CommandFailure ManagedRT.managed_agent_stopping(runner,f.object[],identity;invoke=f.options.invoke,which=f.options.which)
            f.properties["ActiveState"]["data"] = "active"
            f.properties["InvocationID"]["data"] = collect(2:17)
            @test_throws CommandFailure ManagedRT.managed_agent_stopping(runner,f.object[],identity;invoke=f.options.invoke,which=f.options.which)
            f.properties["InvocationID"]["data"] = collect(1:16)
            for (key,value) in (("MainPID",getpid()+1),("MainPID",true),("Id","unrelated.service"),
                    ("ActiveState","inactive"),("Type","simple"),("ExitType","cgroup"),("KillMode","process"),
                    ("KillSignal",15),("SendSIGKILL",false),("FinalKillSignal",2),("TimeoutStopUSec",typemax(UInt64)),
                    ("Restart","always"),("ControlGroup","/host/unrelated"),("InvocationID",zeros(Int,16)),
                    ("InvocationID",[1,2]),("ExecStopPost",[]),("ExecStart",[]),
                    ("RemainAfterExit",true),("TimeoutStopFailureMode","abort"),("ExecStop",["unowned-command"]))
                previous = f.properties[key]["data"]
                f.properties[key]["data"] = value
                @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...)
                f.properties[key]["data"] = previous
            end
            f.properties["InvocationID"]["data"] = collect(2:17)
            @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...,previous=identity)
            f.properties["InvocationID"]["data"] = collect(1:16)
            f.properties["ExecStopPost"]["data"][1][3] = true
            @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...)
            f.properties["ExecStopPost"]["data"][1][3] = false
            f.properties["MainPID"]["type"] = "s"
            @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...)
            f.properties["MainPID"]["type"] = "u"
            @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...,read_cgroup=()->"0::/other")
            f.object[] = "/org/freedesktop/systemd1/../../unrelated"
            @test_throws CommandFailure verify_managed_agent(runner,path,config;f.options...)
            @test !occursin(directory,sprint(showerror,CommandFailure(:managed_agent_unverified)))
            @test all(args->!("start" in args || "stop" in args),f.history)
        finally
            close(runner)
        end
    end
end

@testset "unit rendering is inert and captures the original recovery target" begin
    mktempdir() do directory
        path,config,data = managed_config(directory)
        unit = agent_service_unit(path)
        @test occursin("Type=exec",unit) && occursin("KillMode=control-group",unit)
        @test occursin("TimeoutStopSec=60",unit) && occursin("Restart=no",unit)
        @test occursin("recover-agent",unit) && occursin("--journal",unit)
        @test occursin(ManagedRT.agent_journal_root(config),unit)
        @test !occursin("managed-fixture-private-secret",unit)
        @test !occursin("broker.password",unit)
        @test !ispath(config.scratch_root)
        @test ManagedRT.systemd_quote("a b%\$\"\\") == "\"a b%%\$\$\\\"\\\\\""
        @test_throws ArgumentError ManagedRT.systemd_quote("bad\nargument")
        original = only(filter(line->startswith(line,"ExecStopPost="),split(unit,'\n')))
        data["agent"]["scratch_root"] = joinpath(directory,"replacement")
        open(path,"w") do io; TOML.print(io,data); end
        @test original != only(filter(line->startswith(line,"ExecStopPost="),split(agent_service_unit(path),'\n')))
        @test occursin(ManagedRT.agent_journal_root(config),original)
        @test !occursin("replacement",original)
    end
end

@testset "agent recovery never races a live owner or restores authority" begin
    mktempdir() do directory
        root = joinpath(directory,"journal")
        @test recover_agent_resources!(root,"worker-a") === nothing
        @test !ispath(root)
        journal = ResourceJournal(root,"worker-a")
        @test_throws ArgumentError recover_agent_resources!(root,"worker-a")
        close(journal)
        @test recover_agent_resources!(root,"worker-a") === nothing
        @test_throws ArgumentError recover_agent_resources!(root,"other-worker")
        journal = ResourceJournal(root,"worker-a")
        @test isempty(resource_receipts(journal))
        close(journal)
        unknown = joinpath(directory,"unknown")
        mkdir(unknown)
        @test_throws ArgumentError recover_agent_resources!(unknown,"worker-a")
        @test isempty(readdir(unknown))
        @test_throws ArgumentError recover_agent_resources!(homedir(),"worker-a")
    end
end
