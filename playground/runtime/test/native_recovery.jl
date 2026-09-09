using Test, UUIDs, LineCableModelsRuntime
const NativeRT = LineCableModelsRuntime

function with_native_fixture(action)
    mktempdir() do directory
        runner = CommandRunner()
        journal = ResourceJournal(joinpath(directory,"journal"),"worker-a")
        machine = Ref(repeat("1",32))
        owner = Ref(1000)
        units = Dict{String,Any}()
        jobs = Any[]
        history = Vector{String}[]
        mutations = Vector{String}[]
        clock = Ref(0.0)
        state = Dict(:list=>true,:properties=>true,:stop=>true,:group_empty=>true,
            :bad_list=>false,:bad_jobs=>false,:change_invocation=>false,:scope_after_stop=>false,:retain_record=>false)
        record(type,data) = Dict{String,Any}("type"=>type,"data"=>data)
        function invoke(arguments)
            push!(history,copy(arguments))
            payload = if "GetMachineId" in arguments
                record("s",[machine[]])
            elseif "GetConnectionUnixUser" in arguments
                record("u",[owner[]])
            elseif "ListUnitsByPatterns" in arguments
                state[:list] || return CommandResult(1,"","private bus failure")
                name = last(arguments)
                records = state[:bad_list] ? Any[[name]] : haskey(units,name) ?
                    Any[Any[name,"owned","loaded",units[name]["ActiveState"]["data"],"running","",units[name]["path"],0,"","/"]] : Any[]
                record("a(ssssssouso)",[records])
            elseif "ListJobs" in arguments
                record("a(usssoo)",[state[:bad_jobs] ? Any[[1,"bad"]] : jobs])
            elseif "get-property" in arguments
                state[:properties] || return CommandResult(1,"","private property failure")
                index = findfirst(==("get-property"),arguments)
                name = only([name for (name,unit) in units if unit["path"] == arguments[index+2]])
                keys = arguments[index+4:end]
                text = join((NativeRT.JSON3.write(units[name][key]) for key in keys),'\n')
                if state[:change_invocation] && "InvocationID" in keys
                    units[name]["InvocationID"]["data"] = collect(2:17)
                end
                return CommandResult(0,text,"")
            else
                error("Unexpected fixture query.")
            end
            return CommandResult(0,NativeRT.JSON3.write(payload),"")
        end
        options = (;uid=1000,machine_id=repeat("1",32),socket_check=_->nothing,invoke,which=_->"/fixture/busctl")
        bus_options = (;uid=1000,invoke,which=_->"/fixture/busctl")
        scope = native_scope(runner;options...)
        function add!(;bound=false,physical=true)
            fence = NativeRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",
                string(uuid4()),string(uuid4()),"fixture","1.0.0",repeat("a",64),1)
            receipt = reserve_resource!(journal,fence,:native,scope)
            name = resource_name(receipt)
            if physical
                units[name] = Dict{String,Any}("path"=>"/org/freedesktop/systemd1/unit/owned_" * replace(string(receipt.id),"-"=>"_"),
                    "Id"=>record("s",name),"Description"=>record("s",native_resource_description(receipt)),
                    "InvocationID"=>record("ay",collect(1:16)),"Transient"=>record("b",true),"ActiveState"=>record("s","active"),
                    "MainPID"=>record("u",2000),"ControlGroup"=>record("s",NativeRT.native_cgroup(receipt,1000)),
                    "Slice"=>record("s","app.slice"),"Type"=>record("s","exec"),"ExitType"=>record("s","main"),
                    "RemainAfterExit"=>record("b",false),"Restart"=>record("s","no"),"KillMode"=>record("s","control-group"),
                    "KillSignal"=>record("i",2),"SendSIGKILL"=>record("b",true),"FinalKillSignal"=>record("i",9),
                    "TimeoutStopUSec"=>record("t",2_000_000),"TimeoutStopFailureMode"=>record("s","terminate"),
                    "ExecStop"=>record("a(sasbttttuii)",[]),"ExecStopPost"=>record("a(sasbttttuii)",[]))
            end
            return bound ? bind_resource!(journal,receipt,bytes2hex(UInt8.(1:16))) : receipt
        end
        function control(arguments)
            push!(mutations,copy(arguments))
            name = last(arguments)
            receipt = only(filter(r->resource_name(r)==name,resource_receipts(journal)))
            pending = haskey(units,name) && all(iszero,units[name]["InvocationID"]["data"]) &&
                units[name]["MainPID"]["data"] == 0 && units[name]["ActiveState"]["data"] == "inactive"
            receipt.physical_id !== nothing || pending || error("Mutation before durable invocation binding.")
            state[:stop] || return false
            if state[:retain_record]
                units[name]["ActiveState"]["data"] = "inactive"
                units[name]["MainPID"]["data"] = 0
                units[name]["ControlGroup"]["data"] = ""
            else
                delete!(units,name)
            end
            filter!(job->job[2] != name,jobs)
            state[:scope_after_stop] && (machine[] = repeat("2",32))
            return true
        end
        cleanup_options = (;options...,control,group_empty=_->state[:group_empty],clock=()->clock[],
            pause=seconds->(clock[]+=seconds),timeout_seconds=0.2)
        try
            action((;runner,journal,units,jobs,state,machine,owner,history,mutations,scope,add!,
                options,bus_options,cleanup_options,directory))
        finally
            close(journal);close(runner)
        end
    end
end

@testset "native queued pre-activation intent can be canceled without a PID" begin
    with_native_fixture() do f
        receipt = f.add!()
        unit = f.units[resource_name(receipt)]
        unit["InvocationID"]["data"] = Int[]
        unit["ActiveState"]["data"] = "inactive"
        unit["MainPID"]["data"] = 0
        unit["ControlGroup"]["data"] = ""
        push!(f.jobs,Any[9,resource_name(receipt),"start","waiting","/owned","/job"])
        @test inspect_native_unit(f.runner,receipt;f.bus_options...).invocation === nothing
        unit["InvocationID"]["data"] = zeros(Int,16)
        @test inspect_native_unit(f.runner,receipt;f.bus_options...).invocation === nothing
        unit["InvocationID"]["data"] = Int[]
        for (key,changed) in (("MainPID",1),("ControlGroup",NativeRT.native_cgroup(receipt,1000)),("ActiveState","active"))
            saved = unit[key]["data"];unit[key]["data"] = changed
            @test_throws CommandFailure inspect_native_unit(f.runner,receipt;f.bus_options...)
            unit[key]["data"] = saved
        end
        @test remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test f.mutations == [["stop","--no-block",resource_name(receipt)]]
        @test isempty(f.jobs) && isempty(f.units) && isempty(resource_receipts(f.journal))
    end
end

@testset "inactive native metadata need not erase another unit's references" begin
    with_native_fixture() do f
        receipt = f.add!(;bound=true)
        f.state[:retain_record] = true
        @test remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test haskey(f.units,resource_name(receipt))
        @test f.units[resource_name(receipt)]["MainPID"]["data"] == 0
        @test isempty(resource_receipts(f.journal))
        @test f.mutations == [["stop","--no-block",resource_name(receipt)]]
        @test remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test f.mutations == [["stop","--no-block",resource_name(receipt)]]
    end
end

@testset "native cleanup scope binds local socket, machine and user" begin
    with_native_fixture() do f
        @test occursin(r"^[a-f0-9]{64}$",native_scope(f.runner;f.options...))
        @test native_scope(f.runner;f.options...) == f.scope
        @test native_scope(f.runner;f.options...,machine_id=strip(" " * repeat("1",32) * " ")) == f.scope
        @test all(args->args[2]=="--address=unix:path=/run/user/1000/bus",f.history)
        f.machine[] = repeat("2",32)
        @test_throws CommandFailure native_scope(f.runner;f.options...)
        f.machine[] = repeat("1",32);f.owner[] = 2000
        @test_throws CommandFailure native_scope(f.runner;f.options...)
        @test_throws CommandFailure native_scope(f.runner;f.options...,uid=0)
        @test_throws CommandFailure native_scope(f.runner;f.options...,socket_check=_->false)
        @test isempty(f.mutations) && isempty(resource_receipts(f.journal))
    end
end

@testset "native crash-gap binding and exact independent cleanup" begin
    with_native_fixture() do f
        first,second = f.add!(),f.add!(;bound=true)
        identity = inspect_native_unit(f.runner,first;f.bus_options...)
        @test identity.unit == resource_name(first) && identity.pid == 2000
        @test !occursin("user.slice",repr(MIME"text/plain"(),identity))
        @test remove_owned_native!(f.journal,f.runner,first;f.cleanup_options...)
        @test length(f.units) == 1 && haskey(f.units,resource_name(second))
        @test only(resource_receipts(f.journal)).id == second.id
        @test f.mutations == [["stop","--no-block",resource_name(first)]]
        @test remove_owned_native!(f.journal,f.runner,first;f.cleanup_options...)
        @test length(f.mutations) == 1
        @test recover_native!(f.journal,f.runner;f.cleanup_options...) === nothing
        @test isempty(f.units) && isempty(resource_receipts(f.journal))
    end
end

@testset "native identity and cleanup policy cannot drift" begin
    with_native_fixture() do f
        receipt = f.add!(;bound=true)
        unit = f.units[resource_name(receipt)]
        for (key,changed) in (("Description","unrelated service"),("Id","unrelated.service"),("Transient",false),
                ("InvocationID",collect(2:17)),("InvocationID",zeros(Int,16)),("InvocationID",[1]),
                ("MainPID",true),("ControlGroup","/user.slice/unrelated"),("Slice","other.slice"),
                ("Type","simple"),("ExitType","cgroup"),("RemainAfterExit",true),("Restart","always"),
                ("KillMode","process"),("KillSignal",15),("SendSIGKILL",false),("FinalKillSignal",2),
                ("TimeoutStopUSec",typemax(UInt64)),("TimeoutStopFailureMode","abort"),
                ("ExecStop",["foreign"]),("ExecStopPost",["foreign"]),("ActiveState","unknown"))
            saved = unit[key]["data"];unit[key]["data"] = changed
            @test_throws CommandFailure remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
            unit[key]["data"] = saved
        end
        f.state[:change_invocation] = true
        @test_throws CommandFailure remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test isempty(f.mutations) && length(resource_receipts(f.journal)) == 1
    end
end

@testset "native completed failed record is reset only after kernel/job cleanup" begin
    with_native_fixture() do f
        receipt = f.add!(;bound=true)
        unit = f.units[resource_name(receipt)]
        unit["ActiveState"]["data"] = "failed"
        unit["MainPID"]["data"] = 0
        unit["ControlGroup"]["data"] = ""
        push!(f.jobs,Any[7,"unrelated.service","start","waiting","/unrelated","/job"])
        @test remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test f.mutations == [["reset-failed",resource_name(receipt)]]
        @test length(f.jobs) == 1 && f.jobs[1][2] == "unrelated.service"
        @test isempty(resource_receipts(f.journal))
    end
end

@testset "native inspection errors and pending kernel/jobs retain ownership" begin
    with_native_fixture() do f
        receipt = f.add!()
        for key in (:list,:properties,:stop)
            f.state[key] = false
            if key == :stop
                @test !remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
            else
                @test_throws CommandFailure remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
            end
            @test length(resource_receipts(f.journal)) == 1
            f.state[key] = true
        end
        empty!(f.units)
        f.state[:bad_list] = true
        @test_throws CommandFailure remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        f.state[:bad_list] = false;f.state[:bad_jobs] = true
        @test_throws CommandFailure remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        f.state[:bad_jobs] = false
        push!(f.jobs,Any[1,resource_name(receipt),"start","waiting","/unit","/job"])
        @test !remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        empty!(f.jobs);f.state[:group_empty] = false
        @test !remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test length(resource_receipts(f.journal)) == 1
        f.state[:group_empty] = true
        @test remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test isempty(resource_receipts(f.journal))
    end
end

@testset "native recovery keeps other backends and rejects changed scopes" begin
    with_native_fixture() do f
        receipt = f.add!(;bound=true)
        absent = f.add!(;physical=false)
        other = reserve_resource!(f.journal,NativeRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"alice","other",
            "worker-a",string(uuid4()),string(uuid4()),"fixture","1.0.0",repeat("a",64),1),:podman,repeat("c",64))
        f.state[:scope_after_stop] = true
        @test_throws CommandFailure remove_owned_native!(f.journal,f.runner,receipt;f.cleanup_options...)
        @test length(resource_receipts(f.journal)) == 3
        f.state[:scope_after_stop] = false;f.machine[] = repeat("1",32)
        @test recover_native!(f.journal,f.runner;f.cleanup_options...) === nothing
        @test only(resource_receipts(f.journal)).id == other.id
    end
end
