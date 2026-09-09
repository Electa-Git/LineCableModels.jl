using Test, UUIDs, LineCableModelsRuntime
const PolicyRT = LineCableModelsRuntime
policy_json(value) = PolicyRT.JSON3.read(PolicyRT.JSON3.write(value))

function policy_fixture(; engine=:podman, kind=:scientific)
    image_ref = "registry.invalid/lcm/$kind@sha256:" * repeat("a",64)
    profile = ProfileDefinition("example",image_ref,repeat("a",64);kind,isolation=:container,
        operations=kind == :scientific ? ("fixture.echo",) : (),
        budget=ResourceBudget(cpus=0.125,memory_bytes=1024^3,pids=128,scratch_bytes=256*1024^2+19))
    fence = PolicyRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",
        string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
    receipt = ResourceReceipt(uuid4(),uuid4(),fence,engine,repeat("d",64),nothing)
    policy = ContainerPolicy(profile,receipt)
    image = Dict{String,Any}("RepoDigests"=>[image_ref],"Id"=>"sha256:" * repeat("c",64),"Os"=>"linux",
        "Config"=>Dict{String,Any}("Env"=>["PATH=/usr/bin:/bin","JULIA_VERSION=1.12.7"],
            "Labels"=>Dict("org.linecablemodels.runtime-image"=>"1","org.linecablemodels.profile"=>"example",
                "org.linecablemodels.profile-kind"=>string(kind))))
    checked_image = PolicyRT.verify_container_image(policy,policy_json(image))
    env = merge(checked_image.environment,PolicyRT.container_environment(policy))
    object = Dict{String,Any}("Id"=>repeat("b",64),"Name"=>"/" * resource_name(receipt),
        "Image"=>image["Id"],"State"=>Dict{String,Any}("Running"=>false,"Pid"=>0),"Mounts"=>[],
        "Config"=>Dict{String,Any}("Labels"=>resource_labels(receipt),"User"=>"1000:1000", "WorkingDir"=>"/tmp",
            "OpenStdin"=>true,"Hostname"=>"lcm-executor","Tty"=>kind == :terminal,"Entrypoint"=>["/usr/local/julia/bin/julia"],
            "Cmd"=>PolicyRT.container_julia_arguments(policy),"Env"=>["$k=$v" for (k,v) in env]),
        "HostConfig"=>Dict{String,Any}("ReadonlyRootfs"=>true,"Privileged"=>false,"NetworkMode"=>"none",
            "IpcMode"=>"private","CgroupnsMode"=>"private","PidMode"=>engine == :podman ? "private" : "",
            "UTSMode"=>engine == :podman ? "private" : "","Memory"=>1024^3,"MemorySwap"=>1024^3,
            "PidsLimit"=>128,"CpuPeriod"=>100000,"CpuQuota"=>12500,"ShmSize"=>65536,
            "CapDrop"=>["ALL"],"SecurityOpt"=>["no-new-privileges"],"LogConfig"=>Dict("Type"=>"none"),
            "RestartPolicy"=>Dict("Name"=>"no"),"Tmpfs"=>Dict("/tmp"=>PolicyRT.container_tmpfs_options(policy)),
            "Ulimits"=>[Dict("Name"=>k,"Soft"=>v,"Hard"=>v) for (k,v) in (("core",0),("msgqueue",0),("nofile",1024))]))
    if engine == :podman
        object["EffectiveCaps"] = nothing
        object["BoundingCaps"] = nothing
        object["HostConfig"]["CapDrop"] = ["CAP_CHOWN","CAP_SETUID"]
        object["HostConfig"]["CgroupMode"] = pop!(object["HostConfig"],"CgroupnsMode")
        object["HostConfig"]["Tmpfs"]["/tmp"] = replace(PolicyRT.container_tmpfs_options(policy),"notmpcopyup"=>"rprivate")
        for u in object["HostConfig"]["Ulimits"]
            u["Name"] = "RLIMIT_" * uppercase(u["Name"])
        end
    end
    return (; policy, profile, fence, receipt, image, checked_image, object)
end

@testset "one container policy for both engines and both process kinds" begin
    for engine in (:podman,:docker), kind in (:scientific,:terminal)
        f = policy_fixture(;engine,kind)
        args = container_create_arguments(f.policy)
        @test PolicyRT.verify_created_container(f.policy,policy_json(f.object),f.checked_image) === nothing
        @test args[1:2] == ["container","create"]
        @test "--cpu-quota=12500" in args
        @test "--memory-swap=1073741824" in args
        @test "--network=none" in args && "--read-only" in args
        @test "--ulimit=msgqueue=0:0" in args
        @test ("--tty" in args) == (kind == :terminal)
        @test ("--read-only-tmpfs=false" in args) == (engine == :podman)
        @test ("--http-proxy=false" in args) == (engine == :podman)
        @test (engine == :podman ? "--pid=private" : "--pid=") in args
        @test PolicyRT.container_tmpfs_bytes(f.policy) + 65536 == 256*1024^2
        @test !any(x->occursin("alice",x),args)
        @test last(args) == (kind == :scientific ? "/opt/lcm/scientific.jl" : "--color=yes")
        @test !any(x->startswith(x,"--volume") || startswith(x,"--mount") || x in ("run","--privileged"),args)
    end
end

@testset "image/config drift never authorizes a container start" begin
    for engine in (:podman,:docker)
        f = policy_fixture(;engine)
        for (key,value) in (("RepoDigests",["unapproved"]),("Id","short"),("Os","windows"))
            object = deepcopy(f.image); object[key] = value
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_container_image(f.policy,policy_json(object))
        end
        image = deepcopy(f.image); image["Id"] = repeat("c",64)
        @test PolicyRT.verify_container_image(f.policy,policy_json(image)).id == f.checked_image.id
        for id in (repeat("c",63),repeat("C",64),"sha512:" * repeat("c",64))
            image["Id"] = id
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_container_image(f.policy,policy_json(image))
        end
        for (key,value) in (("Volumes",Dict("/data"=>Dict())),("Healthcheck",Dict("Test"=>["CMD","probe"])),
                ("OnBuild",["RUN something"]),("Env",["AWS_SECRET_ACCESS_KEY=secret"]),("Env",["PATH=/bin","PATH=/wrong"]),
                ("Labels",Dict()),("Env",["*"]),("Env",["PROXY=private"]))
            object = deepcopy(f.image); object["Config"][key] = value
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_container_image(f.policy,policy_json(object))
        end
        for (key,value) in (("Privileged",true),("ReadonlyRootfs",false),("NetworkMode","host"),("PidMode","host"),
                ("IpcMode","shareable"),(engine == :podman ? "CgroupMode" : "CgroupnsMode","host"),("MemorySwap",-1),("CpuQuota",0),("PidsLimit",-1),
                ("ShmSize",1024^3),("CapAdd",["SYS_ADMIN"]),("CapDrop",[]),("Devices",[Dict("PathOnHost"=>"/dev/sda")]),
                ("Binds",["/home:/host"]),("GroupAdd",["0"]),("SecurityOpt",["seccomp=unconfined"]),
                ("Tmpfs",Dict("/tmp"=>"rw")),("Ulimits",[]),("LogConfig",Dict("Type"=>"json-file")),
                ("RestartPolicy",Dict("Name"=>"always")),("ReadonlyRootfs",1))
            object = deepcopy(f.object); object["HostConfig"][key] = value
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
        for (key,value) in (("User","0"),("Hostname","unexpected"),("Cmd",["-e","exit()"]),("Entrypoint",["sh"]),("Tty",true),
                ("Env",["AWS_SECRET_ACCESS_KEY=secret"]),("WorkingDir","/host"),("OpenStdin",1))
            object = deepcopy(f.object); object["Config"][key] = value
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
        for (key,value) in (("Name","unrelated"),("Image","sha256:" * repeat("f",64)),
                ("State",Dict("Running"=>true,"Pid"=>23)),("Mounts",[Dict("Type"=>"bind","Destination"=>"/tmp")]))
            object = deepcopy(f.object); object[key] = value
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
        object = deepcopy(f.object); push!(object["Config"]["Env"],"container=podman")
        if engine == :podman
            @test PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image) === nothing
            for key in ("EffectiveCaps","BoundingCaps")
                object = deepcopy(f.object); object[key] = ["CAP_CHOWN"]
                @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
                delete!(object,key)
                @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
            end
            object = deepcopy(f.object); object["HostConfig"]["Tmpfs"]["/tmp"] *= ",tmpcopyup"
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        else
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
        object = deepcopy(f.object); object["Config"]["Entrypoint"] = "/usr/local/julia/bin/julia"
        if engine == :podman
            @test PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image) === nothing
        else
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
        for entry in ("/usr/local/julia/bin/julia -e exit()", "sh", "", nothing,
                ["/usr/local/julia/bin/julia", "--interactive"])
            object = deepcopy(f.object); object["Config"]["Entrypoint"] = entry
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
        object = deepcopy(f.object)
        push!(object["HostConfig"]["Ulimits"],Dict("Name"=>"RLIMIT_NPROC","Soft"=>255428,"Hard"=>255428))
        if engine == :podman
            @test PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image) === nothing
            for (key,value) in (("Name","RLIMIT_MEMLOCK"),("Soft",-1),("Soft",true),("Hard",0),("Hard",typemax(Int64)))
                changed = deepcopy(object); changed["HostConfig"]["Ulimits"][end][key] = value
                @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(changed),f.checked_image)
            end
        else
            @test_throws PolicyRT.ExecutionCore.IsolationError PolicyRT.verify_created_container(f.policy,policy_json(object),f.checked_image)
        end
    end
end

@testset "checked acquisition records before create and keeps failed resources recoverable" begin
    for engine in (:podman,:docker)
        f = policy_fixture(;engine)
        mktempdir() do directory
            runner = CommandRunner()
            journal = ResourceJournal(joinpath(directory,"journal"),"worker-a")
            host = ContainerHostCheck(ContainerEngine(engine,string(engine),true),[string(engine)],true,())
            info = engine == :docker ? Dict("ID"=>"fixture-daemon","OSType"=>"linux") :
                Dict("host"=>Dict("os"=>"linux","serviceIsRemote"=>false),"store"=>Dict("graphRoot"=>"/owned/store"))
            options = (probe=args->(true,PolicyRT.JSON3.write(info)),machine_id=repeat("1",32),storage_identity=p->(p,"1","2"))
            created = Ref(false); changed = Ref(false); history = Vector{String}[]
            function invoke(args)
                push!(history,copy(args))
                if "image" in args
                    return CommandResult(0,PolicyRT.JSON3.write(f.image),"")
                elseif "create" in args
                    @test length(resource_receipts(journal)) == 1
                    @test only(resource_receipts(journal)).physical_id === nothing
                    created[] = true
                    return CommandResult(0,repeat("b",64),"")
                elseif "ls" in args
                    return CommandResult(0,"","")
                elseif "inspect" in args
                    receipt = only(resource_receipts(journal))
                    object = deepcopy(f.object)
                    object["Name"] = resource_name(receipt)
                    object["Config"]["Labels"] = resource_labels(receipt)
                    changed[] && (object["HostConfig"]["ReadonlyRootfs"] = false)
                    return CommandResult(0,PolicyRT.JSON3.write(object),"")
                end
                error("unexpected command")
            end
            try
                unavailable = ContainerHostCheck(host.engine,host.command,true,(:cpu_controller_missing,))
                @test_throws PolicyRT.ExecutionCore.IsolationError create_owned_container!(journal,runner,unavailable,f.profile,f.fence;options...,invoke)
                @test isempty(history) && isempty(resource_receipts(journal))
                receipt = create_owned_container!(journal,runner,host,f.profile,f.fence;options...,invoke)
                @test receipt.physical_id == repeat("b",64) && created[]
                @test !any(args->"start" in args || "run" in args || "pull" in args,history)
                @test_throws PolicyRT.ExecutionCore.IsolationError create_owned_container!(journal,runner,host,f.profile,f.fence;options...,invoke)
                @test count(args->"create" in args,history) == 1
                # Fixture-only forgetting: no actual container was allocated.
                forget_resource!(journal,receipt)
                changed[] = true
                @test_throws PolicyRT.ExecutionCore.IsolationError create_owned_container!(journal,runner,host,f.profile,f.fence;options...,invoke)
                @test only(resource_receipts(journal)).physical_id === nothing
                @test only(resource_receipts(journal)).fence == f.fence
            finally
                close(journal); close(runner)
            end
        end
    end
end
