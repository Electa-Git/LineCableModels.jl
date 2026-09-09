using Test, UUIDs, LineCableModelsRuntime
const PhysicalRT = LineCableModelsRuntime

# Opt-in physical gate: uses installed approved images and the real policy,
# journal, entry guard and PTY transport. It does not claim broker/lease or
# managed-agent integration, which have separate end-to-end gates.
function physical_terminal_gate()
    engine = get(ENV,"LCM_TEST_CONTAINER_RUNTIME","")
    engine in ("docker","podman") || error("Select an explicit test container engine.")
    image = get(ENV,"LCM_TEST_TERMINAL_IMAGE","")
    occursin(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*@sha256:[a-f0-9]{64}$",image) ||
        error("Supply an installed digest-pinned terminal image; this gate never pulls.")
    directory = mktempdir(;prefix="lcm-physical-terminal-",cleanup=false)
    runner = CommandRunner(;timeout_seconds=20)
    journal = ResourceJournal(joinpath(directory,"journal"),"physical-terminal";capacity=2)
    handles = Any[]
    host = check_container_host(runner;requested=engine)
    isempty(host.failures) || error("Mandatory host controls unavailable.")
    profile = ProfileDefinition("julia-terminal",image,last(split(image,"@sha256:"));
        kind=:terminal,isolation=:container,
        budget=ResourceBudget(cpus=0.5,memory_bytes=512*1024^2,pids=64,scratch_bytes=8*1024^2))
    command(args) = PhysicalRT.scoped_container_command(runner,host,args)
    function checked(args)
        result = command(args)
        result.exitcode == 0 || error("Physical inspection failed: " * result.diagnostic)
        return result.output
    end
    inventory() = sort(split(strip(checked(["container","ls","--all","--no-trunc","--format","{{.ID}}"])),
        '\n';keepempty=false))
    baseline = inventory()
    object(receipt) = PhysicalRT.JSON3.read(checked(["container","inspect","--format","{{json .}}",receipt.physical_id]))
    function text(process)
        lock(process.lock) do
            String(PhysicalRT.read_terminal(process.output,0,length(process.output.data)).bytes)
        end
    end
    function expect(process,needle;seconds=30)
        result = timedwait(()->occursin(needle,text(process)) || process.cleanup_complete,seconds;pollint=0.025)
        if result != :ok || !occursin(needle,text(process))
            path = joinpath(directory,"last-terminal.txt")
            write(path,text(process))
            chmod(path,0o600)
            error("Physical terminal did not produce expected evidence; diagnostics: " * directory)
        end
        return nothing
    end
    function send(process,code)
        PhysicalRT.write_terminal!(process,collect(codeunits(code * "\r")))
    end
    function acquire(label;flood=false)
        fence = PhysicalRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"physical-owner",label,
            "physical-terminal",string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
        receipt = PhysicalRT.create_owned_container!(journal,runner,host,profile,fence)
        process = PhysicalRT.TerminalProcess(;limits=PhysicalRT.TerminalIOLimits(
            lifetime_seconds=300,output_bytes=flood ? 8192 : 1024^2,
            output_bytes_per_second=flood ? 8192 : 1024^2))
        handle = (;receipt,process,label)
        push!(handles,handle)
        PhysicalRT.start_terminal!(process,PhysicalRT.terminal_attach_command(host,receipt);columns=100,rows=30)
        expect(process,"\x1elcm-terminal-ready:" * string(receipt.id) * "\x1f";seconds=60)
        policy = ContainerPolicy(profile,receipt)
        approved = PhysicalRT.verify_container_image(policy,PhysicalRT.JSON3.read(
            checked(["image","inspect","--format","{{json .}}",image])))
        @test PhysicalRT.verify_created_container(policy,object(receipt),approved;running=true) === nothing
        send(process,"println(\"KERNEL_\", \"READY:\", LCM_CONTAINER_ISOLATION.cpus, \":\", LCM_CONTAINER_ISOLATION.memory_bytes)")
        expect(process,"KERNEL_READY:0.5:536870912")
        return handle
    end
    function retire(handle)
        close(handle.process)
        @test remove_owned_container!(journal,runner,host,handle.receipt)
        filter!(h->h !== handle,handles)
    end
    function alive(handle,round)
        send(handle.process,"println(\"SURVIVOR_\", \"ROUND:\", $round)")
        expect(handle.process,"SURVIVOR_ROUND:$round")
        @test object(handle.receipt).State.Running === true
    end
    try
        survivor = acquire("survivor")
        send(survivor.process,"private_counter=41; println(\"RESULT_\", \"VALUE:\", private_counter+1)")
        expect(survivor.process,"RESULT_VALUE:42")
        @test PhysicalRT.resize_terminal!(survivor.process,120,40) === nothing
        other = acquire("separate")
        send(other.process,"println(\"SEPARATE_\", \"STATE:\", isdefined(Main,:private_counter)); println(\"SOCKET_\", \"PRESENT:\", ispath(\"/var/run/docker.sock\"))")
        expect(other.process,"SEPARATE_STATE:false")
        expect(other.process,"SOCKET_PRESENT:false")
        retire(other)
        alive(survivor,1)

        flood = acquire("flood";flood=true)
        send(flood.process,"for i in 1:100000; println(repeat(\"x\",4096)); end")
        @test timedwait(()->flood.process.cleanup_complete,30;pollint=0.025) == :ok
        @test flood.process.failure == :output_rate_limit
        @test 0 < length(flood.process.output.data) <= 8192
        retire(flood)
        alive(survivor,2)

        memory = acquire("memory")
        # Keep the parent alive to read the kernel's cgroup event counters.
        # Rootless engine OOMKilled metadata is not the effective-limit evidence.
        allocation = ["/usr/local/julia/bin/julia","--startup-file=no","--history-file=no",
            "-e","held=fill(UInt8(17),1024^3); println(\"ALLOCATION_UNEXPECTEDLY_FINISHED\")"]
        send(memory.process,"victim=run(pipeline(ignorestatus(Cmd($(repr(allocation))));stdout=devnull,stderr=devnull)); println(\"MEMORY_\",\"EXIT:\",victim.exitcode,\":\",victim.termsignal); println(\"MEMORY_\",\"EVENTS:\",replace(strip(read(\"/sys/fs/cgroup/memory.events\",String)),Char(10)=>';'))")
        expect(memory.process,"MEMORY_EVENTS:";seconds=45)
        evidence = text(memory.process)
        oom = match(r"MEMORY_EVENTS:[^\r\n]*oom_kill (\d+)",evidence)
        @test oom !== nothing && parse(Int,oom.captures[1]) > 0
        @test occursin("MEMORY_EXIT:0:9",evidence) || occursin("MEMORY_EXIT:137:0",evidence)
        retire(memory)
        alive(survivor,3)

        tasks = acquire("tasks")
        send(tasks.process,"children=Any[]; for i in 1:100; try push!(children,run(pipeline(`/bin/sleep 30`;stdout=devnull,stderr=devnull);wait=false)) catch; break; end; end; println(\"PID_\", \"EVENTS:\", strip(read(\"/sys/fs/cgroup/pids.events\",String))); println(\"PID_\", \"CHILDREN:\",length(children))")
        expect(tasks.process,"PID_EVENTS:max ";seconds=30)
        count = match(r"PID_EVENTS:max (\d+)",text(tasks.process))
        @test count !== nothing && parse(Int,count.captures[1]) > 0
        retire(tasks)
        alive(survivor,4)
        retire(survivor)
        @test isempty(resource_receipts(journal))
        @test inventory() == baseline
    finally
        for handle in handles
            try close(handle.process) catch end
        end
        try
            recover_containers!(journal,runner,host)
            isempty(resource_receipts(journal)) || error("Physical cleanup still owns resources: " * directory)
        finally
            close(journal); close(runner)
        end
        println("Physical terminal diagnostics: ",directory)
    end
end

@testset "real container terminal guard, resource limits and exact cleanup" begin
    physical_terminal_gate()
end
