using Test, UUIDs, LineCableModelsRuntime
const TerminalDriverRT = LineCableModelsRuntime

# Explicit synthetic admission fixture. No image, cgroup, engine or managed
# service proof is supplied, so real acquisition must fail before any command.
function terminal_partition_fixture(f)
    driver=f.driver
    terminal=ProfileDefinition("terminal","registry.invalid/lcm@sha256:"*repeat("b",64),repeat("b",64);
        kind=:terminal,isolation=:container)
    register!(driver.profiles,terminal)
    driver.config=AgentConfig(driver.config.worker_id,driver.config.endpoint,driver.profiles,
        driver.config.scratch_root;capacity=2,container_runtime=:podman)
    fence=TerminalDriverRT.Protocol.AssignmentFence((name==:lease_id ? string(uuid4()) :
        name==:role ? "terminal" : name==:profile_id ? terminal.id :
        name==:fingerprint ? terminal.fingerprint : getfield(f.fence,name)
        for name in fieldnames(typeof(f.fence)))...)
    driver.recovered=true
    return (;terminal,fence,science=ManagedScientificView(driver),terminals=ManagedTerminalView(driver))
end

@testset "agent root composes partitions without duplicating lease authority" begin
    with_container_driver_fixture() do f
        p=terminal_partition_fixture(f)
        resources=ManagedAgentResources(f.driver)
        agent=AgentService(f.driver.config,resources)
        try
            @test TerminalDriverRT.RequiredInterfaces.check_interface_implemented(AbstractAgentResources,ManagedAgentResources)
            @test installed_profiles(resources) === f.driver.profiles
            @test resources.science.driver.parent === resources.terminals.driver.parent === f.driver
            @test resources.science.ledger === resources.terminals.ledger === agent.ledger
            @test agent.science.resources === agent.jobs.resources === resources.science
            @test agent.science.ledger === agent.jobs.ledger === agent.ledger
            @test agent.ledger.capacity==2
            @test recover_owned!(resources) === nothing
            @test resources.science.recovered && resources.terminals.recovered
            @test_throws CommandFailure executor_for!(resources.science.driver,f.profile,f.fence)
            @test_throws CommandFailure terminal_for!(resources.terminals.driver,p.terminal,p.fence)
            @test length(f.driver.handles)==2
            @test release_owned!(resources,p.fence)
            @test haskey(f.driver.handles,f.fence.lease_id) && length(f.driver.handles)==1
            @test release_owned!(resources,f.fence)
            @test isempty(f.driver.handles) && !f.driver.closed
        finally
            close(agent)
        end
        @test agent.cleanup_complete && f.driver.cleanup_complete
        @test resources.science.closed && resources.terminals.closed
        @test close(resources) === nothing
    end
end

@testset "managed partitions share physical capacity and exact cleanup" begin
    @test ManagedResourceDriver === ManagedScientificDriver
    @test TerminalDriverRT.RequiredInterfaces.check_interface_implemented(AbstractScientificDriver,ManagedScientificView)
    @test TerminalDriverRT.RequiredInterfaces.check_interface_implemented(AbstractTerminalDriver,ManagedTerminalView)
    with_container_driver_fixture() do f
        p=terminal_partition_fixture(f)
        @test only(values(installed_profiles(p.science).definitions)) === f.profile
        @test only(values(installed_profiles(p.terminals).definitions)) === p.terminal
        @test isempty(f.driver.handles) && isempty(resource_receipts(f.driver.journal))
        @test recover_owned!(p.science) === nothing
        @test recover_owned!(p.terminals) === nothing
        @test_throws ArgumentError executor_for!(p.science,p.terminal,p.fence)
        @test_throws ArgumentError terminal_for!(p.terminals,f.profile,f.fence)
        @test isempty(f.driver.handles)
        @test_throws CommandFailure executor_for!(p.science,f.profile,f.fence)
        @test_throws CommandFailure terminal_for!(p.terminals,p.terminal,p.fence)
        @test length(f.driver.handles)==2
        @test isempty(resource_receipts(f.driver.journal)) && isempty(f.driver.runner.active)
        replacement=TerminalDriverRT.Protocol.AssignmentFence((name==:lease_id ? string(uuid4()) :
            getfield(p.fence,name) for name in fieldnames(typeof(p.fence)))...)
        @test_throws CapacityUnavailable terminal_for!(p.terminals,p.terminal,replacement)
        @test_throws ArgumentError release_owned!(p.science,p.fence)
        @test_throws ArgumentError release_owned!(p.terminals,f.fence)
        @test length(f.driver.handles)==2
        # A passive retained PTY proves partition retirement dispatch without
        # bypassing physical launch guards or touching an existing container.
        process=TerminalDriverRT.TerminalProcess()
        f.driver.handles[p.fence.lease_id].supervisor=process
        @test terminal_for!(p.terminals,p.terminal,p.fence) === process
        @test close(p.terminals) === nothing
        @test process.cleanup_complete && length(f.driver.handles)==1
        @test haskey(f.driver.handles,f.fence.lease_id) && !f.driver.closed
        @test_throws ArgumentError ResourceJournal(f.driver.journal.root,f.driver.config.worker_id)
        @test_throws ArgumentError terminal_for!(p.terminals,p.terminal,p.fence)
        @test_throws ArgumentError recover_owned!(p.terminals)
        @test close(p.terminals) === nothing
        @test close(p.science) === nothing && isempty(f.driver.handles)
        @test !f.driver.closed && close(f.driver) === nothing
        @test f.driver.cleanup_complete
    end
end

@testset "partition shutdown joins admitted acquisition and leaves its sibling live" begin
    with_container_driver_fixture() do f
        p=terminal_partition_fixture(f)
        entered,continue_acquisition=Channel{Nothing}(1),Channel{Nothing}(1)
        acquisition=@async TerminalDriverRT.with_managed_view(p.terminals) do
            put!(entered,nothing)
            take!(continue_acquisition)
            TerminalDriverRT.managed_handle!(f.driver,p.fence)
        end
        take!(entered)
        shutdown=@async close(p.terminals)
        @test timedwait(()->p.terminals.lifecycle.closed,2)==:ok
        @test !istaskdone(shutdown)
        @test_throws ArgumentError terminal_for!(p.terminals,p.terminal,p.fence)
        @test_throws CommandFailure executor_for!(p.science,f.profile,f.fence)
        @test haskey(f.driver.handles,f.fence.lease_id)
        put!(continue_acquisition,nothing)
        @test fetch(acquisition).fence==p.fence
        @test fetch(shutdown) === nothing
        @test p.terminals.lifecycle.active==0
        @test !haskey(f.driver.handles,p.fence.lease_id)
        @test haskey(f.driver.handles,f.fence.lease_id)
        @test close(p.science) === nothing
    end
end

@testset "terminal attachment names only the owned physical container" begin
    for engine in (:podman,:docker)
        f=policy_fixture(;engine,kind=:terminal)
        host=ContainerHostCheck(ContainerEngine(engine,string(engine),true),["/fixed/"*string(engine)],true,())
        receipt=ResourceReceipt(f.receipt.id,f.receipt.journal_id,f.fence,engine,f.receipt.scope,repeat("b",64))
        command=TerminalDriverRT.terminal_attach_command(host,receipt)
        @test collect(command)==[host.command;"container";"start";"--attach";"--interactive";
            "--detach-keys=";receipt.physical_id]
        @test command.env !== nothing
        @test !any(x->occursin("alice",x),collect(command))
        other=engine==:podman ? :docker : :podman
        @test_throws ArgumentError TerminalDriverRT.terminal_attach_command(
            ContainerHostCheck(ContainerEngine(other,string(other),true),host.command,true,()),receipt)
        @test_throws ArgumentError TerminalDriverRT.terminal_attach_command(host,f.receipt)
    end
end
