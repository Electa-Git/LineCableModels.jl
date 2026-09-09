using Test, UUIDs, LineCableModelsRuntime
const DriverRT = LineCableModelsRuntime

function with_container_driver_fixture(action)
    mktempdir() do directory
        path,config,_ = managed_config(directory)
        runner = CommandRunner()
        journal = ResourceJournal(DriverRT.agent_journal_root(config),config.worker_id)
        identity = ManagedAgentIdentity(DriverRT.agent_unit_name(config.worker_id),repeat("a",32),
            "/fixture/invalid-cgroup",path,journal.root,config.worker_id)
        driver = ContainerScientificDriver(config,identity,config.profiles,Dict{String,Symbol}(),
            runner,journal,nothing,nothing,nothing,Dict{String,DriverRT.ContainerScientificHandle}(),
            false,false,false,ReentrantLock(),ReentrantLock())
        profile = only(values(config.profiles.definitions))
        fence = DriverRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main",
            config.worker_id,string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
        try
            action((;driver,profile,fence,directory))
        finally
            # This fixture opens no physical resources. Its synthetic receipts
            # are removed explicitly in each test before closing the driver.
            close(driver)
        end
    end
end

@testset "container scientific adapter preserves admission and partial ownership" begin
    @test ManagedScientificDriver === ContainerScientificDriver
    @test DriverRT.RequiredInterfaces.check_interface_implemented(AbstractScientificDriver,ContainerScientificDriver)
    with_container_driver_fixture() do f
        (;driver,profile,fence) = f
        @test isempty(driver.handles) && isempty(resource_receipts(driver.journal))
        @test installed_profiles(driver) === driver.config.profiles
        @test !occursin(f.directory,repr(MIME"text/plain"(),driver))
        @test_throws ArgumentError DriverRT.container_driver_profile(driver,profile,fence)
        # Internal fixture only; public startup requires actual manager proof.
        driver.recovered = true
        @test DriverRT.container_driver_profile(driver,profile,fence) === nothing
        other = DriverRT.Protocol.AssignmentFence((field == :worker_id ? "foreign" : getfield(fence,field)
            for field in fieldnames(typeof(fence)))...)
        @test_throws ArgumentError verify_executor!(driver,profile,other)
        @test isempty(driver.handles)
        @test_throws CommandFailure executor_for!(driver,profile,fence)
        handle = only(values(driver.handles))
        @test handle.fence == fence && handle.receipt === nothing && handle.supervisor === nothing
        @test release_owned!(driver,fence)
        @test isempty(driver.handles)
        @test release_owned!(driver,fence)
        @test close(driver) === nothing
        @test driver.closed && driver.cleanup_complete
        @test close(driver) === nothing
        @test_throws ArgumentError executor_for!(driver,profile,fence)
    end
end

@testset "shared scientific owner rejects native admission without native host proof" begin
    with_container_driver_fixture() do f
        (;driver,fence) = f
        native = ProfileDefinition("fixture","/fixture/project",repeat("a",64);operations=["fixture.echo"])
        driver.profiles.definitions[native.id] = native
        driver.recovered = true
        @test_throws CommandFailure verify_executor!(driver,native,fence)
        @test_throws CommandFailure executor_for!(driver,native,fence)
        @test isempty(resource_receipts(driver.journal))
        handle = only(values(driver.handles))
        @test handle.supervisor === nothing && handle.receipt === nothing
        @test release_owned!(driver,fence)
        @test isempty(driver.handles) && isempty(driver.runner.active)
    end
end

# Exercise the same internal native acquisition method without substituting
# production authority hooks or starting a command. Only host proof is synthetic;
# the installed scientific source fingerprint and journal operations are real.
struct NativeAcquisitionFixture
    host::NativeHostCheck
    manager::ManagedAgentIdentity
    journal::ResourceJournal
end
DriverRT.verified_native_driver_host(f::NativeAcquisitionFixture) = f.host

@testset "native command acquisition retains intent before any process starts" begin
    mktempdir() do directory
        project = normpath(joinpath(@__DIR__,"..","..","worker","profiles","line-parameters"))
        fingerprint = native_environment_fingerprint(project)
        profile = ProfileDefinition("fixture",project,fingerprint.digest;operations=["fixture.echo"])
        journal = ResourceJournal(joinpath(directory,"journal"),"worker-a";capacity=1)
        manager = ManagedAgentIdentity("lcm-agent-worker-a.service",repeat("b",32),"/private/agent",
            "/private/config",journal.root,"worker-a")
        fence = DriverRT.Protocol.AssignmentFence(string(uuid4()),string(uuid4()),"alice","main","worker-a",
            string(uuid4()),string(uuid4()),profile.id,string(profile.version),profile.fingerprint,1)
        handle = DriverRT.ManagedScientificHandle(fence,nothing,nothing,false,ReentrantLock())
        host = NativeHostCheck(repeat("a",64),"/usr/bin/systemd-run",1000,1000,())
        fixture = NativeAcquisitionFixture(host,manager,journal)
        try
            @test DriverRT.acquire_scientific_backend!(Val(:trusted_process),fixture,handle,profile,fence) === nothing
            @test handle.receipt.backend == :native && handle.receipt.physical_id === nothing
            @test only(resource_receipts(journal)).id == handle.receipt.id
            @test handle.supervisor.process === nothing && handle.supervisor.generation == 0
            @test handle.supervisor.discard_stderr
            @test "--unit=" * resource_name(handle.receipt) in collect(handle.supervisor.command)
            @test !any(m->nameof(m) in (:LineCableModels,:PowerImpedance),values(Base.loaded_modules))
        finally
            # No process/service was started; retire only this fixture's intent.
            handle.supervisor === nothing || DriverRT.ExecutionCore.stop_executor!(handle.supervisor)
            for receipt in resource_receipts(journal); forget_resource!(journal,receipt); end
            close(journal)
        end
    end
end

@testset "unresolved physical receipts retain driver ownership through close" begin
    with_container_driver_fixture() do f
        (;driver,profile,fence) = f
        driver.recovered = true
        # A synthetic native intent has a foreign scope and cannot be verified;
        # it must not be silently forgotten or translated into a container ID.
        receipt = reserve_resource!(driver.journal,fence,:native,repeat("b",64))
        supervisor = DriverRT.ScientificProcess(;command=`false`)
        handle = DriverRT.ContainerScientificHandle(fence,nothing,supervisor,false,ReentrantLock())
        driver.handles[fence.lease_id] = handle
        try
            @test !release_owned!(driver,fence)
            @test handle.closing && haskey(driver.handles,fence.lease_id)
            @test only(resource_receipts(driver.journal)).id == receipt.id
            @test_throws ArgumentError close(driver)
            @test driver.closed && !driver.cleanup_complete
            @test_throws ArgumentError ResourceJournal(driver.journal.root,driver.config.worker_id)
        finally
            # No physical resource was created by this synthetic receipt.
            forget_resource!(driver.journal,receipt)
        end
        @test close(driver) === nothing
        @test driver.cleanup_complete && isempty(driver.handles)
    end
end
