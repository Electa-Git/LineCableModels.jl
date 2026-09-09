include("agent_resource_fixture.jl")

@testset "agent configuration and required resource ownership remain passive" begin
    mktempdir() do dir
        data = control_config_data()
        delete!(data, "workers")
        data["agent"] = Dict{String,Any}("worker_id"=>"worker-a", "scratch_root"=>"owned-agent",
            "container_runtime"=>"podman")
        secret = joinpath(dir, "coordinator.password")
        write(secret, "agent-config-fixture-secret"); chmod(secret, 0o600)
        path = joinpath(dir, "agent.toml")
        function parse_agent(value)
            open(path, "w") do io
                TOML.print(io, value)
            end
            read_agent_config(path)
        end
        config = parse_agent(data)
        @test config.worker_id == "worker-a"
        @test config.container_runtime == :podman
        @test config.scratch_root == joinpath(dir, "owned-agent")
        @test !ispath(config.scratch_root)
        @test config.artifacts===nothing
        stored=deepcopy(data)
        stored["artifacts"]=Dict("backend"=>"filesystem","root"=>"private-results")
        @test parse_agent(stored).artifacts.root==joinpath(dir,"private-results")
        @test !ispath(joinpath(dir,"private-results"))
        parse_agent(data)
        output = joinpath(dir, "check-agent.txt")
        open(output, "w") do io
            redirect_stdout(io) do
                runtime_cli(["runtime", "check-agent", "--config", path])
            end
        end
        @test occursin("No profile was prepared", read(output, String))
        @test !ispath(config.scratch_root)
        @test !occursin(dir, repr(config))
        @test !occursin("fixture-secret", repr(config))
        resources = AgentResourceFixture(config.profiles)
        @test_throws ArgumentError AgentService(config, IncompleteAgentResources())
        agent = AgentService(config, resources)
        @test resources.recovered == 0 && agent.task === nothing && agent.connector === nothing
        @test isempty(agent.ledger.leases) && agent.state == :offline
        other = ProfileRegistry()
        register!(other, ProfileDefinition("foreign", "/unapproved", repeat("b",64); operations=("system.echo",)))
        @test_throws ArgumentError AgentService(config, AgentResourceFixture(other))
        disabled = AgentService(config, AgentResourceFixture(ProfileRegistry()))
        @test isempty(disabled.ledger.profiles.definitions)
        close(disabled)
        resources.recovery_fails = true
        @test_throws ArgumentError start_agent!(agent)
        @test agent.task === nothing && agent.link.control === nothing
        close(agent)
        @test resources.closed
        @test_throws ArgumentError start_agent!(agent)
        for runtime in ("auto", "podman", "docker")
            data["agent"]["container_runtime"] = runtime
            @test parse_agent(data).container_runtime == Symbol(runtime)
        end
        data["agent"]["container_runtime"] = "browser-command"
        @test_throws ArgumentError parse_agent(data)
        data["agent"]["container_runtime"] = "auto"
        data["agent"]["scratch_root"] = "."
        @test_throws ArgumentError parse_agent(data)
        data["agent"]["scratch_root"] = "owned-agent"
        data["agent"]["capacity"] = true
        @test_throws ArgumentError parse_agent(data)
    end
end

@testset "failed scheduler cannot bypass owned resource retirement" begin
    mktempdir() do dir
        profiles = ProfileRegistry()
        register!(profiles, ProfileDefinition("fixture", "/approved", repeat("a",64); operations=("system.echo",)))
        config = AgentConfig("worker-a", BrokerEndpoint("tls://broker.invalid", "/private/password"),
            profiles, joinpath(dir, "owned"))
        for failure in (InterruptException(), ErrorException("scheduler fixture failure"))
            resources = AgentResourceFixture(profiles)
            agent = AgentService(config, resources)
            agent.task = @async throw(failure)
            wait(agent.task; throw=false)
            @test istaskfailed(agent.task)
            @test close(agent) === nothing
            @test resources.closed && agent.cleanup_complete && agent.state == :stopped
            @test close(agent) === nothing
        end
    end
end

@testset "agent cleanup never blocks control or claims an unresolved resource is released" begin
    mktempdir() do dir
        profiles = ProfileRegistry()
        register!(profiles, ProfileDefinition("line-parameters", "/approved", repeat("a",64); operations=("system.echo",)))
        config = AgentConfig("worker-a", BrokerEndpoint("tls://broker.invalid", "/private/password"),
            profiles, joinpath(dir, "owned"))
        resources = AgentResourceFixture(profiles)
        clock = Ref(0.0)
        agent = AgentService(config, resources; clock=()->clock[])
        P = RT.Protocol
        coordinator_id = string(uuid4())
        probe = P.WorkerProbe("2.0", "worker-a", coordinator_id, string(uuid4()))
        @test receive_probe!(agent.ledger, probe)
        fence = P.AssignmentFence(string(uuid4()), string(uuid4()), "alice", "main",
            "worker-a", agent.ledger.boot_id, coordinator_id, "line-parameters", "1.0.0", repeat("a",64), 1)
        grant = P.LeaseControl("2.0", string(uuid4()), "grant", fence, 1, 10_000)
        @test handle_lease_control!(agent.ledger, grant).accepted
        release = P.LeaseControl("2.0", string(uuid4()), "release", fence, 2, 0)
        resources.release_allowed[] = false
        @test handle_lease_control!(agent.ledger, release) === nothing
        tick_agent!(agent) # first specialization; cleanup becomes independently pending
        yield()
        @test length(agent.cleanups) == 1
        @test isempty(resources.released) && !agent_lease_usable(agent.ledger, fence)
        for _ in 1:5
            @test (@elapsed tick_agent!(agent)) < 0.1
            @test receive_probe!(agent.ledger, P.WorkerProbe("2.0", "worker-a", coordinator_id, string(uuid4())))
        end
        @test only(values(agent.ledger.leases)).state == :closing
        @test length(agent.cleanups) == 1
        resources.release_allowed[] = true
        @test timedwait(() -> istaskdone(only(values(agent.cleanups))), 2) == :ok
        tick_agent!(agent)
        @test only(values(agent.ledger.leases)).state == :closed
        @test isempty(agent.cleanups) && isempty(agent.cleanup_retry)
        @test handle_lease_control!(agent.ledger, release).accepted
        @test !handle_lease_control!(agent.ledger, grant).accepted
        next_fence = P.AssignmentFence(string(uuid4()), fence.run_id, fence.owner, fence.role,
            fence.worker_id, fence.worker_boot, fence.coordinator_id, fence.profile_id,
            fence.profile_version, fence.fingerprint, 2)
        @test handle_lease_control!(agent.ledger,
            P.LeaseControl("2.0", string(uuid4()), "grant", next_fence, 1, 10_000)).accepted
        resources.release_allowed[] = false
        @test handle_lease_control!(agent.ledger,
            P.LeaseControl("2.0", string(uuid4()), "release", next_fence, 2, 0)) === nothing
        tick_agent!(agent)
        yield()
        first_close = @async close(agent)
        second_close = @async close(agent)
        wait(second_close)
        @test istaskdone(first_close) && agent.cleanup_complete
        @test agent.state == :stopped && resources.closed

        failed_resource = AgentResourceFixture(profiles)
        failed_resource.release_fails[] = true
        failing = AgentService(config, failed_resource; clock=()->clock[])
        @test receive_probe!(failing.ledger, probe)
        failed_fence = P.AssignmentFence(string(uuid4()), string(uuid4()), "alice", "main",
            "worker-a", failing.ledger.boot_id, coordinator_id, "line-parameters", "1.0.0", repeat("a",64), 1)
        @test handle_lease_control!(failing.ledger,
            P.LeaseControl("2.0", string(uuid4()), "grant", failed_fence, 1, 10_000)).accepted
        @test handle_lease_control!(failing.ledger,
            P.LeaseControl("2.0", string(uuid4()), "release", failed_fence, 2, 0)) === nothing
        tick_agent!(failing); yield(); tick_agent!(failing)
        @test only(values(failing.ledger.leases)).state == :closing
        @test isempty(failing.cleanups) && length(failing.cleanup_retry) == 1
        clock[] = 0.5
        tick_agent!(failing)
        @test isempty(failing.cleanups) # retry rate is bounded, not one task per tick
        @test_throws ArgumentError close(failing)
        @test !failing.cleanup_complete
        @test_throws ArgumentError close(failing)
    end
end
