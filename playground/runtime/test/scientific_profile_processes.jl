using Test, Dates, UUIDs, LineCableModelsRuntime
const RT = LineCableModelsRuntime
const P = RT.Protocol
lease_command(fence; action="grant", revision=1, duration_ms=60000) =
    P.LeaseControl("2.0", string(uuid4()), action, fence, revision, duration_ms)
include("scientific_driver_fixture.jl")

@testset "real scientific profiles execute under independent live leases" begin
    mktempdir() do directory
        profiles = ProfileRegistry()
        commands = Dict{String,Cmd}()
        for (id, operations) in (("line-parameters", ("line.frequency_scan",)),
                ("power-flow", ("powerflow.prepare", "impedance.evaluate")))
            project = normpath(joinpath(@__DIR__, "..", "..", "worker", "profiles", id))
            fingerprint = native_environment_fingerprint(project)
            register!(profiles, ProfileDefinition(id, project, fingerprint.digest; operations,
                budget=ResourceBudget(memory_bytes=4*1024^3, prepare_seconds=600, job_seconds=120)))
            commands[id] = fixture_command(project, joinpath(@__DIR__, "scientific_real_child.jl"), id)
        end
        driver = ScientificDriverFixture(profiles, commands,
            Dict{String,Tuple{P.AssignmentFence,RT.ExecutionCore.ExecutorSupervisor}}(), 0, 0, false, true, false)
        resources = ScientificResources(driver)
        config = AgentConfig("worker-a", BrokerEndpoint("tls://broker.invalid", "/unused-password"),
            profiles, joinpath(directory, "owned"); capacity=2)
        agent = AgentService(config, resources)
        recover_owned!(resources)
        coordinator = string(uuid4())
        probe() = P.WorkerProbe("2.0", "worker-a", coordinator, string(uuid4()))
        receive_probe!(agent.ledger, probe())
        run_id = string(uuid4())
        fences = Dict(id => P.AssignmentFence(string(uuid4()), run_id, "scientific-test", id,
            "worker-a", agent.ledger.boot_id, coordinator, id, "1.0.0", profile.fingerprint, 1)
            for (id, profile) in profiles.definitions)
        for fence in values(fences)
            @test handle_lease_control!(agent.ledger, lease_command(fence)).accepted
        end
        maintaining = Ref(true)
        probes = Ref(0)
        maintenance = @async while maintaining[]
            receive_probe!(agent.ledger, probe())
            probes[] += 1
            for lease in values(agent.ledger.leases)
                if lease.state == :active && lease.expires_at - agent.ledger.clock() < 30
                    @assert handle_lease_control!(agent.ledger, lease_command(lease.fence;
                        action="renew", revision=lease.command.revision+1)).accepted
                end
            end
            sleep(1)
        end
        inputs = Dict{String,Any}("frequencies_hz"=>[50.0], "separation_m"=>0.5, "depth_m"=>1.0,
            "earth_resistivity_ohm_m"=>100.0, "line_length_m"=>1000.0)
        flow_inputs = Dict("specification"=>Dict("earth_resistivity_ohm_m"=>100.0))
        try
            @test all(scientific_status(resources, fence).preparation == :cold for fence in values(fences))
            line = fetch(prepare_assigned!(resources, fences["line-parameters"], inputs))
            @test line["cache_status"] == "miss"
            power_task = prepare_assigned!(resources, fences["power-flow"], flow_inputs)
            independent = @elapsed result = fetch(execute_assigned!(resources,
                assigned_fixture_job(resources,fences["line-parameters"], "line.frequency_scan", inputs;timeout=Second(120)))).value
            @test result["frequencies_hz"] == [50.0]
            @test haskey(result, "series_impedance_ohm_per_m")
            @test independent < 30
            @test !istaskdone(power_task)
            power = fetch(power_task)
            @test power["evidence"]["preparation_kind"] == "solved_and_linearized_model"
            @test scientific_status(resources, fences["power-flow"]).preparation == :ready
            @test probes[] > 5 && !istaskfailed(maintenance)
            @test all(agent_lease_usable(agent.ledger, fence) for fence in values(fences))
            @test all(!haskey(driver.processes[id][2].command.env === nothing ? Dict() :
                Dict(split(entry,'=';limit=2)[1]=>true for entry in driver.processes[id][2].command.env), "NATS_CONNECT_URL")
                for id in keys(driver.processes))
            previous = scientific_status(resources, fences["power-flow"])
            warm = @elapsed repeated = fetch(prepare_assigned!(resources, fences["power-flow"], flow_inputs))
            @test repeated["cache_status"] == "hit"
            @test scientific_status(resources, fences["power-flow"]).executor_id == previous.executor_id
            release = lease_command(fences["power-flow"]; action="release",
                revision=agent.ledger.leases[(run_id,"power-flow")].command.revision+1, duration_ms=0)
            @test handle_lease_control!(agent.ledger, release) === nothing
            @test release_owned!(resources, fences["power-flow"])
            @test complete_agent_cleanup!(agent.ledger, fences["power-flow"]).accepted
            @test fetch(execute_assigned!(resources, assigned_fixture_job(resources,fences["line-parameters"],
                "line.frequency_scan", inputs;timeout=Second(120)))).value["frequencies_hz"] == [50.0]
            @test !any(id.name in ("LineCableModels", "PowerImpedance", "Bonito") for id in keys(Base.loaded_modules))
            println("Lease-owned independent line request seconds=$independent; cached power preparation seconds=$warm")
        finally
            maintaining[] = false
            wait(maintenance)
            close(agent)
        end
        @test isempty(driver.processes) && isempty(resources.handles)
    end
end
