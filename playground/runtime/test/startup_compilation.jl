@testset "startup compilation does not acquire authority or execute hooks" begin
    mktempdir() do directory
        (;driver,resources,agent,fences) = scientific_fixture(directory;capacity=1)
        store = RuntimeStore(joinpath(directory,"compilation.sqlite"))
        trust = WorkerTrust("worker-a","credential-worker-a",("fixture",))
        control = ControlService(ControlConfig(agent.config.endpoint,driver.profiles,[trust]),
            store,ApplicationRegistry())
        try
            @test RT.compile_runtime_paths(agent) === nothing
            @test driver.verified == driver.recovered == 0
            @test isempty(driver.processes) && isempty(resources.handles)
            @test length(agent.ledger.leases) == 1 && agent_lease_usable(agent.ledger,only(fences))
            @test agent.task === agent.connector === agent.science.connection === agent.jobs.connection === nothing
            @test isempty(agent.jobs.flights) && agent.jobs.task===agent.jobs.connector===nothing
            @test RT.compile_runtime_paths(control) === nothing
            @test control.task === control.connector === control.link.control === nothing
            @test control.jobs.connection === control.science.connection === nothing
            @test isempty(control.jobs.flights) && isempty(control.jobs.retry_at)
            @test isempty(control.coordinator.flights) && isempty(control.events.records)
            @test isempty(list_assignments(store,Principal("operator";administrator=true)))
            @test isempty(list_runs(store,Principal("operator";administrator=true)))
        finally
            close(agent); close(control); close(store)
        end
    end
end

@testset "terminal compilation remains passive before announcing availability" begin
    with_terminal_sessions(;capacity=1) do f
        verified,recovered = f.driver.verified,f.driver.recovered
        @test RT.compile_runtime_paths(f.agent) === nothing
        @test f.driver.verified == verified && f.driver.recovered == recovered
        @test isempty(f.driver.handles) && isempty(f.resources.handles)
        @test f.agent.task === f.agent.connector === f.agent.terminals.connection === nothing
    end
end
