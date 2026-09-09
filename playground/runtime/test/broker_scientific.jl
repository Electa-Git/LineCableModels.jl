# Real registered consumers and numerical adapters over protected HTTP/TLS.
# The finite native driver is test-only; this does not certify OS isolation.
using Test, UUIDs, Dates, Sockets, LineCableModelsRuntime
import JSON3
const RT = LineCableModelsRuntime
const P = RT.Protocol
directory, port, artifact_port = ARGS
certs = joinpath(directory, "certs")
endpoint(id) = BrokerEndpoint("tls://127.0.0.1:$port", joinpath(directory, id * ".password");
    ca_file=joinpath(certs,"ca.pem"), certificate_file=joinpath(certs,"worker-cert.pem"),
    key_file=joinpath(certs,"worker-key.pem"), server_name="localhost")
artifact_location(id) = S3RuntimeArtifacts("https://127.0.0.1:$artifact_port", "lcm-runtime-private", "runtime-v1",
    joinpath(directory,"artifact-$id.toml"); ca_file=joinpath(certs,"ca.pem"))
include("scientific_driver_fixture.jl")
include(joinpath(@__DIR__,"..","..","src","scientific","StudyCases.jl"))

@testset "registered scientific consumers execute through real protected transport" begin
    profiles = ProfileRegistry()
    commands = Dict{String,Cmd}()
    for (id, operations) in (("line-parameters",("line.frequency_scan",)),
            ("power-flow",("powerflow.prepare","impedance.evaluate")))
        project = normpath(joinpath(@__DIR__,"..","..","worker","profiles",id))
        register!(profiles, ProfileDefinition(id,project,native_environment_fingerprint(project).digest;
            operations, budget=ResourceBudget(memory_bytes=4*1024^3,prepare_seconds=600,job_seconds=120)))
        commands[id] = fixture_command(project,joinpath(@__DIR__,"scientific_real_child.jl"),id)
    end
    store = RuntimeStore(joinpath(directory,"runtime.sqlite"))
    applications = default_applications()
    owner, operator = Principal("scientific-browser"), Principal("operator";administrator=true)
    trusts = [WorkerTrust("worker-a","credential-worker-a",("line-parameters",);capacity=2),
        WorkerTrust("worker-b","credential-worker-b",("power-flow",);capacity=2)]
    agents, drivers = AgentService[], ScientificDriverFixture[]
    for trust in trusts
        enroll_worker!(store,operator,trust)
        set_registration_state!(store,operator,trust.worker_id,:approved;expected_revision=1)
        installed = ProfileRegistry()
        foreach(id->register!(installed,profiles.definitions[id]),trust.profiles)
        driver = ScientificDriverFixture(installed,commands,
            Dict{String,Tuple{P.AssignmentFence,RT.ExecutionCore.ExecutorSupervisor}}(),0,0,false,true,false)
        push!(drivers,driver)
        push!(agents,AgentService(AgentConfig(trust.worker_id,endpoint(trust.worker_id),installed,
            joinpath(directory,trust.worker_id); capacity=2,artifacts=artifact_location(trust.worker_id)),ScientificResources(driver)))
    end
    control = ControlService(ControlConfig(endpoint("coordinator"),profiles,trusts;
        artifacts=artifact_location("coordinator")),store,applications)
    supervisor = UIHostSupervisor(store,applications,joinpath(directory,"hosts");
        limits=RunLimits(max_runs=2,max_runs_per_owner=2,startup_seconds=120,shutdown_seconds=1))
    listener = Sockets.listen(ip"127.0.0.1",0)
    webport = Int(last(Sockets.getsockname(listener))); close(listener)
    origin = "http://127.0.0.1:$webport"
    site = PublishedSite(normpath(joinpath(@__DIR__,"..","..","_site")))
    server = start_gateway(supervisor,LocalIdentity(origin,owner);port=webport,site,control)
    browser = nothing
    running = Ref(true)
    fault = nothing
    try
        foreach(start_agent!,agents); start_control!(control)
        @test timedwait(()->all(haskey(control.inventory.presence,w.worker_id) &&
            control.inventory.presence[w.worker_id].report.sequence>=2 for w in trusts) &&
            control.jobs.state==:online && all(a->a.jobs.state==:online,agents),30)==:ok
        @test all(isempty, (driver.processes for driver in drivers))
        for trust in trusts, stream in (RT.job_stream(trust.worker_id),RT.result_stream(trust.worker_id))
            info = RT.job_request(control.jobs.connection,RT.NATS.JetStream.StreamInfo,"\$JS.API.STREAM.INFO.$stream")
            @test info.config.name == stream
        end
        start_ui!(supervisor,owner,"ichqp-showcase")
        start_ui!(supervisor,owner,"cable-study")
        # A local test signal may stop only this fixture's power-flow agent.
        # No diagnostic/fault endpoint is exposed through the application.
        fault = @async while running[]
            if isfile(joinpath(directory,"stop-power-worker"))
                close(agents[2])
                write(joinpath(directory,"power-worker-stopped"),"stopped\n")
                break
            end
            sleep(0.1)
        end
        script = joinpath(@__DIR__,"scientific_live_browser.mjs")
        browser = Base.run(pipeline(`node $script $origin $directory`;stdout=stdout,stderr=stderr);wait=false)
        # Four independent cold profiles, one post-cancellation replacement,
        # and bounded browser checks. Individual runtime deadlines are unchanged.
        @test timedwait(()->process_exited(browser),2450;pollint=0.2)==:ok
        if !process_exited(browser)
            kill(browser,Base.SIGTERM)
            timedwait(()->process_exited(browser),10)==:ok || kill(browser,Base.SIGKILL)
        end
        wait(browser)
        @test success(browser)
        if !success(browser)
            # Operator-side fixture diagnostics only; never expose raw failures
            # to browser clients or change the production scheduler's behavior.
            for trust in trusts
                try
                    RT.ensure_worker_streams!(control.jobs.connection,trust)
                    println("Post-failure stream check: ",trust.worker_id," compatible")
                catch error
                    showerror(stderr,error,catch_backtrace());println(stderr)
                end
            end
            for run in list_runs(store,owner), receipt in list_jobs(store,owner,run.id)
                println("Post-failure receipt: ",receipt.state," · ",receipt.job.request.job_id)
                try
                    owned_job_result(control.jobs,owner,UUID(receipt.job.request.job_id))
                catch error
                    showerror(stderr,error,catch_backtrace());println(stderr)
                end
            end
            error("Scientific browser failed; retained diagnostics identify the failed gate")
        end
        evidence = JSON3.read(read(joinpath(directory,"scientific-results.json"),String),Dict{String,Any})
        for (name,case) in (("parameters",StudyCases.LineParameters()),("corridor",StudyCases.CorridorImpedance()))
            left,right = evidence["deck"][name],evidence["workbench"][name]
            expected = StudyCases.inputs(case;minimum_frequency_hz=100,maximum_frequency_hz=150,frequency_points=2)
            @test left["parameters"] == right["parameters"] == expected
            @test left["receipt"]["input_hash"] == right["receipt"]["input_hash"] ==
                P.input_hash(StudyCases.operation(case),expected)
            @test left["receipt"]["run_id"] != right["receipt"]["run_id"]
            @test left["receipt"]["executor_id"] != right["receipt"]["executor_id"]
            a,b = StudyCases.result_series(case,left["value"]),StudyCases.result_series(case,right["value"])
            @test a.frequency == b.frequency
            # Existing logarithmic sampling may return 150 + one Float64 ULP.
            @test a.frequency ≈ [100.0,150.0]
            @test all(x.values ≈ y.values for (x,y) in zip(a.curves,b.curves))
            if name == "parameters"
                for quantity in ("reactance","conductance","susceptance")
                    a,b = StudyCases.result_series(case,left["value"],quantity),StudyCases.result_series(case,right["value"],quantity)
                    @test all(x.values ≈ y.values for (x,y) in zip(a.curves,b.curves))
                end
            else
                for entry in (left,right)
                    curves = entry["value"]["curves"]
                    @test [curves[key]["corridor_length_m"] for key in ("base_minus_error","base","base_plus_error")] ==
                        [95000.0,100000.0,105000.0]
                    @test curves["base_minus_error"]["magnitude_db_ohm"] != curves["base_plus_error"]["magnitude_db_ohm"]
                end
            end
        end
        @test !any(id.name in ("Bonito","LineCableModelsPlayground","LineCableModels","PowerImpedance")
            for id in keys(Base.loaded_modules))
    finally
        running[] = false
        if browser !== nothing && !process_exited(browser)
            kill(browser,Base.SIGTERM)
            timedwait(()->process_exited(browser),10)==:ok || kill(browser,Base.SIGKILL)
            wait(browser)
        end
        fault === nothing || wait(fault)
        close(server)
        close(supervisor)
        foreach(close,agents)
        close(control)
        close(store)
    end
    @test all(driver->isempty(driver.processes),drivers)
    @test isempty(supervisor.handles)
end
