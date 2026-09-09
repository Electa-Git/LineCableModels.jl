# Opt-in consumer gate against separately provisioned physical agents. The
# caller owns their installation/cleanup; this process never launches an engine.
using Test, UUIDs, Sockets, LineCableModelsRuntime
import JSON3, HTTP
const RT = LineCableModelsRuntime
include(joinpath(@__DIR__,"..","..","src","scientific","StudyCases.jl"))
include("terminal_live_browser.jl")

function physical_consumers(config_file,directory,stop_command_file)
    mode=get(ENV,"LCM_PHYSICAL_CONSUMERS","full")
    mode in ("full","terminal") || error("Physical consumer mode must be full or terminal")
    directory=abspath(directory)
    !ispath(directory) && !islink(directory) || error("Use a new diagnostic directory; do not replay old fault markers")
    config = read_control_config(config_file)
    Set(keys(config.workers))==Set(["worker-a","worker-b"]) || error("Expected the two owned acceptance agents")
    all(p->p.isolation==:container,values(config.profiles.definitions)) || error("Physical container profiles required")
    config.artifacts isa S3RuntimeArtifacts || error("Independent private S3 artifact transport required")
    # Operator-owned, private argv file; never an HTTP/browser-supplied command.
    RT.broker_file(stop_command_file;private=true,max_bytes=8192)
    stop_args = JSON3.read(read(stop_command_file,String),Vector{String})
    !isempty(stop_args) && isabspath(first(stop_args)) || error("Explicit owned stop command required")
    mkpath(directory;mode=0o700)
    applications = default_applications()
    owner,operator = Principal("scientific-browser"),Principal("operator";administrator=true)
    store = RuntimeStore(joinpath(directory,"runtime.sqlite"))
    control = ControlService(config,store,applications)
    supervisor = UIHostSupervisor(store,applications,joinpath(directory,"hosts");
        limits=RunLimits(max_runs=2,max_runs_per_owner=2,startup_seconds=120,shutdown_seconds=2))
    listener=Sockets.listen(ip"127.0.0.1",0)
    port=Int(last(Sockets.getsockname(listener)));close(listener)
    origin="http://127.0.0.1:$port"
    site=PublishedSite(normpath(joinpath(@__DIR__,"..","..","_site")))
    server=start_gateway(supervisor,LocalIdentity(origin,owner);port,site,control)
    browser=nothing; fault=nothing
    try
        for trust in values(config.workers)
            enroll_worker!(store,operator,trust)
            set_registration_state!(store,operator,trust.worker_id,:approved;expected_revision=1)
        end
        start_control!(control)
        @test timedwait(180;pollint=0.1) do
            lock(control.inventory.lock) do
                all(RT.presence_state(control.inventory,get(control.inventory.presence,id,nothing))==:online
                    for id in (mode=="terminal" ? ("worker-a",) : keys(config.workers)))
            end && control.jobs.state==:online
        end == :ok
        @test all(isempty(list_jobs(store,owner,run.id)) for run in list_runs(store,owner))
        write(joinpath(directory,"gateway.json"),JSON3.write((;origin)))
        if mode=="full"
            driver=joinpath(@__DIR__,"physical_fault_driver.mjs")
            fault=run(pipeline(`node $driver $directory $stop_command_file`;stdout=stdout,stderr=stderr);wait=false)
            deck=start_ui!(supervisor,owner,"ichqp-showcase")
            start_ui!(supervisor,owner,"cable-study")
            script=joinpath(@__DIR__,"scientific_live_browser.mjs")
            command=addenv(`node $script $origin $directory`,"LCM_TEST_KEEP_CONSUMERS"=>"1")
            browser=run(pipeline(command;stdout=stdout,stderr=stderr);wait=false)
            ready=joinpath(directory,"scientific-browser-ready")
            @test timedwait(()->isfile(ready)||process_exited(browser),2450;pollint=0.2)==:ok
            isfile(ready) || error("Physical consumer browser failed; retained diagnostics identify the gate")
            @test timedwait(()->process_exited(fault),125)==:ok
            process_exited(fault) || error("Owned remote stop did not finish")
            wait(fault);@test success(fault)
            success(fault) || error("Owned remote stop failed")
            @test isfile(joinpath(directory,"power-worker-stopped"))
        else
            deck=reserve_run!(store,owner,applications.definitions["ichqp-showcase"])
        end
        # Exercise the unchanged Chrome/xterm renderer against the remote
        # container after the scientific worker-loss and artifact checks.
        terminal_live_browser(control,supervisor,owner,deck,directory) do
            lease=reserve_assignment!(control.coordinator.assignments,owner,deck.id,
                "terminal","julia-terminal";placement=PinnedPlacement("worker-a"))
            id=UUID(lease.fence.lease_id)
            grant_assignment!(control.coordinator,owner,id)
            @test timedwait(()->assignment_usable(control.coordinator,owner,id),5)==:ok
        end
        if browser!==nothing
            write(joinpath(directory,"terminal-browser-finished"),"finished\n")
            @test timedwait(()->process_exited(browser),15)==:ok
            process_exited(browser) || error("Scientific browser did not close")
            wait(browser);@test success(browser)
        end
    finally
        if browser!==nothing && !process_exited(browser)
            kill(browser,Base.SIGTERM)
            timedwait(()->process_exited(browser),10)==:ok || kill(browser,Base.SIGKILL)
            wait(browser)
        end
        if fault!==nothing && !process_exited(fault)
            kill(fault,Base.SIGTERM)
            timedwait(()->process_exited(fault),10)==:ok || kill(fault,Base.SIGKILL)
            wait(fault)
        end
        close(server);close(supervisor);close(control);close(store)
    end
    @test isempty(supervisor.handles)
    # Numerical comparison belongs after live ownership is closed; compiling
    # test-only projections must not consume a live lease's ACK budget.
    mode=="full" && verify_physical_consumer_results(directory)
    @test !any(id.name in ("Bonito","LineCableModelsPlayground","LineCableModels","PowerImpedance") for id in keys(Base.loaded_modules))
end

function verify_physical_consumer_results(directory)
    data=JSON3.read(read(joinpath(directory,"scientific-results.json"),String),Dict{String,Any})
    for (name,case) in (("parameters",StudyCases.LineParameters()),("corridor",StudyCases.CorridorImpedance()))
        a,b=data["deck"][name],data["workbench"][name]
        expected=StudyCases.inputs(case;minimum_frequency_hz=100,maximum_frequency_hz=150,frequency_points=2)
        @test a["parameters"]==b["parameters"]==expected
        @test a["receipt"]["run_id"]!=b["receipt"]["run_id"]
        @test a["receipt"]["executor_id"]!=b["receipt"]["executor_id"]
        @test a["receipt"]["input_hash"]==b["receipt"]["input_hash"]==RT.Protocol.input_hash(StudyCases.operation(case),expected)
        left,right=StudyCases.result_series(case,a["value"]),StudyCases.result_series(case,b["value"])
        @test left.frequency≈[100.0,150.0] && left.frequency==right.frequency
        @test all(x.values≈y.values for (x,y) in zip(left.curves,right.curves))
    end
end

@testset "registered consumers across a physical worker boundary" begin
    length(ARGS)==3 || error("Expected private control configuration, new diagnostic directory, and owned stop argv file")
    physical_consumers(ARGS...)
end
