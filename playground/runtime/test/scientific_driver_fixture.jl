# Deliberately test-only native driver: no container/OS-budget certification.
mutable struct ScientificDriverFixture <: AbstractScientificDriver
    profiles::ProfileRegistry
    commands::Dict{String,Cmd}
    processes::Dict{String,Tuple{RT.Protocol.AssignmentFence,RT.ExecutionCore.ExecutorSupervisor}}
    recovered::Int
    verified::Int
    fail_verify::Bool
    allow_release::Bool
    closed::Bool
end
RT.installed_profiles(driver::ScientificDriverFixture) = driver.profiles
RT.recover_owned!(driver::ScientificDriverFixture) = (driver.recovered += 1; nothing)
function RT.verify_executor!(driver::ScientificDriverFixture, profile::ProfileDefinition, _::RT.Protocol.AssignmentFence)
    driver.verified += 1
    driver.fail_verify && throw(ArgumentError("private source path must not be exposed"))
    verify_native_environment(profile)
    return nothing
end
function RT.executor_for!(driver::ScientificDriverFixture, profile::ProfileDefinition, fence::RT.Protocol.AssignmentFence)
    if haskey(driver.processes, fence.lease_id)
        previous, supervisor = driver.processes[fence.lease_id]
        previous == fence || error("fixture fence changed")
        return supervisor
    end
    supervisor = RT.ExecutionCore.ExecutorSupervisor(profile.environment;
        command=driver.commands[profile.id], startup_timeout_seconds=120)
    driver.processes[fence.lease_id] = (fence, supervisor)
    return supervisor
end
function RT.release_owned!(driver::ScientificDriverFixture, fence::RT.Protocol.AssignmentFence)
    driver.allow_release || return false
    previous = get(driver.processes, fence.lease_id, nothing)
    if previous !== nothing
        previous[1] == fence || error("fixture cleanup fence changed")
        RT.ExecutionCore.stop_executor!(previous[2])
        delete!(driver.processes, fence.lease_id)
    end
    return true
end
function Base.close(driver::ScientificDriverFixture)
    for (fence, _) in collect(values(driver.processes))
        RT.release_owned!(driver, fence) || error("fixture release refused")
    end
    driver.closed = true
    return nothing
end

function fixture_command(project, child, args...)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    setenv(`setpriv --pdeathsig KILL $julia --startup-file=no --history-file=no --compiled-modules=existing --project=$project $child $args`,
        ["PATH"=>ENV["PATH"], "JULIA_DEPOT_PATH"=>join(DEPOT_PATH,':'),
         "JULIA_LOAD_PATH"=>"@:@stdlib", "JULIA_NUM_THREADS"=>"1", "OPENBLAS_NUM_THREADS"=>"1",
         "LCM_EXECUTOR_PARENT_PID"=>string(getpid())])
end

function scientific_fixture(directory; capacity=2, ttl=300)
    project = normpath(joinpath(@__DIR__, "..", "..", "worker", "core"))
    fingerprint = native_environment_fingerprint(project)
    profiles = ProfileRegistry()
    register!(profiles, ProfileDefinition("fixture", project, fingerprint.digest;
        operations=("fixture.echo", "fixture.delay", "fixture.evict", "fixture.fail"),
        budget=ResourceBudget(prepare_seconds=30, job_seconds=30)))
    commands = Dict("fixture"=>fixture_command(project, joinpath(@__DIR__, "scientific_child.jl"), string(ttl)))
    driver = ScientificDriverFixture(profiles, commands,
        Dict{String,Tuple{RT.Protocol.AssignmentFence,RT.ExecutionCore.ExecutorSupervisor}}(), 0, 0, false, true, false)
    resources = ScientificResources(driver)
    config = AgentConfig("worker-a", BrokerEndpoint("tls://broker.invalid", "/unused-password"),
        profiles, joinpath(directory, "owned"); capacity)
    clock = Ref(0.0)
    agent = AgentService(config, resources; clock=()->clock[])
    coordinator = string(uuid4())
    probe = RT.Protocol.WorkerProbe("2.0", "worker-a", coordinator, string(uuid4()))
    receive_probe!(agent.ledger, probe)
    fences = [RT.Protocol.AssignmentFence(string(uuid4()), string(uuid4()), "owner-$i", "main",
        "worker-a", agent.ledger.boot_id, coordinator, "fixture", "1.0.0", fingerprint.digest, 1) for i in 1:capacity]
    for fence in fences
        handle_lease_control!(agent.ledger, lease_command(fence))
    end
    return (; driver, resources, agent, clock, fences)
end

function assigned_fixture_job(resources,fence, operation="fixture.echo", parameters=Dict("value"=>3); timeout=Dates.Second(30), engine_constraint=nothing)
    request = RT.Protocol.new_job_request(operation, parameters; session_id=fence.run_id, timeout, engine_constraint)
    target = try
        prepared_execution(resources,fence)
    catch error
        error isa AccessDenied || rethrow()
        # Explicit synthetic target for pre-preparation rejection tests only.
        RT.Protocol.PreparedExecution(string(uuid4()),1,repeat("b",64))
    end
    RT.Protocol.AssignedJob("2.0", fence, request,target)
end
