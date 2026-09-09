using Dates
include("scientific_driver_fixture.jl")
struct IncompleteScientificDriver <: AbstractScientificDriver end

@testset "scientific ownership is passive, explicit and independently prepared" begin
    @test_throws ArgumentError ScientificResources(IncompleteScientificDriver())
    mktempdir() do directory
        (; driver, resources, agent, clock, fences) = scientific_fixture(directory)
        first, second = fences
        try
            @test isempty(driver.processes) && driver.verified == 0
            @test_throws AccessDenied prepare_assigned!(resources, first, Dict{String,Any}())
            @test recover_owned!(resources) === nothing
            @test recover_owned!(resources) === nothing && driver.recovered == 1
            @test_throws ArgumentError bind_agent!(resources, AgentLeaseLedger("worker-a", driver.profiles))
            @test scientific_status(resources, first).preparation == :cold
            @test_throws AccessDenied execute_assigned!(resources, assigned_fixture_job(resources,first))
            @test isempty(driver.processes)
            @test_throws AccessDenied prepare_assigned!(resources,
                change_runtime_record(first; owner="foreign"), Dict{String,Any}())
            a = prepare_assigned!(resources, first, Dict("seconds"=>0.15))
            @test prepare_assigned!(resources, first, Dict("seconds"=>0.15)) === a
            @test_throws AccessDenied prepare_assigned!(resources, first, Dict("seconds"=>0.2))
            b = prepare_assigned!(resources, second, Dict{String,Any}())
            @test fetch(a)["cache_status"] == "miss"
            @test fetch(b)["cache_status"] == "miss"
            status = scientific_status(resources, first)
            @test status.preparation == :ready && status.phase == :idle
            @test status.generation == 1 && status.preparation_key !== nothing
            @test scientific_status(resources, second).executor_id != status.executor_id
            @test driver.processes[first.lease_id][2].process !== driver.processes[second.lease_id][2].process
            job = assigned_fixture_job(resources,first)
            @test job.execution==prepared_execution(resources,first)
            @test_throws AccessDenied execute_assigned!(resources,change_runtime_record(job;
                execution=change_runtime_record(job.execution;executor_generation=2)))
            @test_throws AccessDenied execute_assigned!(resources,change_runtime_record(job;
                execution=change_runtime_record(job.execution;preparation_key=repeat("b",64))))
            result_task = execute_assigned!(resources, job)
            @test fetch(result_task).value["value"] == 3
            @test fetch(result_task).schema_version=="1.2"
            @test execute_assigned!(resources, job) === result_task
            altered = change_runtime_record(job; request=change_runtime_record(job.request; priority="high"))
            @test_throws AccessDenied execute_assigned!(resources, altered)
            @test fetch(refresh_preparation!(resources, first))["ready"]
            @test_throws AccessDenied execute_assigned!(resources, assigned_fixture_job(resources,first, "foreign.operation"))
            @test_throws AccessDenied execute_assigned!(resources, assigned_fixture_job(resources,first; engine_constraint="wrong"))
            failure = execute_assigned!(resources, assigned_fixture_job(resources,first, "fixture.fail", Dict{String,Any}()))
            @test_throws TaskFailedException fetch(failure)
            @test scientific_status(resources, first).failure == "operation-rejected"
            @test scientific_status(resources, first).preparation == :ready
            @test !occursin("private", sprint(showerror, failure.exception))
            @test !occursin("secret diagnostic", sprint(showerror, failure.exception))
            @test !occursin("/private/fixture", sprint(showerror, TaskFailedException(failure)))
            @test fetch(execute_assigned!(resources, assigned_fixture_job(resources,first, "fixture.evict", Dict{String,Any}()))).value["evicted"]
            @test scientific_status(resources, first).preparation == :cold
            @test_throws AccessDenied execute_assigned!(resources, assigned_fixture_job(resources,first))
            @test fetch(prepare_assigned!(resources, first, Dict("seconds"=>0.15)))["cache_status"] == "miss"
            # Once both are prepared, a slow role cannot serialize its sibling.
            slow_job = assigned_fixture_job(resources,first, "fixture.delay", Dict("seconds"=>10.0))
            slow = execute_assigned!(resources, slow_job)
            quick = @elapsed fetch(execute_assigned!(resources, assigned_fixture_job(resources,second)))
            @test quick < 1.5
            @test !istaskdone(slow)
            tick_agent!(agent) # compile control pass before measuring responsiveness
            @test (@elapsed tick_agent!(agent)) < 0.1
            @test !cancel_assigned!(resources, first, string(uuid4()))
            @test cancel_assigned!(resources, first, slow_job.request.job_id)
            @test_throws TaskFailedException fetch(slow)
            @test scientific_status(resources, first).preparation != :ready
            @test scientific_status(resources, first).failure == "canceled"
            @test scientific_status(resources, second).preparation == :ready
            previous_id = scientific_status(resources, first).executor_id
            @test fetch(prepare_assigned!(resources, first, Dict{String,Any}()))["cache_status"] == "miss"
            @test scientific_status(resources, first).executor_id != previous_id
            @test_throws AccessDenied execute_assigned!(resources,job) # prior prepared child is gone
            @test driver.verified >= 8
            @test !any(id.name in ("LineCableModels", "PowerImpedance", "Bonito") for id in keys(Base.loaded_modules))
        finally
            close(agent)
        end
        @test driver.closed && isempty(driver.processes) && isempty(resources.handles)
    end
end

@testset "closed agent can finish previously unresolved owned cleanup" begin
    mktempdir() do directory
        (; driver, resources, agent, clock, fences) = scientific_fixture(directory; capacity=1)
        fence = only(fences)
        try
            recover_owned!(resources)
            fetch(prepare_assigned!(resources, fence, Dict{String,Any}()))
            driver.allow_release = false
            @test_throws ArgumentError close(agent)
            @test agent.closed && !agent.cleanup_complete
            @test haskey(resources.handles, fence.lease_id)
            @test_throws AccessDenied prepare_assigned!(resources, fence, Dict{String,Any}())
            driver.allow_release = true
            @test close(agent) === nothing
            @test agent.cleanup_complete && isempty(resources.handles) && isempty(driver.processes)
            @test close(agent) === nothing
        finally
            driver.allow_release = true
            close(agent)
        end
    end
end

@testset "lease loss cancels execution and cleanup cannot release capacity early" begin
    mktempdir() do directory
        (; driver, resources, agent, clock, fences) = scientific_fixture(directory; capacity=1)
        fence = only(fences)
        try
            recover_owned!(resources)
            fetch(prepare_assigned!(resources, fence, Dict{String,Any}()))
            process = driver.processes[fence.lease_id][2].process
            running = execute_assigned!(resources, assigned_fixture_job(resources,fence, "fixture.delay", Dict("seconds"=>10.0)))
            yield()
            driver.allow_release = false
            clock[] = 6.0 # original five-second grant expired, independently of the browser
            @test_throws TaskFailedException fetch(running)
            @test !process_running(process)
            @test_throws AccessDenied scientific_status(resources, fence)
            @test !release_owned!(resources, fence)
            @test haskey(resources.handles, fence.lease_id)
            @test only(values(agent.ledger.leases)).state == :closing
            driver.allow_release = true
            @test release_owned!(resources, fence)
            @test isempty(resources.handles) && isempty(driver.processes)
            complete_agent_cleanup!(agent.ledger, fence)
            @test only(values(agent.ledger.leases)).state == :closed
        finally
            driver.allow_release = true
            close(agent)
        end
    end
end

@testset "preflight failure never starts or revives a prepared process" begin
    mktempdir() do directory
        (; driver, resources, agent, clock, fences) = scientific_fixture(directory; capacity=1)
        fence = only(fences)
        try
            recover_owned!(resources)
            driver.fail_verify = true
            @test_throws TaskFailedException fetch(prepare_assigned!(resources, fence, Dict{String,Any}()))
            @test isempty(driver.processes)
            @test scientific_status(resources, fence).preparation == :failed
            driver.fail_verify = false
            fetch(prepare_assigned!(resources, fence, Dict{String,Any}()))
            process = driver.processes[fence.lease_id][2].process
            driver.fail_verify = true
            @test_throws TaskFailedException fetch(execute_assigned!(resources, assigned_fixture_job(resources,fence)))
            @test !process_running(process) && isempty(driver.processes)
            @test scientific_status(resources, fence).preparation != :ready
        finally
            driver.fail_verify = false
            close(agent)
        end
    end
end
