using Logging

function test_result(job_id; value=1, cache_status="miss")
    return JobResult(
        PROTOCOL_VERSION,
        job_id,
        "system.echo",
        "1.0",
        repeat("a", 64),
        "0.2.0",
        repeat("b", 64),
        "worker-test",
        cache_status,
        utc_timestamp(),
        utc_timestamp(),
        Dict{String,Any}("value" => value),
        nothing,
        nothing,
        String[]
    )
end

@testset "job lifecycle" begin
    @test LineCableModelsPlayground.missing_durable_result_error(
        LineCableModelsPlayground.NATS.NATSError(404, "Message Not Found")
    )
    @test !LineCableModelsPlayground.missing_durable_result_error(
        LineCableModelsPlayground.NATS.NATSError(503, "Unavailable")
    )
    handle = JobHandle()
    @test handle.input_revision == 0
    mark_dirty!(handle)
    @test handle.input_revision == 1
    @test handle.dirty_inputs[]
    handle = JobHandle()
    job_id = "5ec88436-5e75-4c0e-9434-0e4abbec0f87"
    result = test_result(job_id)
    LineCableModelsPlayground.apply_result!(handle, result)
    @test handle.state[] == :ready
    @test handle.last_successful_result[] == result

    mark_dirty!(handle)
    @test handle.state[] == :dirty
    @test handle.last_successful_result[] == result

    progress = JobEvent(
        PROTOCOL_VERSION,
        "JobProgress",
        job_id,
        2,
        utc_timestamp(),
        "worker-test",
        0.5,
        "solving",
        "Halfway",
        Dict{String,Any}()
    )
    LineCableModelsPlayground.apply_event!(handle, progress)
    @test handle.state[] == :running
    @test handle.progress[] == 0.5
    @test handle.worker[] == "worker-test"

    stale = JobEvent(
        PROTOCOL_VERSION,
        "JobProgress",
        job_id,
        1,
        utc_timestamp(),
        "worker-test",
        0.1,
        "stale",
        "Out of order",
        Dict{String,Any}()
    )
    LineCableModelsPlayground.apply_event!(handle, stale)
    @test handle.progress[] == 0.5

    pending = JobHandle()
    pending.job_id[] = job_id
    pending.pending_input_hash[] = repeat("a", 64)
    LineCableModelsPlayground.set_job_state!(pending, :running, "running")
    mark_dirty!(pending)
    @test pending.dirty_inputs[]
    LineCableModelsPlayground.apply_result!(pending, result)
    @test pending.state[] == :dirty
    @test pending.last_successful_result[] == result

    client = BrokerClient(autostart=false)
    rendered = Ref{Union{Nothing,JobResult}}(nothing)
    panel = JobPanel(
        client;
        operation="system.echo",
        parameters=Dict("message" => "persistent renderer"),
        result_handler=value -> (rendered[] = value)
    )
    @test isnothing(rendered[])
    LineCableModelsPlayground.apply_result!(panel.handle, result)
    @test rendered[] === result

    failing_panel = JobPanel(
        client;
        operation="system.echo",
        parameters=Dict("message" => "renderer failure"),
        result_handler=_ -> error("render exploded")
    )
    with_logger(NullLogger()) do
        LineCableModelsPlayground.apply_result!(failing_panel.handle, result)
    end
    @test failing_panel.handle.state[] == :ready
    @test failing_panel.handle.last_successful_result[] === result
    @test any(
        entry -> occursin("Result rendering failed", entry.text),
        failing_panel.handle.messages[]
    )
end
