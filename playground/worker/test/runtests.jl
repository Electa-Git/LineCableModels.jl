using Dates
using LineCableModelsWorker
using LineCableModelsPlaygroundProtocol
using Logging
using Test

const Worker = LineCableModelsWorker

function test_context()
    events = Any[]
    context = Worker.ExecutionContext(
        "test-job",
        Worker.CancellationToken(),
        (progress, stage, message) -> push!(events, (progress, stage, message)),
        message -> push!(events, message),
        String[],
        Worker.PreparedResourceCache(),
        nothing
    )
    return context, events
end

@testset "closed operation registry" begin
    registry = Worker.default_registry()
    @test all(name -> haskey(registry.operations, name), (
        "system.echo",
        "system.delay",
        "system.executor_delay",
        "system.executor_warning",
        "system.fail",
        "system.progress",
        "cable.geometry_summary",
        "cable.constants",
        "cable.frequency_sweep",
        "line.frequency_scan",
        "powerflow.prepare",
        "impedance.evaluate",
    ))
    @test_throws Worker.PermanentOperationError Worker.registered_operation(
        registry,
        "system.eval"
    )
    @test occursin(r"^[0-9a-f]{64}$", Worker.environment_fingerprint())
    @test occursin(r"^0\.2\.0\+source\.[0-9a-f]{12}$", Worker.declared_engine_version())
    @test !isdefined(Worker, :nats_client_url)
    @test Worker.retryable_runtime_setup_error(
        Worker.NATS.JetStream.ApiError(code=404, description="stream not found")
    )
    @test Worker.retryable_runtime_setup_error(
        Worker.NATS.NATSError(503, "temporarily unavailable")
    )
    @test !Worker.retryable_runtime_setup_error(
        Worker.NATS.JetStream.ApiError(code=403, description="forbidden")
    )
    setup_attempts = Ref(0)
    recovered = Logging.with_logger(Logging.NullLogger()) do
        Worker.retry_runtime_setup(Threads.Atomic{Bool}(false); wait=_ -> nothing) do
            setup_attempts[] += 1
            setup_attempts[] < 3 && throw(Worker.NATS.JetStream.ApiError(
                code=404,
                description="stream not found"
            ))
            :ready
        end
    end
    @test recovered == :ready
    @test setup_attempts[] == 3
    @test_throws Worker.NATS.JetStream.ApiError Worker.retry_runtime_setup(
        Threads.Atomic{Bool}(false);
        wait=_ -> nothing
    ) do
        throw(Worker.NATS.JetStream.ApiError(code=403, description="forbidden"))
    end

    validator = parameters -> parameters
    executor = (_, parameters) -> parameters
    @test_throws ArgumentError Worker.OperationSpec(
        "system.echo",
        validator,
        executor;
        schema_version="latest"
    )
    @test_throws ArgumentError Worker.OperationSpec(
        "system.echo",
        validator,
        executor;
        timeout_seconds=24 * 60 * 60 + 1
    )
    @test_throws ArgumentError Worker.OperationSpec(
        "system.echo",
        validator,
        executor;
        capability="invalid capability"
    )
end

@testset "supervised cancellation and recovery" begin
    registry = Worker.default_registry()
    spec = Worker.registered_operation(registry, "system.executor_delay")
    supervisor = Worker.ExecutorSupervisor()
    canceled_context, _ = test_context()
    timer = Timer(0.2) do _
        Worker.cancel!(canceled_context.token)
    end
    try
        @test_throws Worker.OperationCanceled Worker.execute_supervised!(
            supervisor,
            spec,
            canceled_context,
            Dict{String,Any}("seconds" => 10.0)
        )
        close(timer)
        canceled_generation = supervisor.generation

        recovered_context, _ = test_context()
        result = Worker.execute_supervised!(
            supervisor,
            spec,
            recovered_context,
            Dict{String,Any}("seconds" => 0.05)
        )
        @test result["elapsed_seconds"] == 0.05
        @test supervisor.generation == canceled_generation + 1

        warning_context, warning_events = test_context()
        warning = Worker.execute_supervised!(
            supervisor,
            Worker.registered_operation(registry, "system.executor_warning"),
            warning_context,
            Dict{String,Any}("message" => "expected scientific warning")
        )
        @test warning["message"] == "expected scientific warning"
        @test warning_context.warnings == ["expected scientific warning"]
        @test any(
            event -> event == "Warning: expected scientific warning",
            warning_events
        )

        deadline_context, _ = test_context()
        deadline_context.deadline = Dates.now(Dates.UTC) + Millisecond(100)
        @test_throws Worker.OperationDeadlineExpired Worker.execute_supervised!(
            supervisor,
            spec,
            deadline_context,
            Dict{String,Any}("seconds" => 10.0)
        )
    finally
        isopen(timer) && close(timer)
        Worker.stop_executor!(supervisor)
    end
end

@testset "executor framing isolates engine output" begin
    frame = Worker.decode_executor_line(
        Worker.EXECUTOR_FRAME_PREFIX * "{\"type\":\"result\",\"result\":{\"x\":1}}"
    )
    @test frame.type == "result"
    @test frame.result.x == 1

    engine_output = Worker.decode_executor_line("[ PowerModels | Warn]: diagnostic")
    @test engine_output.type == "engine_output"
    @test engine_output.message == "[ PowerModels | Warn]: diagnostic"
end

@testset "diagnostic operations" begin
    registry = Worker.default_registry()
    context, events = test_context()
    echo = Worker.registered_operation(registry, "system.echo")
    result = Worker.execute_operation(
        echo,
        context,
        Dict{String,Any}("message" => "hello")
    )
    @test result["echo"]["message"] == "hello"
    @test !isempty(events)

    failure = Worker.registered_operation(registry, "system.fail")
    @test_throws Worker.PermanentOperationError Worker.execute_operation(
        failure,
        context,
        Dict{String,Any}("message" => "expected")
    )
end

@testset "prepared resources are single-flight" begin
    cache = Worker.PreparedResourceCache(ttl_seconds=10, max_entries=2)
    builds = Threads.Atomic{Int}(0)
    builder = () -> begin
        builds[] += 1
        sleep(0.05)
        return :resource
    end
    tasks = [@async Worker.prepare_resource!(builder, cache, "same-key") for _ in 1:4]
    @test fetch.(tasks) == fill(:resource, 4)
    @test builds[] == 1
    @test Worker.prepared_status(cache, "same-key") == :hot

    failure_builds = Threads.Atomic{Int}(0)
    failing_builder = () -> begin
        failure_builds[] += 1
        error("expected preparation failure")
    end
    @test_throws TaskFailedException Worker.prepare_resource!(
        failing_builder,
        cache,
        "failed-key"
    )
    @test Worker.prepared_status(cache, "failed-key") == :failed
    @test_throws TaskFailedException Worker.prepare_resource!(
        failing_builder,
        cache,
        "failed-key"
    )
    @test failure_builds[] == 1

    bounded = Worker.PreparedResourceCache(ttl_seconds=10, max_entries=2)
    Worker.prepare_resource!(() -> :one, bounded, "one")
    Worker.prepare_resource!(() -> :two, bounded, "two")
    Worker.prepare_resource!(() -> :three, bounded, "three")
    @test length(bounded.entries) == 2

    expiring = Worker.PreparedResourceCache(ttl_seconds=0.02, max_entries=1)
    Worker.prepare_resource!(() -> :brief, expiring, "brief")
    sleep(0.03)
    @test Worker.prepared_status(expiring, "brief") == :cold
end

@testset "result and artifact stores" begin
    mktempdir() do directory
        result_cache = Worker.ResultCache(joinpath(directory, "results"))
        request = LineCableModelsPlaygroundProtocol.new_job_request(
            "system.echo",
            Dict("message" => "cache");
            session_id="worker-test"
        )
        key = Worker.result_cache_key(request, "1.0", "1.0", repeat("f", 64))
        cached = Worker.CachedResult(Dict("value" => 7), ["warning"])
        Worker.cache_put!(result_cache, key, cached)
        loaded = Worker.cache_get(result_cache, key)
        @test loaded.result == Dict("value" => 7)
        @test loaded.warnings == ["warning"]

        artifacts = Worker.ArtifactStore(
            joinpath(directory, "artifacts");
            max_inline_bytes=16
        )
        small = Worker.split_result(artifacts, Dict("x" => 1))
        @test !isnothing(small.inline)
        @test isnothing(small.artifact)
        large = Worker.split_result(artifacts, Dict("payload" => repeat("x", 100)))
        @test isnothing(large.inline)
        @test large.artifact.retrieval_reference ==
              "/artifacts/sha256/$(large.artifact.sha256)"
        @test isfile(joinpath(artifacts.directory, large.artifact.sha256))
        @test isfile(joinpath(
            artifacts.directory,
            "$(large.artifact.sha256).metadata.json"
        ))
        @test large.artifact.artifact_id == "sha256:$(large.artifact.sha256)"
        @test_throws ArgumentError Worker.store_artifact!(
            artifacts,
            UInt8[0x01];
            media_type="text/plain\r\nX-Injected: true"
        )

        @test_throws ArgumentError Worker.S3ArtifactStore(
            "http://127.0.0.1:9000",
            "artifacts",
            "writer",
            "secret"
        )
        remote = Worker.S3ArtifactStore(
            "http://127.0.0.1:9000",
            "artifacts",
            "writer",
            "secret";
            allow_insecure=true,
            prefix="test"
        )
        @test Worker.artifact_object_key(remote, "sha256", repeat("a", 64)) ==
              "test/sha256/$(repeat("a", 64))"
        @test_throws ArgumentError Worker.S3ArtifactStore(
            "https://objects.example.invalid",
            "artifacts",
            "writer",
            "secret";
            prefix="../escape"
        )

        attempts = Ref(0)
        @test Worker.with_s3_transport_retry(attempts=3) do
            attempts[] += 1
            attempts[] < 3 && error("transient transport failure")
            :recovered
        end == :recovered
        @test attempts[] == 3
        @test_throws ArgumentError Worker.with_s3_transport_retry(
            () -> nothing;
            attempts=0
        )
    end
end

@testset "scientific cache projection" begin
    registry = Worker.default_registry()
    specification = Dict{String,Any}(
        "case_id" => "ohl_ugc_transition_v1",
        "earth_resistivity_ohm_m" => 100.0,
        "corridor_length_m" => 100_000.0,
    )
    base = Dict{String,Any}(
        "specification" => specification,
        "prepared_resource_key" => "worker-a-resource",
        "ugc_share" => 0.5,
        "corridor_length_m" => 100_000.0,
        "length_error_percent" => 5.0,
        "minimum_frequency_hz" => 100.0,
        "maximum_frequency_hz" => 2_500.0,
        "frequency_points" => 10,
    )
    other_hint = deepcopy(base)
    other_hint["prepared_resource_key"] = "worker-b-resource"
    spec = Worker.registered_operation(registry, "impedance.evaluate")
    @test spec.cache_key_inputs(base) == spec.cache_key_inputs(other_hint)

    short_specification = deepcopy(specification)
    short_specification["corridor_length_m"] = 25_000.0
    long_specification = deepcopy(specification)
    long_specification["corridor_length_m"] = 250_000.0
    short_preparation = Worker.validate_powerflow_spec(Dict(
        "specification" => short_specification
    ))
    long_preparation = Worker.validate_powerflow_spec(Dict(
        "specification" => long_specification
    ))
    @test Worker.powerflow_resource_key(short_preparation["specification"]) ==
          Worker.powerflow_resource_key(long_preparation["specification"])

    longer = deepcopy(base)
    longer["corridor_length_m"] = 250_000.0
    @test spec.cache_key_inputs(base) != spec.cache_key_inputs(longer)
end

@testset "real engine stays behind supervised executor" begin
    registry = Worker.default_registry()
    spec = Worker.registered_operation(registry, "cable.geometry_summary")
    context, _ = test_context()
    supervisor = Worker.ExecutorSupervisor()
    try
        result = Worker.execute_supervised!(
            supervisor,
            spec,
            context,
            Dict{String,Any}()
        )
        @test result["cable_id"] == "playground-coaxial"
        @test result["terminals"] == ["core", "sheath"]
        @test result["outer_radius_m"] ≈ 0.016
    finally
        Worker.stop_executor!(supervisor)
    end
end

include("scientific.jl")
