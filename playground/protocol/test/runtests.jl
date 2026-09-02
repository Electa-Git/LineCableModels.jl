using Dates
using JSON3
using LineCableModelsPlaygroundProtocol
using Test

@testset "wire contract" begin
    request = new_job_request(
        "system.echo",
        Dict("message" => "hello", "values" => [1, 2.0]);
        session_id="session-1",
        timeout=Second(30)
    )
    @test decode_job_request(encode_message(request)) == request
    @test input_hash("system.echo", Dict("b" => 2, "a" => 1)) ==
          input_hash("system.echo", Dict("a" => 1, "b" => 2))
    @test operation_subject_token("line.frequency_scan") == "line_frequency_scan"
    @test_throws ArgumentError operation_subject_token(
        "line." * repeat("a", 124)
    )

    bad_version = JobRequest(
        "99.0",
        request.job_id,
        request.operation,
        request.session_id,
        request.submitted_at,
        request.deadline,
        request.input_hash,
        request.parameters,
        nothing,
        "normal"
    )
    @test_throws ArgumentError validate(bad_version)
    @test_throws ArgumentError new_job_request(
        "eval" , Dict("code" => "1 + 1"); session_id="session-1"
    )
    @test_throws ArgumentError new_job_request(
        "system.echo", Dict("callable" => identity); session_id="session-1"
    )
    @test_throws ArgumentError new_job_request(
        "system.echo", Dict("message" => "late");
        session_id="session-1",
        timeout=Second(0)
    )
    @test_throws ArgumentError new_job_request(
        "system.echo", Dict("message" => "too long");
        session_id="session-1",
        timeout=Hour(25)
    )
    @test_throws ArgumentError new_job_request(
        "system.echo", Dict("message" => "priority");
        session_id="session-1",
        priority="urgent"
    )
    @test_throws ArgumentError new_job_request(
        "a.$(repeat("b", 127))", Dict{String,Any}();
        session_id="session-1"
    )
    @test_throws ArgumentError normalize_wire(repeat("x", 128 * 1024 + 1))
    @test_throws ArgumentError normalize_wire(Dict(repeat("k", 257) => 1))
    @test_throws ArgumentError normalize_wire(Dict{Any,Any}(:same => 1, "same" => 2))
    @test_throws ArgumentError normalize_wire(Dict(1 => "non-string key"))
    @test_throws ArgumentError normalize_wire(Dict("value" => Inf))
    @test_throws ArgumentError decode_job_request(fill(UInt8(' '), 256 * 1024 + 1))
    @test_throws Exception decode_job_request("{")
    @test parse_utc_timestamp("2026-01-02T03:04:05.006Z") ==
          DateTime(2026, 1, 2, 3, 4, 5, 6)
    @test_throws ArgumentError parse_utc_timestamp("2026-01-02T03:04:05Z")
    @test_throws ArgumentError parse_utc_timestamp("2026-01-02T03:04:05.006+00:00")

    golden_path = joinpath(@__DIR__, "golden", "job-request-v1.json")
    golden = decode_job_request(read(golden_path, String))
    @test encode_message(golden) == chomp(read(golden_path, String))

    tampered = JobRequest(
        request.protocol_version,
        request.job_id,
        request.operation,
        request.session_id,
        request.submitted_at,
        request.deadline,
        request.input_hash,
        Dict{String,Any}("message" => "tampered"),
        request.engine_constraint,
        request.priority
    )
    @test_throws ArgumentError validate(tampered)

    invalid_priority = JobRequest(
        request.protocol_version,
        request.job_id,
        request.operation,
        request.session_id,
        request.submitted_at,
        request.deadline,
        request.input_hash,
        request.parameters,
        request.engine_constraint,
        "urgent"
    )
    @test_throws ArgumentError validate(invalid_priority)
end

@testset "events and results" begin
    event = JobEvent(
        PROTOCOL_VERSION,
        "JobProgress",
        "5ec88436-5e75-4c0e-9434-0e4abbec0f87",
        3,
        utc_timestamp(),
        "worker-1",
        0.5,
        "waiting",
        "Halfway",
        Dict{String,Any}()
    )
    @test decode_job_event(encode_message(event)) == event
    @test !isterminal(event)
    invalid_progress = JobEvent(
        event.protocol_version,
        event.event_type,
        event.job_id,
        event.event_sequence,
        event.timestamp,
        event.worker_id,
        1.1,
        event.stage,
        event.message,
        event.details
    )
    @test_throws ArgumentError validate(invalid_progress)
    invalid_stage = JobEvent(
        event.protocol_version,
        event.event_type,
        event.job_id,
        event.event_sequence,
        event.timestamp,
        event.worker_id,
        event.progress,
        "not a stage",
        event.message,
        event.details
    )
    @test_throws ArgumentError validate(invalid_stage)

    result = JobResult(
        PROTOCOL_VERSION,
        event.job_id,
        "system.echo",
        "1.0",
        repeat("a", 64),
        "diagnostic",
        repeat("b", 64),
        "worker-1",
        "miss",
        utc_timestamp(),
        utc_timestamp(),
        Dict{String,Any}("message" => "hello"),
        nothing,
        nothing,
        String[]
    )
    @test decode_job_result(encode_message(result)) == result

    artifact = ArtifactReference(
        "sha256:$(repeat("c", 64))",
        "application/json",
        512,
        repeat("c", 64),
        "s3",
        "/artifacts/sha256/$(repeat("c", 64))"
    )
    @test validate(artifact) == artifact
    @test_throws ArgumentError validate(ArtifactReference(
        artifact.artifact_id,
        "not-a-media-type",
        artifact.size,
        artifact.sha256,
        artifact.storage_backend,
        artifact.retrieval_reference
    ))
    @test_throws ArgumentError validate(ArtifactReference(
        artifact.artifact_id,
        artifact.media_type,
        artifact.size,
        artifact.sha256,
        artifact.storage_backend,
        "https://objects.example.invalid/$(artifact.sha256)"
    ))
    @test_throws ArgumentError validate(ArtifactReference(
        artifact.artifact_id,
        artifact.media_type,
        artifact.size,
        artifact.sha256,
        artifact.storage_backend,
        "/artifacts/sha256/$(repeat("d", 64))"
    ))

    failure = FailureInfo(
        "invalid_input",
        "The request is invalid",
        "c3cc35f7-4a0c-46ab-a7ab-1ae7169ca920",
        false
    )
    @test validate(failure) == failure
    @test_throws ArgumentError validate(FailureInfo(
        failure.category,
        "",
        failure.diagnostic_id,
        failure.retryable
    ))

    invalid_payload_count = JobResult(
        result.protocol_version,
        result.job_id,
        result.operation,
        result.schema_version,
        result.input_hash,
        result.engine_version,
        result.environment_fingerprint,
        result.worker_id,
        result.cache_status,
        result.started_at,
        result.completed_at,
        result.inline_result,
        artifact,
        nothing,
        result.warnings
    )
    @test_throws ArgumentError validate(invalid_payload_count)
    excessive_warning = JobResult(
        result.protocol_version,
        result.job_id,
        result.operation,
        result.schema_version,
        result.input_hash,
        result.engine_version,
        result.environment_fingerprint,
        result.worker_id,
        result.cache_status,
        result.started_at,
        result.completed_at,
        result.inline_result,
        nothing,
        nothing,
        [repeat("w", 4 * 1024 + 1)]
    )
    @test_throws ArgumentError validate(excessive_warning)

    heartbeat = WorkerHeartbeat(
        PROTOCOL_VERSION,
        "worker-1",
        utc_timestamp(),
        ["system.echo", "cable.constants"],
        "0.1.0",
        "0.2.0",
        2
    )
    @test decode_worker_heartbeat(encode_message(heartbeat)) == heartbeat
    @test_throws ArgumentError validate(WorkerHeartbeat(
        heartbeat.protocol_version,
        heartbeat.worker_id,
        heartbeat.timestamp,
        ["system.echo", "system.echo"],
        heartbeat.worker_version,
        heartbeat.engine_version,
        heartbeat.available_slots
    ))
    @test_throws ArgumentError validate(WorkerHeartbeat(
        heartbeat.protocol_version,
        heartbeat.worker_id,
        heartbeat.timestamp,
        heartbeat.capabilities,
        heartbeat.worker_version,
        heartbeat.engine_version,
        1_025
    ))

    for (name, decoder) in (
            "job-event-v1.json" => decode_job_event,
            "job-result-v1.json" => decode_job_result,
            "worker-heartbeat-v1.json" => decode_worker_heartbeat,
        )
        golden_path = joinpath(@__DIR__, "golden", name)
        golden = decoder(read(golden_path, String))
        @test encode_message(golden) == chomp(read(golden_path, String))
    end
end
