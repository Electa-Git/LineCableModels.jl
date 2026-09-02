using UUIDs
using NATS

function await_test(predicate, timeout_seconds, description)
    timedwait(predicate, timeout_seconds; pollint=0.05) == :ok || error(
        "Timed out waiting for $description"
    )
end

function await_terminal(handle; timeout_seconds=20)
    outcome = timedwait(
        () -> handle.state[] in LineCableModelsPlayground.JOB_TERMINAL_STATES,
        timeout_seconds;
        pollint=0.05
    )
    outcome == :ok || error(
        "Timed out waiting for terminal job state " *
        "(job=$(handle.job_id[]), state=$(handle.state[]), stage=$(handle.stage[]))"
    )
    return handle
end

function submit_test_job(
        client,
        operation,
        parameters;
        timeout=Second(30),
        priority="normal"
    )
    request = new_job_request(
        operation,
        parameters;
        session_id="integration-$(uuid4())",
        timeout,
        priority
    )
    handle = JobHandle()
    submit!(client, handle, request)
    return handle, request
end

function start_test_worker(sequence)
    launcher = normpath(joinpath(@__DIR__, "..", "lcm"))
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    command = `$launcher worker start --id integration-worker-$sequence`
    worker_environment = Pair{String,String}[
        "LCM_JULIA" => julia,
        "NATS_CONNECT_URL" => ENV["NATS_TEST_WORKER_URL"],
        "LCM_WORKER_DATA" => ENV["NATS_TEST_WORKER_DATA"],
        "LCM_HEARTBEAT_SECONDS" => "0.5",
        "LCM_ACK_GRACE_SECONDS" => "1"
    ]
    for (target, source) in (
        "NATS_TLS_CA_PATH" => "NATS_TEST_WORKER_TLS_CA_PATH",
        "NATS_TLS_CERT_PATH" => "NATS_TEST_WORKER_TLS_CERT_PATH",
        "NATS_TLS_KEY_PATH" => "NATS_TEST_WORKER_TLS_KEY_PATH",
        "NATS_TLS_SERVER_NAME" => "NATS_TEST_WORKER_TLS_SERVER_NAME",
    )
        haskey(ENV, source) && push!(worker_environment, target => ENV[source])
    end
    command = addenv(command, worker_environment...)
    log_path = joinpath(ENV["NATS_TEST_WORKER_DATA"], "worker-$sequence.log")
    log = open(log_path, "a")
    process = run(pipeline(command; stdout=log, stderr=log); wait=false)
    return (; process, log, log_path)
end

function stop_test_worker(worker; signal=Base.SIGINT, timeout_seconds=8)
    if process_running(worker.process)
        kill(worker.process, signal)
        timedwait(
            () -> !process_running(worker.process),
            timeout_seconds;
            pollint=0.05
        )
    end
    process_running(worker.process) && kill(worker.process, Base.SIGKILL)
    try
        wait(worker.process)
    catch
    end
    isopen(worker.log) && close(worker.log)
    return nothing
end

function direct_child_pids(process)
    path = "/proc/$(getpid(process))/task/$(getpid(process))/children"
    isfile(path) || return Int[]
    contents = strip(read(path, String))
    isempty(contents) && return Int[]
    return parse.(Int, split(contents))
end

process_exists(pid::Integer) = isdir("/proc/$pid")

function pause_test_broker()
    runtime = get(ENV, "NATS_TEST_CONTAINER_RUNTIME", "docker")
    container = ENV["NATS_TEST_CONTAINER"]
    run(Cmd([runtime, "pause", container]))
    return nothing
end

function unpause_test_broker()
    runtime = get(ENV, "NATS_TEST_CONTAINER_RUNTIME", "docker")
    container = ENV["NATS_TEST_CONTAINER"]
    run(Cmd([runtime, "unpause", container]))
    return nothing
end

@testset "broker lifecycle" begin
    client = BrokerClient(
        url=ENV["NATS_TEST_PUBLISHER_URL"],
        heartbeat_ttl=2
    )
    worker = nothing
    replacement_worker = nothing
    try
        await_test(() -> client.connection_status[] == :online, 10, "broker connection")

        queued_message = "queued-$(uuid4())"
        queued_handle, _ = submit_test_job(
            client,
            "system.echo",
            Dict("message" => queued_message);
            timeout=Second(45)
        )
        @test queued_handle.state[] == :unavailable
        @test occursin("no compatible worker", queued_handle.stage[])

        queued_cancel_handle, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 2.0);
            timeout=Second(45)
        )
        cancel!(client, queued_cancel_handle)
        @test queued_cancel_handle.state[] == :cancelling

        worker = start_test_worker(1)
        await_test(
            () -> any(
                operations -> "system.echo" in operations,
                values(client.available_capabilities[])
            ),
            30,
            "worker capability heartbeat"
        )
        await_terminal(queued_handle; timeout_seconds=30)
        @test queued_handle.state[] == :ready
        @test queued_handle.last_successful_result[].inline_result["echo"]["message"] ==
              queued_message
        await_terminal(queued_cancel_handle; timeout_seconds=30)
        @test queued_cancel_handle.state[] == :canceled

        blocker, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 1.0)
        )
        await_test(() -> blocker.state[] == :running, 10, "priority blocker")
        normal_handle, _ = submit_test_job(
            client,
            "system.delay",
            Dict("seconds" => 0.75);
            priority="normal"
        )
        high_handle, _ = submit_test_job(
            client,
            "system.echo",
            Dict("message" => "high-priority");
            priority="high"
        )
        await_terminal(high_handle; timeout_seconds=20)
        @test high_handle.state[] == :ready
        @test !(normal_handle.state[] in LineCableModelsPlayground.JOB_TERMINAL_STATES)
        await_terminal(normal_handle; timeout_seconds=20)
        @test normal_handle.state[] == :ready

        unique_message = "cache-$(uuid4())"
        first_handle, _ = submit_test_job(
            client,
            "system.echo",
            Dict("message" => unique_message)
        )
        await_terminal(first_handle)
        @test first_handle.state[] == :ready
        @test first_handle.last_successful_result[].cache_status == "miss"

        cached_handle, _ = submit_test_job(
            client,
            "system.echo",
            Dict("message" => unique_message)
        )
        await_terminal(cached_handle)
        @test cached_handle.state[] == :ready
        @test cached_handle.last_successful_result[].cache_status == "hit"

        progress_handle, _ = submit_test_job(
            client,
            "system.progress",
            Dict("steps" => 5, "interval_seconds" => 0.03)
        )
        await_terminal(progress_handle)
        @test progress_handle.state[] == :ready
        @test progress_handle.progress[] == 1.0
        @test any(
            entry -> occursin("diagnostic step", entry.text),
            progress_handle.messages[]
        )

        warning_handle, _ = submit_test_job(
            client,
            "system.executor_warning",
            Dict("message" => "integration scientific warning")
        )
        await_terminal(warning_handle)
        @test warning_handle.state[] == :ready
        @test warning_handle.last_successful_result[].warnings ==
              ["integration scientific warning"]
        @test any(
            entry -> occursin("Warning: integration scientific warning", entry.text),
            warning_handle.messages[]
        )

        dirty_request = new_job_request(
            "system.delay",
            Dict("seconds" => 0.2);
            session_id="dirty-during-submit-$(uuid4())",
            timeout=Second(30)
        )
        dirty_handle = JobHandle()
        submit_async!(client, dirty_handle, dirty_request)
        mark_dirty!(dirty_handle)
        await_test(
            () -> dirty_handle.state[] == :dirty &&
                !isnothing(dirty_handle.last_successful_result[]),
            20,
            "dirty result after asynchronous submission"
        )
        @test dirty_handle.state[] == :dirty
        @test dirty_handle.dirty_inputs[]
        @test !isnothing(dirty_handle.last_successful_result[])

        failed_handle, _ = submit_test_job(
            client,
            "system.fail",
            Dict("category" => "expected_failure", "message" => "expected")
        )
        await_terminal(failed_handle)
        @test failed_handle.state[] == :failed
        @test failed_handle.failure[].category == "expected_failure"

        canceled_handle, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 4.0);
            timeout=Second(15)
        )
        await_test(() -> canceled_handle.state[] == :running, 10, "supervised delay")
        cancel!(client, canceled_handle)
        await_terminal(canceled_handle)
        @test canceled_handle.state[] == :canceled
        @test canceled_handle.failure[].category == "canceled"

        timeout_handle, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 8.0);
            timeout=Second(15)
        )
        await_terminal(timeout_handle; timeout_seconds=15)
        @test timeout_handle.state[] == :failed
        @test timeout_handle.failure[].category == "timeout"

        deadline_handle, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 4.0);
            timeout=Second(1)
        )
        await_terminal(deadline_handle; timeout_seconds=10)
        @test deadline_handle.state[] == :expired
        @test deadline_handle.failure[].category == "deadline_expired"

        crash_handle, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 1.0);
            timeout=Second(90)
        )
        await_test(() -> crash_handle.state[] == :running, 10, "job before worker crash")
        await_test(
            () -> !isempty(direct_child_pids(worker.process)),
            5,
            "supervised executor child"
        )
        executor_pids = direct_child_pids(worker.process)
        stop_test_worker(worker; signal=Base.SIGKILL)
        worker = nothing
        await_test(
            () -> all(!process_exists(pid) for pid in executor_pids),
            5,
            "executor exit after daemon crash"
        )
        @test all(!process_exists(pid) for pid in executor_pids)
        await_test(
            () -> isempty(client.available_capabilities[]),
            6,
            "crashed worker heartbeat expiry"
        )
        replacement_worker = start_test_worker(2)
        await_test(
            () -> haskey(
                client.available_capabilities[],
                "integration-worker-2"
            ),
            30,
            "replacement worker heartbeat"
        )
        await_terminal(crash_handle; timeout_seconds=45)
        @test crash_handle.state[] == :ready

        recovery_handle, _ = submit_test_job(
            client,
            "system.echo",
            Dict("message" => "after-executor-and-worker-recovery")
        )
        await_terminal(recovery_handle)
        @test recovery_handle.state[] == :ready

        duplicate_request = new_job_request(
            "system.echo",
            Dict("message" => "deduplicated-$(uuid4())");
            session_id="duplicate-integration",
            timeout=Second(30)
        )
        connection = LineCableModelsPlayground.active_connection(client)
        duplicate_subject = LineCableModelsPlayground.job_subject(
            duplicate_request.operation
        )
        payload = encode_message(duplicate_request)
        first_ack = NATS.JetStream.stream_publish(
            connection,
            duplicate_subject,
            (payload, ["Nats-Msg-Id" => duplicate_request.job_id])
        )
        duplicate_ack = NATS.JetStream.stream_publish(
            connection,
            duplicate_subject,
            (payload, ["Nats-Msg-Id" => duplicate_request.job_id])
        )
        @test first_ack.duplicate in (false, nothing)
        @test duplicate_ack.duplicate == true

        surviving_handle, _ = submit_test_job(
            client,
            "system.executor_delay",
            Dict("seconds" => 1.0);
            timeout=Second(45)
        )
        await_test(
            () -> surviving_handle.state[] == :running,
            10,
            "in-flight job before publisher restart"
        )
        surviving_job_id = surviving_handle.job_id[]
        close!(client)
        client = BrokerClient(
            url=ENV["NATS_TEST_PUBLISHER_URL"],
            heartbeat_ttl=2
        )
        await_test(
            () -> client.connection_status[] == :online,
            10,
            "publisher restart broker connection"
        )
        surviving_attachment = JobHandle()
        reattach!(client, surviving_attachment, surviving_job_id)
        await_terminal(surviving_attachment; timeout_seconds=30)
        @test surviving_attachment.state[] == :ready
        @test surviving_attachment.last_successful_result[].inline_result[
            "elapsed_seconds"
        ] == 1.0

        if haskey(ENV, "NATS_TEST_CONTAINER")
            paused = false
            try
                pause_test_broker()
                paused = true
                await_test(
                    () -> client.connection_status[] != :online,
                    10,
                    "publisher broker-outage state"
                )
                offline_handle = JobHandle()
                offline_request = new_job_request(
                    "system.echo",
                    Dict("message" => "must-not-buffer");
                    session_id="offline-integration",
                    timeout=Second(30)
                )
                @test_throws LineCableModelsPlayground.BrokerUnavailable submit!(
                    client,
                    offline_handle,
                    offline_request
                )
            finally
                paused && unpause_test_broker()
            end
            await_test(
                () -> client.connection_status[] == :online,
                20,
                "publisher broker reconnection"
            )
            await_test(
                () -> haskey(
                    client.available_capabilities[],
                    "integration-worker-2"
                ),
                20,
                "worker heartbeat after broker reconnection"
            )
            reconnected_handle, _ = submit_test_job(
                client,
                "system.echo",
                Dict("message" => "after-broker-reconnect")
            )
            await_terminal(reconnected_handle)
            @test reconnected_handle.state[] == :ready
        end

        completed_job_id = first_handle.job_id[]
        close!(client)
        replacement = BrokerClient(
            url=ENV["NATS_TEST_PUBLISHER_URL"],
            heartbeat_ttl=2
        )
        try
            await_test(
                () -> replacement.connection_status[] == :online,
                10,
                "replacement publisher connection"
            )
            attached = JobHandle()
            reattach!(replacement, attached, completed_job_id)
            await_terminal(attached)
            @test attached.state[] == :ready
            @test attached.last_successful_result[].inline_result["echo"]["message"] ==
                  unique_message
        finally
            close!(replacement)
        end

        graceful_executor_pids = direct_child_pids(replacement_worker.process)
        @test !isempty(graceful_executor_pids)
        stop_test_worker(replacement_worker; signal=Base.SIGINT)
        @test success(replacement_worker.process)
        await_test(
            () -> all(!process_exists(pid) for pid in graceful_executor_pids),
            5,
            "executor exit after graceful daemon shutdown"
        )
        replacement_worker = nothing
    finally
        client.closed[] || close!(client)
        !isnothing(worker) && stop_test_worker(worker)
        !isnothing(replacement_worker) && stop_test_worker(replacement_worker)
    end
end


@testset "broker authorization" begin
    function forbidden_publish(url, subject)
        connection = NATS.connect(
            url;
            retry_on_init_fail=false,
            send_enqueue_when_disconnected=false
        )
        try
            return NATS.JetStream.stream_publish(connection, subject, "{}")
        finally
            try
                NATS.drain(connection)
            catch
            end
        end
    end

    @test_throws NATS.NATSError forbidden_publish(
        ENV["NATS_TEST_PUBLISHER_URL"],
        "lcm.results.v1.$(uuid4())"
    )
    @test_throws NATS.NATSError forbidden_publish(
        ENV["NATS_TEST_WORKER_URL"],
        "lcm.jobs.v1.system_echo"
    )
    if haskey(ENV, "NATS_TEST_INVALID_URL")
        @test_throws Exception NATS.connect(
            ENV["NATS_TEST_INVALID_URL"];
            retry_on_init_fail=false,
            send_enqueue_when_disconnected=false
        )
    end
end

if get(ENV, "NATS_TEST_TLS", "0") == "1"
    @testset "broker mutual TLS" begin
        common = (
            tls_ca_path=ENV["NATS_TLS_CA_PATH"],
            retry_on_init_fail=false,
            send_enqueue_when_disconnected=false,
        )

        @test_throws Exception NATS.connect(
            ENV["NATS_TEST_PUBLISHER_URL"];
            common...,
            tls_cert_path=nothing,
            tls_key_path=nothing,
            tls_server_name=ENV["NATS_TLS_SERVER_NAME"]
        )

        @test_throws Exception NATS.connect(
            ENV["NATS_TEST_PUBLISHER_URL"];
            common...,
            tls_cert_path=ENV["NATS_TLS_CERT_PATH"],
            tls_key_path=ENV["NATS_TLS_KEY_PATH"],
            tls_server_name="certificate-name-mismatch.invalid"
        )

        @test_throws Exception NATS.connect(
            ENV["NATS_TEST_PUBLISHER_URL"];
            common...,
            tls_cert_path=ENV["NATS_TLS_CERT_PATH"],
            tls_key_path=nothing,
            tls_server_name=ENV["NATS_TLS_SERVER_NAME"]
        )
    end
end
