function admin_connection()
    url = get(
        ENV,
        "NATS_ADMIN_URL",
        get(ENV, "NATS_CONNECT_URL", "nats://127.0.0.1:4222")
    )
    return NATS.connect(
        url;
        name="linecablemodels-runtime-admin",
        retry_on_init_fail=false,
        send_enqueue_when_disconnected=false
    )
end

function job_stream_configuration()
    return NATS.JetStream.StreamConfiguration(
        name=JOB_STREAM,
        description="Durable LineCableModels calculation requests",
        subjects=["lcm.jobs.v1.>"],
        retention=:workqueue,
        max_age=Int64(Dates.value(Dates.Nanosecond(Dates.Hour(24)))),
        max_msg_size=Int32(256 * 1024),
        storage=:file,
        num_replicas=1,
        duplicate_window=Int64(Dates.value(Dates.Nanosecond(Dates.Minute(10))))
    )
end

function result_stream_configuration()
    return NATS.JetStream.StreamConfiguration(
        name=RESULT_STREAM,
        description="Durable LineCableModels terminal results",
        subjects=["lcm.results.v1.>"],
        retention=:limits,
        max_age=Int64(Dates.value(Dates.Nanosecond(Dates.Day(7)))),
        max_msgs_per_subject=1,
        max_msg_size=Int32(256 * 1024),
        storage=:file,
        num_replicas=1,
        duplicate_window=Int64(Dates.value(Dates.Nanosecond(Dates.Minute(10)))),
        allow_direct=true
    )
end

function initialize_runtime()
    connection = admin_connection()
    try
        jobs = NATS.JetStream.stream_update_or_create(
            connection,
            job_stream_configuration()
        )
        results = NATS.JetStream.stream_update_or_create(
            connection,
            result_stream_configuration()
        )
        println("JetStream runtime initialized")
        println("  $(jobs.config.name): $(jobs.config.retention) retention")
        println("  $(results.config.name): $(results.config.retention) retention")
    finally
        NATS.drain(connection)
    end
    return nothing
end

function runtime_status()
    connection = admin_connection()
    try
        println("NATS runtime: connected")
        for stream_name in (JOB_STREAM, RESULT_STREAM)
            try
                info = NATS.JetStream.stream_info(connection, stream_name)
                println(
                    "  $stream_name: $(info.state.messages) messages · " *
                    "$(info.state.consumer_count) consumers"
                )
            catch error
                println("  $stream_name: unavailable ($(sprint(showerror, error)))")
            end
        end
    finally
        NATS.drain(connection)
    end
    return nothing
end
