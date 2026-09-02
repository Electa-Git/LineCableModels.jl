const JOB_STREAM = "LCM_JOBS"
const RESULT_STREAM = "LCM_RESULTS"

function job_subject(operation, priority::AbstractString="normal")
    priority in PRIORITIES || throw(ArgumentError("unsupported job priority: $priority"))
    return "lcm.jobs.v1.$priority.$(operation_subject_token(operation))"
end
result_subject(job_id) = "lcm.results.v1.$job_id"
event_subject(job_id) = "lcm.events.v1.$job_id"
cancel_subject(job_id) = "lcm.control.v1.cancel.$job_id"
heartbeat_subject(worker_id) = "lcm.workers.v1.heartbeat.$worker_id"
capability_subject(worker_id) = "lcm.workers.v1.capabilities.$worker_id"

redact_nats_url(url::AbstractString) = replace(
    string(url),
    r"//([^:/@]+):([^@]+)@" => s"//\1:***@"
)

"""
    WorkerConfig(; nats_url, worker_id, capacity, heartbeat_seconds,
                 ack_grace_seconds, data_directory)

Configure one scientific worker daemon.
"""
struct WorkerConfig
    "NATS connection URL used only by the server-side worker."
    nats_url::String

    "Stable worker identity advertised through heartbeats."
    worker_id::String

    "Maximum number of concurrently claimed jobs."
    capacity::Int

    "Heartbeat interval \\[s\\]."
    heartbeat_seconds::Float64

    "Additional JetStream acknowledgement time beyond each operation limit \\[s\\]."
    ack_grace_seconds::Float64

    "Directory containing persistent result and artifact data."
    data_directory::String
end

function WorkerConfig(;
        nats_url::AbstractString=get(ENV, "NATS_CONNECT_URL", "nats://127.0.0.1:4222"),
        worker_id::AbstractString=get(ENV, "LCM_WORKER_ID", "$(Sockets.gethostname())-$(getpid())"),
        capacity::Integer=parse(Int, get(ENV, "LCM_WORKER_CAPACITY", "1")),
        heartbeat_seconds::Real=parse(Float64, get(ENV, "LCM_HEARTBEAT_SECONDS", "2")),
        ack_grace_seconds::Real=parse(
            Float64,
            get(ENV, "LCM_ACK_GRACE_SECONDS", "60")
        ),
        data_directory::AbstractString=get(
            ENV,
            "LCM_WORKER_DATA",
            joinpath(homedir(), ".local", "share", "linecablemodels-worker")
        )
    )
    capacity > 0 || throw(ArgumentError("worker capacity must be positive"))
    heartbeat_seconds > 0 || throw(ArgumentError("heartbeat interval must be positive"))
    ack_grace_seconds > 0 || throw(ArgumentError(
        "acknowledgement grace interval must be positive"
    ))
    return WorkerConfig(
        string(nats_url),
        string(worker_id),
        Int(capacity),
        Float64(heartbeat_seconds),
        Float64(ack_grace_seconds),
        abspath(data_directory)
    )
end

mutable struct WorkerRuntime
    config::WorkerConfig
    registry::OperationRegistry
    connection::Any
    result_cache::ResultCache
    prepared_cache::PreparedResourceCache
    artifacts::AbstractArtifactStore
    cancellations::CancellationRegistry
    environment_fingerprint::String
    engine_version::String
    sequence::Dict{String,Int}
    sequence_lock::ReentrantLock
    active_jobs::Threads.Atomic{Int}
    stopped::Threads.Atomic{Bool}
end

function environment_fingerprint()
    worker_root = dirname(Base.active_project())
    protocol_root = normpath(joinpath(worker_root, "..", "protocol"))
    candidates = String[
        joinpath(worker_root, "Project.toml"),
        joinpath(worker_root, "Manifest.toml"),
        joinpath(protocol_root, "Project.toml"),
    ]
    for source_root in (joinpath(worker_root, "src"), joinpath(protocol_root, "src"))
        isdir(source_root) || continue
        append!(candidates, Iterators.flatten(
            joinpath.(root, files)
            for (root, _, files) in walkdir(source_root)
        ))
    end
    files = sort!(filter(path -> isfile(path), candidates))
    material = join((
        relpath(path, worker_root) * "\0" * read(path, String)
        for path in files
    ), '\0')
    return bytes2hex(SHA.sha256(material))
end

function declared_engine_version()
    project = normpath(joinpath(@__DIR__, "..", "..", "..", "Project.toml"))
    document = TOML.parsefile(project)
    version = string(get(document, "version", "unknown"))
    source_root = joinpath(dirname(project), "src")
    source_files = sort!(filter(
        path -> endswith(path, ".jl"),
        collect(Iterators.flatten(
            joinpath.(root, files)
            for (root, _, files) in walkdir(source_root)
        ))
    ))
    material = join((
        relpath(path, dirname(project)) * "\0" * read(path, String)
        for path in vcat([project], source_files)
    ), '\0')
    digest = bytes2hex(SHA.sha256(material))
    return get(ENV, "LCM_ENGINE_VERSION", "$version+source.$(digest[1:12])")
end

function WorkerRuntime(config::WorkerConfig, registry::OperationRegistry, connection)
    mkpath(config.data_directory)
    return WorkerRuntime(
        config,
        registry,
        connection,
        ResultCache(joinpath(config.data_directory, "results")),
        PreparedResourceCache(),
        configured_artifact_store(config.data_directory),
        CancellationRegistry(),
        environment_fingerprint(),
        declared_engine_version(),
        Dict{String,Int}(),
        ReentrantLock(),
        Threads.Atomic{Int}(0),
        Threads.Atomic{Bool}(false)
    )
end

function next_sequence!(runtime::WorkerRuntime, job_id::String)
    return lock(runtime.sequence_lock) do
        sequence = get(runtime.sequence, job_id, 0) + 1
        runtime.sequence[job_id] = sequence
        sequence
    end
end

function emit_event!(
        runtime::WorkerRuntime,
        request::JobRequest,
        event_type::AbstractString;
        progress=nothing,
        stage=nothing,
        message=nothing,
        details=Dict{String,Any}()
    )
    event = JobEvent(
        PROTOCOL_VERSION,
        string(event_type),
        request.job_id,
        next_sequence!(runtime, request.job_id),
        utc_timestamp(),
        runtime.config.worker_id,
        isnothing(progress) ? nothing : Float64(progress),
        isnothing(stage) ? nothing : string(stage),
        isnothing(message) ? nothing : string(message),
        normalize_wire(details)
    )
    NATS.publish(runtime.connection, event_subject(request.job_id), encode_message(event))
    return event
end

function heartbeat(runtime::WorkerRuntime)
    slots = max(0, runtime.config.capacity - runtime.active_jobs[])
    return WorkerHeartbeat(
        PROTOCOL_VERSION,
        runtime.config.worker_id,
        utc_timestamp(),
        capabilities(runtime.registry),
        string(pkgversion(LineCableModelsWorker)),
        runtime.engine_version,
        slots
    )
end

function heartbeat_loop!(runtime::WorkerRuntime, queue)
    while !runtime.stopped[]
        try
            message = encode_message(heartbeat(runtime))
            NATS.publish(
                runtime.connection,
                heartbeat_subject(runtime.config.worker_id),
                message
            )
            NATS.publish(
                runtime.connection,
                capability_subject(runtime.config.worker_id),
                message
            )
        catch error
            if error isa InterruptException
                runtime.stopped[] = true
                isopen(queue) && close(queue)
                return nothing
            end
            runtime.stopped[] || @debug "Could not publish worker heartbeat" exception=error
        end
        sleep(runtime.config.heartbeat_seconds)
    end
    return nothing
end

function parse_deadline(value::AbstractString)
    text = endswith(value, 'Z') ? chop(value) : string(value)
    return DateTime(text, dateformat"yyyy-mm-ddTHH:MM:SS.sss")
end

function request_expired(request::JobRequest)
    return Dates.now(Dates.UTC) > parse_deadline(request.deadline)
end

function delivery_attempt(message)
    reply = message.reply_to
    isnothing(reply) && return 1
    tokens = split(reply, '.')
    length(tokens) >= 5 || return 1
    return something(tryparse(Int, tokens[5]), 1)
end

function durable_result_exists(runtime::WorkerRuntime, job_id::String)
    try
        NATS.JetStream.stream_message_get(
            runtime.connection,
            RESULT_STREAM,
            result_subject(job_id)
        )
        return true
    catch error
        if (error isa NATS.NATSError || error isa NATS.JetStream.ApiError) &&
                error.code == 404
            return false
        end
        rethrow()
    end
end

function publish_result!(runtime::WorkerRuntime, result::JobResult)
    NATS.JetStream.stream_publish(
        runtime.connection,
        result_subject(result.job_id),
        (encode_message(result), ["Nats-Msg-Id" => result.job_id])
    )
    return result
end

function successful_result(
        runtime,
        request,
        spec,
        cached,
        cache_status,
        started_at
    )
    payload = split_result(runtime.artifacts, cached.result)
    return JobResult(
        PROTOCOL_VERSION,
        request.job_id,
        request.operation,
        spec.schema_version,
        request.input_hash,
        runtime.engine_version,
        runtime.environment_fingerprint,
        runtime.config.worker_id,
        cache_status,
        started_at,
        utc_timestamp(),
        payload.inline,
        payload.artifact,
        nothing,
        cached.warnings
    )
end

function failed_result(runtime, request, spec, failure, started_at)
    return JobResult(
        PROTOCOL_VERSION,
        request.job_id,
        request.operation,
        spec.schema_version,
        request.input_hash,
        runtime.engine_version,
        runtime.environment_fingerprint,
        runtime.config.worker_id,
        "bypass",
        started_at,
        utc_timestamp(),
        nothing,
        nothing,
        failure,
        String[]
    )
end

function cancellation_context(runtime, request)
    token = activate_cancellation!(runtime.cancellations, request.job_id)
    progress_callback = (fraction, stage, message) -> emit_event!(
        runtime,
        request,
        "JobProgress";
        progress=fraction,
        stage,
        message
    )
    log_callback = message -> emit_event!(
        runtime,
        request,
        "JobLog";
        stage="log",
        message
    )
    context = ExecutionContext(
        request.job_id,
        token,
        progress_callback,
        log_callback,
        String[],
        runtime.prepared_cache,
        parse_deadline(request.deadline)
    )
    return context
end

function process_job!(
        runtime::WorkerRuntime,
        message,
        spec::OperationSpec,
        supervisor::ExecutorSupervisor
    )
    request = try
        decode_job_request(NATS.payload(message))
    catch error
        @error "Discarding invalid durable job" exception=(error, catch_backtrace())
        NATS.JetStream.consumer_ack(runtime.connection, message, "+TERM")
        return nothing
    end
    if request.operation != spec.name
        @error "Consumer received mismatched operation" expected=spec.name actual=request.operation
        NATS.JetStream.consumer_ack(runtime.connection, message, "+TERM")
        return nothing
    end
    durable_result_exists(runtime, request.job_id) && begin
        NATS.JetStream.consumer_ack(runtime.connection, message)
        return nothing
    end

    started_at = utc_timestamp()
    context = cancellation_context(runtime, request)
    runtime.active_jobs[] += 1
    emit_event!(runtime, request, "JobClaimed"; progress=0.0, stage="claimed")
    try
        if request_expired(request)
            failure = FailureInfo(
                "deadline_expired",
                "The calculation deadline expired before execution",
                string(uuid4()),
                false
            )
            publish_result!(runtime, failed_result(runtime, request, spec, failure, started_at))
            emit_event!(runtime, request, "JobExpired"; stage="expired",
                message=failure.message)
            NATS.JetStream.consumer_ack(runtime.connection, message)
            return nothing
        end

        if !isnothing(request.engine_constraint) &&
                request.engine_constraint != runtime.engine_version
            throw(PermanentOperationError(
                "incompatible_engine",
                "No worker matching the requested engine version accepted this job"
            ))
        end

        emit_event!(runtime, request, "JobStarted"; progress=0.0, stage="starting")
        cache_inputs = Base.invokelatest(spec.cache_key_inputs, request.parameters)
        key = result_cache_key(
            request,
            spec.schema_version,
            runtime.engine_version,
            runtime.environment_fingerprint;
            parameters=cache_inputs
        )
        cached = spec.cache_policy == :content ? cache_get(runtime.result_cache, key) : nothing
        cache_status = isnothing(cached) ? "miss" : "hit"
        if isnothing(cached)
            result = if spec.execution_mode == :supervised
                execute_supervised!(supervisor, spec, context, request.parameters)
            else
                execute_operation(spec, context, request.parameters)
            end
            cached = CachedResult(result, copy(context.warnings))
            spec.cache_policy == :content && cache_put!(runtime.result_cache, key, cached)
        end
        result = successful_result(
            runtime,
            request,
            spec,
            cached,
            spec.cache_policy == :none ? "bypass" : cache_status,
            started_at
        )
        publish_result!(runtime, result)
        emit_event!(runtime, request, "JobSucceeded"; progress=1.0,
            stage="complete", message="Calculation complete")
        NATS.JetStream.consumer_ack(runtime.connection, message)
    catch error
        if error isa OperationCanceled
            failure = FailureInfo(
                "canceled",
                "Calculation canceled",
                string(uuid4()),
                false
            )
            publish_result!(runtime, failed_result(runtime, request, spec, failure, started_at))
            emit_event!(runtime, request, "JobCanceled"; stage="canceled",
                message=failure.message)
            NATS.JetStream.consumer_ack(runtime.connection, message)
        elseif error isa OperationDeadlineExpired
            failure = FailureInfo(
                "deadline_expired",
                "The calculation deadline expired during execution",
                string(uuid4()),
                false
            )
            publish_result!(runtime, failed_result(runtime, request, spec, failure, started_at))
            emit_event!(runtime, request, "JobExpired"; stage="expired",
                message=failure.message)
            NATS.JetStream.consumer_ack(runtime.connection, message)
        elseif error isa InterruptException
            try
                NATS.JetStream.consumer_ack(runtime.connection, message, "-NAK")
            catch
            end
            rethrow()
        elseif error isa PermanentOperationError
            failure = FailureInfo(
                error.category,
                error.message,
                string(uuid4()),
                false
            )
            publish_result!(runtime, failed_result(runtime, request, spec, failure, started_at))
            emit_event!(runtime, request, "JobFailed"; stage="failed",
                message=failure.message)
            NATS.JetStream.consumer_ack(runtime.connection, message)
        else
            diagnostic_id = string(uuid4())
            @error "Retryable worker failure" diagnostic_id request.job_id exception=(
                error,
                catch_backtrace()
            )
            if delivery_attempt(message) >= 3
                failure = FailureInfo(
                    "worker_failure",
                    "Calculation failed after bounded retries",
                    diagnostic_id,
                    false
                )
                publish_result!(runtime,
                    failed_result(runtime, request, spec, failure, started_at))
                emit_event!(runtime, request, "JobFailed"; stage="failed",
                    message=failure.message,
                    details=Dict("diagnostic_id" => diagnostic_id))
                NATS.JetStream.consumer_ack(runtime.connection, message)
            else
                emit_event!(runtime, request, "JobLog"; stage="retrying",
                    message="Worker attempt failed; JetStream will redeliver",
                    details=Dict("diagnostic_id" => diagnostic_id))
                NATS.JetStream.consumer_ack(runtime.connection, message, "-NAK")
            end
        end
    finally
        runtime.active_jobs[] -= 1
        finish_cancellation!(runtime.cancellations, request.job_id)
    end
    return nothing
end

function ensure_consumer!(
        runtime::WorkerRuntime,
        spec::OperationSpec,
        priority::AbstractString
    )
    priority in PRIORITIES || throw(ArgumentError("unsupported job priority: $priority"))
    token = operation_subject_token(spec.name)
    config = NATS.JetStream.ConsumerConfiguration(
        durable_name="lcm_v1_$(priority)_$token",
        name="lcm_v1_$(priority)_$token",
        description="LineCableModels $priority priority workers for $(spec.name)",
        ack_policy=:explicit,
        ack_wait=Int64(ceil(
            (spec.timeout_seconds + runtime.config.ack_grace_seconds) * 1e9
        )),
        max_deliver=3,
        filter_subjects=[job_subject(spec.name, priority)],
        max_ack_pending=runtime.config.capacity
    )
    return NATS.JetStream.consumer_create(runtime.connection, config, JOB_STREAM)
end

function retryable_runtime_setup_error(error)
    return (error isa NATS.NATSError || error isa NATS.JetStream.ApiError) &&
           error.code in (404, 408, 503)
end

function retry_runtime_setup(setup, stopped; wait=sleep)
    attempts = 0
    while !stopped[]
        try
            return setup()
        catch error
            error isa InterruptException && rethrow()
            retryable_runtime_setup_error(error) || rethrow()
            attempts += 1
            if attempts == 1 || attempts % 20 == 0
                @warn "Worker runtime is not initialized yet; retrying" attempts=attempts exception=error
            end
            wait(min(2.0, 0.1 * 2.0^min(attempts - 1, 5)))
        end
    end
    throw(InterruptException())
end

function initialize_consumers!(runtime::WorkerRuntime, registry::OperationRegistry)
    return retry_runtime_setup(runtime.stopped) do
        Dict(
            (name, priority) => ensure_consumer!(runtime, spec, priority)
            for (name, spec) in registry.operations
            for priority in PRIORITIES
        )
    end
end

function next_scheduled_job!(runtime, consumers, registry, cursors)
    operation_names = sort!(collect(keys(registry.operations)))
    for priority in ("high", "normal")
        isempty(operation_names) && continue
        first_index = mod(cursors[priority], length(operation_names)) + 1
        for offset in 0:(length(operation_names) - 1)
            index = mod(first_index - 1 + offset, length(operation_names)) + 1
            name = operation_names[index]
            consumer = consumers[(name, priority)]
            message = NATS.JetStream.consumer_next(
                runtime.connection,
                consumer;
                no_wait=true,
                no_throw=true
            )
            NATS.has_error_status(message) && continue
            cursors[priority] = index
            return (message, registry.operations[name])
        end
    end
    return nothing
end

function scheduling_loop!(runtime, consumers, registry, queue, capacity_limiter)
    cursors = Dict(priority => 0 for priority in PRIORITIES)
    while !runtime.stopped[]
        acquired = false
        try
            Base.acquire(capacity_limiter)
            acquired = true
            runtime.stopped[] && break
            item = next_scheduled_job!(runtime, consumers, registry, cursors)
            if isnothing(item)
                Base.release(capacity_limiter)
                acquired = false
                sleep(0.05)
                continue
            end
            put!(queue, item)
            acquired = false # execution loop releases this slot
        catch error
            if error isa InterruptException
                runtime.stopped[] = true
                isopen(queue) && close(queue)
                return nothing
            end
            runtime.stopped[] || @warn "Job scheduling failed" exception=error
            sleep(0.5)
        finally
            acquired && Base.release(capacity_limiter)
        end
    end
    return nothing
end

function execution_loop!(runtime, queue, capacity_limiter)
    supervisor = ExecutorSupervisor()
    try
        while !runtime.stopped[]
            item = try
                take!(queue)
            catch error
                error isa InvalidStateException && !isopen(queue) && break
                rethrow()
            end
            message, spec = item
            try
                process_job!(runtime, message, spec, supervisor)
            catch error
                if error isa InterruptException
                    runtime.stopped[] = true
                    isopen(queue) && close(queue)
                    return nothing
                end
                rethrow()
            finally
                Base.release(capacity_limiter)
            end
        end
    finally
        stop_executor!(supervisor)
    end
    return nothing
end

function run_worker(config::WorkerConfig=WorkerConfig())
    registry = default_registry()
    connection = NATS.connect(
        config.nats_url;
        name="linecablemodels-worker-$(config.worker_id)",
        ping_interval=2.0,
        max_pings_out=1,
        retry_on_init_fail=true,
        send_enqueue_when_disconnected=false
    )
    while NATS.status(connection) != NATS.CONNECTED
        sleep(0.25)
    end
    runtime = WorkerRuntime(config, registry, connection)
    NATS.subscribe(connection, cancel_subject("*")) do message
        job_id = strip(String(NATS.payload(message)))
        tryparse(UUID, job_id) === nothing && begin
            @warn "Ignoring malformed cancellation request"
            return nothing
        end
        request_cancellation!(runtime.cancellations, job_id)
        return nothing
    end
    consumers = initialize_consumers!(runtime, registry)
    queue = Channel{Any}(config.capacity)
    capacity_limiter = Base.Semaphore(config.capacity)
    errormonitor(@async heartbeat_loop!(runtime, queue))
    errormonitor(@async scheduling_loop!(
        runtime,
        consumers,
        registry,
        queue,
        capacity_limiter
    ))
    workers = [errormonitor(@async execution_loop!(runtime, queue, capacity_limiter))
        for _ in 1:config.capacity]
    println(
        "LineCableModels worker $(config.worker_id) connected to " *
        redact_nats_url(config.nats_url)
    )
    try
        wait(first(workers))
    finally
        runtime.stopped[] = true
        isopen(queue) && close(queue)
        try
            NATS.drain(connection)
        catch
        end
    end
    return nothing
end
