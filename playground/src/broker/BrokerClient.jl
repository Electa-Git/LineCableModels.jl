struct BrokerUnavailable <: Exception
    message::String
end

Base.showerror(io::IO, error::BrokerUnavailable) = print(io, error.message)

"""
    BrokerClient(; url=get(ENV, "NATS_CONNECT_URL", "nats://127.0.0.1:4222"),
                 heartbeat_ttl=10.0, autostart=true)

Create a non-blocking publisher-side connection manager. When `autostart` is
true, connection establishment occurs in a background task and cannot delay
page construction or server startup.
"""
mutable struct BrokerClient
    "NATS server URL."
    url::String

    "Observable connection state."
    connection_status::Observable{Symbol}

    "Observable active worker capabilities keyed by worker identity."
    available_capabilities::Observable{Dict{String,Vector{String}}}

    "Live worker heartbeats."
    workers::Dict{String,WorkerPresence}

    "Session-local handles retained for reconnection and reattachment."
    handles::Dict{String,JobHandle}

    "Underlying NATS connection."
    connection::Base.RefValue{Any}

    "Connection supervisor task."
    supervisor::Base.RefValue{Union{Nothing,Task}}

    "Heartbeat expiry [s]."
    heartbeat_ttl::Float64

    "Closed flag."
    closed::Threads.Atomic{Bool}

    "Lock protecting connections, workers, and handles."
    lock::ReentrantLock
end

function BrokerClient(;
        url::AbstractString=get(ENV, "NATS_CONNECT_URL", "nats://127.0.0.1:4222"),
        heartbeat_ttl::Real=10.0,
        autostart::Bool=true
    )
    heartbeat_ttl > 0 || throw(ArgumentError("heartbeat_ttl must be positive"))
    client = BrokerClient(
        string(url),
        Observable(:offline),
        Observable(Dict{String,Vector{String}}()),
        Dict{String,WorkerPresence}(),
        Dict{String,JobHandle}(),
        Ref{Any}(nothing),
        Ref{Union{Nothing,Task}}(nothing),
        Float64(heartbeat_ttl),
        Threads.Atomic{Bool}(false),
        ReentrantLock()
    )
    autostart && start!(client)
    return client
end

function prune_workers!(client::BrokerClient)
    cutoff = time() - client.heartbeat_ttl
    lock(client.lock) do
        filter!(pair -> last(pair).received_at >= cutoff, client.workers)
        client.available_capabilities[] = worker_capabilities(client.workers)
    end
    client.connection_status[] == :online && refresh_handle_availability!(client)
    return client
end

function register_heartbeat!(client::BrokerClient, message)
    heartbeat = decode_worker_heartbeat(NATS.payload(message))
    lock(client.lock) do
        client.workers[heartbeat.worker_id] = WorkerPresence(heartbeat, time())
        client.available_capabilities[] = worker_capabilities(client.workers)
    end
    refresh_handle_availability!(client)
    return nothing
end

function refresh_handle_availability!(client::BrokerClient)
    snapshot = lock(client.lock) do
        (copy(client.workers), collect(values(client.handles)))
    end
    workers, handles = snapshot
    for handle in handles
        request = handle.request
        isnothing(request) && continue
        compatible = supports_operation(
            workers,
            request.operation;
            engine_constraint=request.engine_constraint
        )
        lock(handle.lock) do
            if handle.state[] == :unavailable && compatible
                set_job_state!(handle, :queued, "Queued · compatible worker discovered")
            elseif handle.state[] == :queued && !compatible
                set_job_state!(handle, :unavailable,
                    "Queued · no compatible worker online")
            end
        end
    end
    return client
end

function mark_pending_offline!(client::BrokerClient)
    handles = lock(client.lock) do
        collect(values(client.handles))
    end
    for handle in handles
        lock(handle.lock) do
            if handle.state[] in (:submitting, :queued, :unavailable)
                set_job_state!(handle, :unavailable,
                    "Broker offline · durable request status unavailable")
            end
        end
    end
    return client
end

function missing_durable_result_error(error)
    return (error isa NATS.JetStream.ApiError || error isa NATS.NATSError) &&
           error.code == 404
end

function fetch_durable_result!(client::BrokerClient, handle::JobHandle)
    connection = active_connection(client)
    job_id = handle.job_id[]
    isnothing(job_id) && return handle
    try
        message = NATS.JetStream.stream_message_get(
            connection,
            RESULT_STREAM,
            result_subject(job_id)
        )
        apply_result!(handle, decode_job_result(NATS.payload(message)))
    catch error
        if missing_durable_result_error(error)
            @debug "No durable result available while reattaching" job_id
        else
            rethrow()
        end
    end
    return handle
end

function transient_durable_result_error(client::BrokerClient, error)
    return client.closed[] || error isa BrokerUnavailable ||
           (error isa NATS.NATSError && error.code in (408, 499))
end

function restore_handles!(client::BrokerClient)
    handles = lock(client.lock) do
        collect(values(client.handles))
    end
    for handle in handles
        handle.state[] in JOB_TERMINAL_STATES && continue
        try
            subscribe_handle!(client, handle)
            fetch_durable_result!(client, handle)
        catch error
            if transient_durable_result_error(client, error)
                @debug "Durable job restoration deferred" job_id=handle.job_id[] exception=error
            else
                client.connection_status[] = :degraded
                @warn "Could not restore durable job state" job_id=handle.job_id[] exception=(
                    error,
                    catch_backtrace()
                )
            end
        end
    end
    return client
end

function monitor_connection!(client::BrokerClient, connection)
    NATS.subscribe(connection, heartbeat_subject()) do message
        try
            register_heartbeat!(client, message)
        catch error
            @warn "Rejected invalid worker heartbeat" exception=(error, catch_backtrace())
        end
    end
    previous = :connecting
    while !client.closed[]
        status = NATS.status(connection)
        if status == NATS.CONNECTED
            client.connection_status[] = :online
            previous == :online || restore_handles!(client)
        elseif status == NATS.CONNECTING
            client.connection_status[] = :connecting
        elseif status == NATS.DISCONNECTED
            client.connection_status[] = :offline
            previous == :offline || mark_pending_offline!(client)
        else
            client.connection_status[] = :degraded
        end
        previous = client.connection_status[]
        prune_workers!(client)
        sleep(0.5)
    end
    return nothing
end

function connection_supervisor!(client::BrokerClient)
    client.connection_status[] = :connecting
    try
        connection = NATS.connect(
            client.url;
            name="linecablemodels-playground-publisher",
            ping_interval=2.0,
            max_pings_out=1,
            retry_on_init_fail=true,
            send_enqueue_when_disconnected=false
        )
        lock(client.lock) do
            client.connection[] = connection
        end
        monitor_connection!(client, connection)
    catch error
        client.closed[] || @warn "NATS publisher connection stopped" exception=(
            error,
            catch_backtrace()
        )
        client.connection_status[] = :offline
    end
    return nothing
end

"""
    start!(client)

Start background broker connection and worker discovery. Returns immediately.
"""
function start!(client::BrokerClient)
    client.closed[] && throw(ArgumentError("cannot restart a closed broker client"))
    lock(client.lock) do
        if isnothing(client.supervisor[]) || istaskdone(client.supervisor[])
            client.supervisor[] = errormonitor(@async connection_supervisor!(client))
        end
    end
    return client
end

function active_connection(client::BrokerClient)
    client.connection_status[] == :online || throw(BrokerUnavailable(
        "Broker is $(client.connection_status[]); calculation was not submitted"
    ))
    connection = client.connection[]
    isnothing(connection) && throw(BrokerUnavailable(
        "Broker connection is unavailable"
    ))
    return connection
end

function unsubscribe_handle!(connection, handle::JobHandle)
    for subscription in handle.subscriptions
        try
            NATS.unsubscribe(connection, subscription)
        catch
        end
    end
    empty!(handle.subscriptions)
    return handle
end

function subscribe_handle!(client::BrokerClient, handle::JobHandle)
    connection = active_connection(client)
    unsubscribe_handle!(connection, handle)
    job_id = handle.job_id[]
    isnothing(job_id) && throw(ArgumentError("job handle has no identifier"))
    event_subscription = NATS.subscribe(
        connection,
        event_subject(job_id)
    ) do message
        try
            event = decode_job_event(NATS.payload(message))
            event.job_id == job_id || return nothing
            apply_event!(handle, event)
            if isterminal(event)
                @async begin
                    sleep(0.05)
                    try
                        fetch_durable_result!(client, handle)
                    catch error
                        if transient_durable_result_error(client, error)
                            @debug "Durable terminal-result load deferred" job_id exception=error
                        else
                            @warn "Could not load durable terminal result" job_id exception=(
                                error,
                                catch_backtrace()
                            )
                        end
                    end
                end
            end
        catch error
            @warn "Rejected invalid job event" job_id exception=(error, catch_backtrace())
        end
        return nothing
    end
    result_subscription = NATS.subscribe(
        connection,
        result_subject(job_id)
    ) do message
        try
            result = decode_job_result(NATS.payload(message))
            result.job_id == job_id || return nothing
            apply_result!(handle, result)
            @async begin
                sleep(0.5)
                try
                    handle.job_id[] == job_id &&
                        unsubscribe_handle!(connection, handle)
                catch
                end
            end
        catch error
            @warn "Rejected invalid job result" job_id exception=(error, catch_backtrace())
        end
        return nothing
    end
    append!(handle.subscriptions, (event_subscription, result_subscription))
    return handle
end

redact_nats_url(url::AbstractString) = replace(
    string(url),
    r"//([^:/@]+):([^@]+)@" => s"//\1:***@"
)

function watch_deadline!(handle::JobHandle, request::JobRequest)
    remaining = Dates.value(
        parse_utc_timestamp(request.deadline) - Dates.now(Dates.UTC)
    ) / 1_000
    remaining > 0 && sleep(remaining)
    expire_handle!(handle, request.job_id)
    return nothing
end

"""
    submit!(client, handle, request)

Subscribe to lifecycle traffic, durably publish `request`, and return after
JetStream accepts it. This method never waits for worker completion.
"""
function submit!(
        client::BrokerClient,
        handle::JobHandle,
        request::JobRequest;
        selector_revision::Union{Nothing,Integer}=nothing
    )
    revision = isnothing(selector_revision) ? lock(handle.lock) do
        handle.input_revision
    end : Int(selector_revision)
    validate(request)
    connection = active_connection(client)
    lock(handle.lock) do
        handle.request = request
        handle.submitted_revision = revision
        handle.job_id[] = request.job_id
        handle.pending_input_hash[] = request.input_hash
        handle.dirty_inputs[] = handle.input_revision != revision
        handle.progress[] = 0.0
        handle.failure[] = nothing
        handle.worker[] = nothing
        handle.last_event_sequence = 0
        set_job_state!(handle, :submitting, "Submitting to JetStream")
    end
    lock(client.lock) do
        client.handles[request.job_id] = handle
    end
    subscribe_handle!(client, handle)
    payload = encode_message(request)
    NATS.JetStream.stream_publish(
        connection,
        job_subject(request.operation, request.priority),
        (payload, ["Nats-Msg-Id" => request.job_id])
    )
    has_worker = lock(client.lock) do
        supports_operation(
            client.workers,
            request.operation;
            engine_constraint=request.engine_constraint
        )
    end
    lock(handle.lock) do
        set_job_state!(
            handle,
            has_worker ? :queued : :unavailable,
            has_worker ? "Queued" : "Queued · no compatible worker online"
        )
    end
    errormonitor(@async watch_deadline!(handle, request))
    return handle
end

"""
    submit_async!(client, handle, request)

Schedule durable submission in a Julia task and return `handle` immediately.
"""
function submit_async!(client::BrokerClient, handle::JobHandle, request::JobRequest)
    revision = lock(handle.lock) do
        handle.input_revision
    end
    @async try
        submit!(client, handle, request; selector_revision=revision)
    catch error
        submission_failed!(handle, error)
    end
    return handle
end

function cancel!(client::BrokerClient, handle::JobHandle)
    job_id = handle.job_id[]
    isnothing(job_id) && return handle
    connection = active_connection(client)
    NATS.publish(connection, cancel_subject(job_id), job_id)
    set_job_state!(handle, :cancelling, "Cancellation requested")
    errormonitor(@async begin
        # Cancellation intentionally uses Core NATS. Repeat while the request is
        # nonterminal so a worker that comes online after the first publication
        # still observes the request.
        while !client.closed[] && handle.job_id[] == job_id &&
                handle.state[] == :cancelling
            sleep(0.5)
            handle.state[] == :cancelling || break
            try
                NATS.publish(active_connection(client), cancel_subject(job_id), job_id)
            catch error
                @debug "Could not repeat cancellation request" job_id exception=error
            end
        end
    end)
    return handle
end

function reattach!(client::BrokerClient, handle::JobHandle, job_id::AbstractString)
    connection = active_connection(client)
    handle.job_id[] = string(job_id)
    lock(client.lock) do
        client.handles[string(job_id)] = handle
    end
    set_job_state!(handle, :queued, "Reattaching to durable job")
    subscribe_handle!(client, handle)
    fetch_durable_result!(client, handle)
    return handle
end

function close!(client::BrokerClient)
    client.closed[] = true
    client.connection_status[] = :closed
    connection = client.connection[]
    if !isnothing(connection)
        try
            NATS.drain(connection)
        catch
        end
    end
    supervisor = client.supervisor[]
    if !isnothing(supervisor) && current_task() !== supervisor
        timedwait(() -> istaskdone(supervisor), 2.0; pollint=0.025)
    end
    return nothing
end
