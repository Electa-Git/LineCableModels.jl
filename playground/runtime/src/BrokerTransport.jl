import NATS, Logging

"""Report a bounded transport failure without exposing broker credentials."""
struct BrokerUnavailable <: Exception end
Base.showerror(io::IO, ::BrokerUnavailable) = print(io, "Runtime broker is unavailable")

"""
    connect_broker(endpoint, identity) -> NATS.Connection

Open an authenticated control connection with finite handshakes, private replies
and no disconnected publication replay. Broker-library diagnostics are suppressed
inside the owned connection tasks because they may include raw authentication or
payload data; the runtime exposes connection state and bounded reason codes.
"""
function connect_broker(endpoint::BrokerEndpoint, identity::AbstractBrokerIdentity)
    broker_file(endpoint.password_file; private=true, max_bytes=4096)
    password = chomp(read(endpoint.password_file, String))
    !isempty(password) && !any(c -> c in ('\r', '\n', '\0'), password) ||
        throw(ArgumentError("invalid broker password file"))
    for (path, private) in ((endpoint.ca_file, false), (endpoint.certificate_file, false), (endpoint.key_file, true))
        path === nothing || broker_file(path; private)
    end
    try
        return Logging.with_logger(Logging.NullLogger()) do
            NATS.connect(endpoint.url; user=broker_user(identity), pass=password,
                auth_token=nothing, jwt=nothing, nkey=nothing, nkey_seed=nothing,
                verbose=false, pedantic=true, name=broker_user(identity),
                tls_required=startswith(endpoint.url, "tls://"),
                tls_ca_path=endpoint.ca_file, tls_cert_path=endpoint.certificate_file,
                tls_key_path=endpoint.key_file, tls_server_name=endpoint.server_name,
                connect_timeout=2.0, retry_on_init_fail=false,
                ignore_advertised_servers=true, retain_servers_order=true,
                inbox_prefix=broker_inbox(identity), send_enqueue_when_disconnected=false,
                send_buffer_limit=262144, send_retry_delays=Float64[],
                ping_interval=1.0, max_pings_out=2, drain_timeout=0.5, drain_poll=0.01,
                reconnect_delays=Base.ExponentialBackOff(n=typemax(Int),
                    first_delay=0.1, max_delay=1.0, jitter=0.1))
        end
    catch
        throw(BrokerUnavailable())
    end
end

"""Carry a validated record with its broker-enforced subject identity."""
struct ControlEnvelope{T<:Protocol.RuntimeRecord}
    "Identity extracted from the exact subscribed subject, not the payload."
    worker_id::String
    "Strictly decoded passive record."
    record::T
end

"""
    BrokerControl(endpoint, identity; worker_ids=(), traffic=:control)

Own only the selected record subscriptions. Per-worker queues and round-robin polling keep
one report stream from monopolizing other workers. Scientific result/log/terminal
traffic uses separate connections; no user callback runs in a NATS handler.
Use a separate instance with traffic=:science for preparation/cancellation/status
records; these never share the heartbeat/lease connection's input queues.
"""
mutable struct BrokerControl{I<:AbstractBrokerIdentity}
    "Broker-authorized role."
    identity::I
    "Dedicated control connection."
    connection::NATS.Connection
    "Exact subject subscriptions, one queue for each worker and record type."
    subscriptions::Vector{Tuple{String,DataType,NATS.Sub}}
    "Next queue to inspect."
    cursor::Int
    "Number of rejected records, saturated at the largest Int."
    rejected::Int
    "Idempotent teardown flag."
    closed::Bool
    "Serialize polling and teardown."
    lock::ReentrantLock
end

control_inputs(::CoordinatorIdentity, ids) = [(id, kind, subject)
    for id in ids for (kind, subject) in
        ((Protocol.WorkerAnnouncement, "lcm.report.v2.$id"),
         (Protocol.LeaseAcknowledgement, "lcm.ack.v2.$id"))]
control_inputs(identity::WorkerIdentity, ids) =
    [(identity.worker_id, Protocol.WorkerProbe, "lcm.control.v2.$(identity.worker_id).probe"),
     (identity.worker_id, Protocol.LeaseControl, "lcm.control.v2.$(identity.worker_id).lease")]

science_inputs(::CoordinatorIdentity, ids) = [(id,Protocol.ScientificReport,"lcm.science.v2.$id.report") for id in ids]
science_inputs(identity::WorkerIdentity, ids) =
    [(identity.worker_id,Protocol.ScientificCommand,"lcm.science.v2.$(identity.worker_id).command")]

function BrokerControl(endpoint::BrokerEndpoint, identity::AbstractBrokerIdentity; worker_ids=(),traffic::Symbol=:control)
    traffic in (:control,:science) || throw(ArgumentError("unsupported runtime record channel"))
    ids = sort!(unique(Protocol.runtime_token.(collect(worker_ids))))
    identity isa CoordinatorIdentity && !(1 <= length(ids) <= 128) &&
        throw(ArgumentError("control connection requires 1:128 provisioned workers"))
    connection = connect_broker(endpoint, identity)
    control = BrokerControl(identity, connection, Tuple{String,DataType,NATS.Sub}[],
        1, 0, false, ReentrantLock())
    try
        Logging.with_logger(Logging.NullLogger()) do
            for (id, kind, subject) in (traffic==:control ? control_inputs(identity,ids) : science_inputs(identity,ids))
                subscription = NATS.subscribe(connection, subject; channel_size=4)
                push!(control.subscriptions, (id, kind, subscription))
            end
            NATS.ping(connection; timeout=2, measure=false)
        end
        return control
    catch
        close(control)
        throw(BrokerUnavailable())
    end
end

control_subject(::CoordinatorIdentity, record::Protocol.WorkerProbe) =
    "lcm.control.v2.$(record.worker_id).probe"
control_subject(::CoordinatorIdentity, record::Protocol.LeaseControl) =
    "lcm.control.v2.$(record.fence.worker_id).lease"
control_subject(::CoordinatorIdentity, record::Protocol.ScientificCommand) =
    "lcm.science.v2.$(record.fence.worker_id).command"
function control_subject(identity::WorkerIdentity, record::Protocol.ScientificReport)
    identity.worker_id == record.fence.worker_id || throw(AccessDenied(403,"Worker identity mismatch"))
    return "lcm.science.v2.$(identity.worker_id).report"
end
function control_subject(identity::WorkerIdentity, record::Protocol.WorkerAnnouncement)
    identity.worker_id == record.worker_id || throw(AccessDenied(403, "Worker identity mismatch"))
    return "lcm.report.v2.$(identity.worker_id)"
end
function control_subject(identity::WorkerIdentity, record::Protocol.LeaseAcknowledgement)
    identity.worker_id == record.fence.worker_id || throw(AccessDenied(403, "Worker identity mismatch"))
    return "lcm.ack.v2.$(identity.worker_id)"
end
control_subject(::AbstractBrokerIdentity, ::Protocol.RuntimeRecord) =
    throw(AccessDenied(403, "Record is not permitted for this broker role"))

"""
    send_control!(control, record)

Validate and publish one transient, role-appropriate control record. Returning
does not mean a grant was acknowledged. Request IDs, revisions and local lease
timers decide whether the receiving party may act.
"""
function send_control!(control::BrokerControl, record::Protocol.RuntimeRecord)
    Protocol.validate(record)
    subject = control_subject(control.identity, record)
    lock(control.lock) do
        !control.closed && NATS.status(control.connection) == NATS.CONNECTED ||
            throw(BrokerUnavailable())
        try
            NATS.publish(control.connection, subject, JSON3.write(record))
        catch
            throw(BrokerUnavailable())
        end
    end
    return nothing
end

record_worker(record::Union{Protocol.WorkerProbe,Protocol.WorkerAnnouncement}) = record.worker_id
record_worker(record::Union{Protocol.LeaseControl,Protocol.LeaseAcknowledgement}) = record.fence.worker_id
record_worker(record::Union{Protocol.ScientificCommand,Protocol.ScientificReport}) = record.fence.worker_id

"""
    poll_control!(control; limit=32) -> Vector{ControlEnvelope}

Read at most limit messages without waiting. Validate exact subject, record
shape and worker identity. Count malformed frames without logging their content.
"""
function poll_control!(control::BrokerControl; limit::Integer=32)
    !(limit isa Bool) && 1 <= limit <= 128 || throw(ArgumentError("invalid control poll bound"))
    return lock(control.lock) do
        result = ControlEnvelope[]
        control.closed && return result
        queues = length(control.subscriptions)
        empty_checks = 0
        consumed = 0
        while empty_checks < queues && consumed < limit
            index = control.cursor
            control.cursor = mod1(index + 1, queues)
            id, kind, sub = control.subscriptions[index]
            message = NATS.next(control.connection, sub; no_wait=true, no_throw=true)
            if message === nothing
                empty_checks += 1
                continue
            end
            empty_checks = 0
            consumed += 1
            try
                message.subject == sub.subject && message.reply_to === nothing &&
                    length(message.payload) <= 262144 || throw(ArgumentError("invalid control frame"))
                record = Protocol.decode_runtime_message(kind, NATS.payload(message))
                record_worker(record) == id || throw(ArgumentError("worker subject mismatch"))
                push!(result, ControlEnvelope(id, record))
            catch
                control.rejected == typemax(Int) || (control.rejected += 1)
            end
        end
        return result
    end
end

function Base.close(control::BrokerControl)
    lock(control.lock) do
        control.closed && return nothing
        control.closed = true
        Logging.with_logger(Logging.NullLogger()) do
            for (_, _, sub) in control.subscriptions
                try
                    NATS.unsubscribe(control.connection, sub)
                catch
                    # The owned connection drain also clears disconnected queues.
                end
            end
            NATS.drain(control.connection)
        end
        empty!(control.subscriptions)
    end
    return nothing
end
