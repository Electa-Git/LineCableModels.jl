"""
    BrokerJobs(endpoint, identity)

Own a scientific job/result connection separate from worker control and terminal
traffic. This transport does not load an engine or establish executor readiness.
"""
mutable struct BrokerJobs{I<:AbstractBrokerIdentity}
    "Server-provisioned role."
    identity::I
    "Dedicated bounded scientific transport connection."
    connection::NATS.Connection
    "Idempotent teardown flag."
    closed::Bool
    "Serialize teardown with local transmission."
    lock::ReentrantLock
end
BrokerJobs(endpoint::BrokerEndpoint, identity::AbstractBrokerIdentity) =
    BrokerJobs(identity, connect_broker(endpoint, identity), false, ReentrantLock())

function matching_result(job::Protocol.AssignedJob,outcome::Protocol.AssignedResult)
    outcome.fence==job.fence && outcome.execution==job.execution &&
        outcome.result.job_id==job.request.job_id && outcome.result.input_hash==job.request.input_hash &&
        outcome.result.operation==job.request.operation
end

function Base.close(transport::BrokerJobs)
    lock(transport.lock) do
        transport.closed && return nothing
        transport.closed = true
        Logging.with_logger(Logging.NullLogger()) do
            NATS.drain(transport.connection)
        end
    end
    return nothing
end

function job_request(transport::BrokerJobs, ::Type{T}, subject::String, data=nothing; timeout=2.0) where T
    transport.closed && throw(BrokerUnavailable())
    try
        return NATS.request(T, transport.connection, subject, data; timeout)
    catch error
        error isa NATS.NATSError && rethrow()
        throw(BrokerUnavailable())
    end
end

function worker_stream_config(worker_id::String, kind::Symbol)
    kind in (:jobs, :results) || throw(ArgumentError("invalid runtime stream kind"))
    jobs = kind == :jobs
    name = jobs ? job_stream(worker_id) : result_stream(worker_id)
    subject = jobs ? "lcm.jobs.v2.$worker_id.>" : "lcm.results.v2.$worker_id.>"
    return NATS.JetStream.StreamConfiguration(name=name, subjects=[subject],
        retention=jobs ? :workqueue : :limits, storage=:file, num_replicas=1,
        max_consumers=jobs ? 1 : -1, max_msgs=jobs ? 256 : 4096,
        max_msgs_per_subject=jobs ? -1 : 1, max_bytes=64 * 1024^2,
        max_age=86_400_000_000_000, max_msg_size=Int32(262144), discard=:new,
        duplicate_window=600_000_000_000, allow_direct=!jobs,
        deny_delete=true, deny_purge=true, allow_rollup_hdrs=false)
end

function compatible_stream(found, expected)
    fields = (:name, :subjects, :retention, :storage, :num_replicas, :max_consumers,
        :max_msgs, :max_msgs_per_subject, :max_bytes, :max_age, :max_msg_size,
        :discard, :duplicate_window, :allow_direct, :deny_delete, :deny_purge, :allow_rollup_hdrs)
    booleans = (:allow_direct, :deny_delete, :deny_purge, :allow_rollup_hdrs)
    equal(field) = field in booleans ?
        something(getfield(found, field), false) == something(getfield(expected, field), false) :
        getfield(found, field) == getfield(expected, field)
    return all(equal, fields) &&
        found.mirror === nothing && found.sources === nothing && found.republish === nothing &&
        found.subject_transform === nothing
end

"""
    ensure_worker_streams!(transport, trust)

Create only the registered worker's bounded v2 job/result streams and fixed
consumer. Existing configuration must match; no stream is silently resized,
purged, remapped or replaced. Legacy v1 streams are untouched.
"""
function ensure_worker_streams!(transport::BrokerJobs{CoordinatorIdentity}, trust::WorkerTrust)
    for kind in (:jobs, :results)
        config = worker_stream_config(trust.worker_id, kind)
        found = job_request(transport, Union{NATS.JetStream.StreamInfo,NATS.JetStream.ApiError},
            "\$JS.API.STREAM.INFO.$(config.name)")
        if found isa NATS.JetStream.ApiError && found.code == 404
            found = job_request(transport, Union{NATS.JetStream.StreamInfo,NATS.JetStream.ApiError},
                "\$JS.API.STREAM.CREATE.$(config.name)", JSON3.write(config))
        end
        found isa NATS.JetStream.StreamInfo || throw(AccessDenied(409, "Runtime stream provisioning failed"))
        compatible_stream(found.config, config) ||
            throw(AccessDenied(409, "Existing runtime stream configuration requires operator reconciliation"))
    end
    stream = job_stream(trust.worker_id)
    config = NATS.JetStream.ConsumerConfiguration(name="agent", durable_name="agent",
        ack_policy=:explicit, ack_wait=30_000_000_000, max_deliver=3,
        filter_subjects=["lcm.jobs.v2.$(trust.worker_id).>"], max_ack_pending=trust.capacity,
        max_waiting=trust.capacity, max_batch=1, max_expires=1_000_000_000, max_bytes=262144)
    found = job_request(transport, Union{NATS.JetStream.ConsumerInfo,NATS.JetStream.ApiError},
        "\$JS.API.CONSUMER.INFO.$stream.agent")
    if found isa NATS.JetStream.ApiError && found.code == 404
        found = job_request(transport, Union{NATS.JetStream.ConsumerInfo,NATS.JetStream.ApiError},
            "\$JS.API.CONSUMER.CREATE.$stream.agent", JSON3.write((stream_name=stream, config=config)))
    end
    found isa NATS.JetStream.ConsumerInfo || throw(AccessDenied(409, "Runtime consumer provisioning failed"))
    for field in (:name, :durable_name, :ack_policy, :ack_wait, :max_deliver, :filter_subjects,
            :max_ack_pending, :max_waiting, :max_batch, :max_expires, :max_bytes, :deliver_subject)
        getfield(found.config, field) == getfield(config, field) ||
            throw(AccessDenied(409, "Existing runtime consumer configuration requires operator reconciliation"))
    end
    return nothing
end

"""
    publish_assigned_job!(transport, coordinator, principal, job)

Authorize the exact acknowledged assignment and allowlisted scientific operation
before durable targeted publication. The caller preserves job_id for retries.
Preparation/queue admission belong to the fixed execution orchestrator and must
precede using this low-level transport; this function does not prepare a model.
"""
function publish_assigned_job!(transport::BrokerJobs{CoordinatorIdentity}, coordinator::LeaseCoordinator,
        principal::Principal, job::Protocol.AssignedJob)
    Protocol.validate(job)
    id = UUID(job.fence.lease_id)
    lease = get_assignment(coordinator.assignments.inventory.store, principal, id)
    lease.fence == job.fence && assignment_usable(coordinator, principal, id) ||
        throw(AccessDenied(409, "Assignment is not usable"))
    profile = coordinator.assignments.inventory.profiles.definitions[job.fence.profile_id]
    profile.kind == :scientific && job.request.operation in profile.operations ||
        throw(AccessDenied(400, "Operation is not permitted for this profile"))
    Protocol.parse_utc_timestamp(job.request.deadline) > now(UTC) ||
        throw(AccessDenied(409, "Scientific request deadline expired"))
    return job_request(transport, NATS.JetStream.PubAck, Protocol.assigned_job_subject(job.fence),
        (Protocol.encode_message(job), ["Nats-Msg-Id"=>Protocol.assigned_job_subject(job.fence) * "." * job.request.job_id]))
end

"""Retain an exact pulled message together with its validated assigned request."""
struct AssignedDelivery
    "Strictly validated passive request."
    job::Protocol.AssignedJob
    "Broker delivery used only for acknowledgement."
    message::NATS.Msg
end

function owned_delivery(transport::BrokerJobs{WorkerIdentity}, message::NATS.Msg)
    reply = something(message.reply_to, "")
    prefix = "\$JS.ACK.$(job_stream(transport.identity.worker_id)).agent."
    startswith(reply, prefix) || throw(AccessDenied(403, "Unexpected job acknowledgement subject"))
    return reply
end

function validate_assigned_delivery(transport::BrokerJobs{WorkerIdentity}, delivery::AssignedDelivery)
    job = Protocol.validate(delivery.job)
    reply = owned_delivery(transport, delivery.message)
    transport.identity.worker_id == job.fence.worker_id &&
        delivery.message.subject == Protocol.assigned_job_subject(job.fence) ||
        throw(AccessDenied(409, "Delivery does not match the assigned worker"))
    decoded = try
        Protocol.decode_runtime_message(Protocol.AssignedJob, NATS.payload(delivery.message))
    catch
        throw(AccessDenied(400, "Invalid scientific delivery payload"))
    end
    decoded == job || throw(AccessDenied(409, "Delivery does not match the submitted scientific request"))
    delivery_count(transport, delivery)
    return reply
end

"""
    progress_assigned_delivery!(transport, ledger, delivery)

Extend only this input's acknowledgement timer while its lease remains live.
This is work-in-progress, not successful completion; result persistence must
still precede the terminal acknowledgement.
"""
function progress_assigned_delivery!(transport::BrokerJobs{WorkerIdentity},ledger::AgentLeaseLedger,delivery::AssignedDelivery)
    agent_lease_usable(ledger,delivery.job.fence) || throw(AccessDenied(409,"Job lease is no longer usable"))
    job_request(transport,NATS.Msg,validate_assigned_delivery(transport,delivery),"+WPI";timeout=0.5)
    return nothing
end

function delivery_count(transport::BrokerJobs{WorkerIdentity},delivery::AssignedDelivery)
    reply=owned_delivery(transport,delivery.message)
    fields=split(reply,'.')
    length(fields)==9 || throw(AccessDenied(403,"Unexpected job delivery metadata"))
    values=tryparse.(UInt64,fields[5:9])
    all(value->value!==nothing,values) && all(value->value>0,values[1:4]) ||
        throw(AccessDenied(403,"Invalid job delivery metadata"))
    return values[1]
end

function terminate_assigned_delivery!(transport::BrokerJobs{WorkerIdentity},delivery::AssignedDelivery)
    job_request(transport,NATS.Msg,owned_delivery(transport,delivery.message),"+TERM";timeout=0.5)
    return nothing
end

"""
    poll_assigned_job!(transport, ledger) -> Union{Nothing,AssignedDelivery}

Perform one finite pull, not an indefinite wait loop. Reject and terminate stale
or foreign authority before returning a request to the executor orchestrator.
No operation runs in the transport. At most one stale request is consumed per
call, leaving control scheduling independent.
"""
function poll_assigned_job!(transport::BrokerJobs{WorkerIdentity}, ledger::AgentLeaseLedger)
    transport.identity.worker_id == ledger.worker_id || throw(AccessDenied(403, "Agent identity mismatch"))
    stream = job_stream(ledger.worker_id)
    message = try
        job_request(transport, NATS.Msg, "\$JS.API.CONSUMER.MSG.NEXT.$stream.agent",
            JSON3.write((batch=1, expires=100_000_000)); timeout=0.3)
    catch error
        error isa NATS.NATSError && error.code in (404, 408) && return nothing
        rethrow()
    end
    reply = owned_delivery(transport, message)
    job = try
        Protocol.decode_runtime_message(Protocol.AssignedJob, NATS.payload(message))
    catch
        job_request(transport, NATS.Msg, reply, "+TERM")
        return nothing
    end
    permitted = message.subject == Protocol.assigned_job_subject(job.fence) &&
        agent_lease_usable(ledger, job.fence)
    if permitted
        profile = ledger.profiles.definitions[job.fence.profile_id]
        permitted = profile.kind == :scientific && job.request.operation in profile.operations &&
            Protocol.parse_utc_timestamp(job.request.deadline) > now(UTC)
    end
    if !permitted
        job_request(transport, NATS.Msg, reply, "+TERM")
        return nothing
    end
    return AssignedDelivery(job, message)
end

"""
    assigned_result(transport, fence, job_id) -> Union{Nothing,AssignedResult}

Read an exact result from its worker-scoped stream. Worker transports may read
only their own identity; public callers must authorize the owned job separately
before invoking the coordinator transport.
"""
function assigned_result(transport::BrokerJobs, fence::Protocol.AssignmentFence, job_id::AbstractString)
    transport.identity isa WorkerIdentity && transport.identity.worker_id != fence.worker_id &&
        throw(AccessDenied(403, "Result worker identity mismatch"))
    subject = Protocol.assigned_result_subject(fence, job_id)
    message = try
        job_request(transport, NATS.Msg, "\$JS.API.DIRECT.GET.$(result_stream(fence.worker_id))",
            JSON3.write((last_by_subj=subject,)))
    catch error
        error isa NATS.NATSError && error.code == 404 && return nothing
        rethrow()
    end
    outcome = Protocol.decode_runtime_message(Protocol.AssignedResult, NATS.payload(message))
    outcome.fence == fence && outcome.result.job_id == job_id ||
        throw(AccessDenied(409, "Stored result does not match its assignment"))
    return outcome
end

"""
    persist_assigned_result!(transport, ledger, delivery, outcome)

Validate provenance and persist a terminal result before acknowledging its input.
A result acknowledgement failure leaves the job unacknowledged for bounded
redelivery. Callers check assigned_result before re-executing a redelivered job;
the transport never promises exactly-once effects.
"""
function persist_assigned_result!(transport::BrokerJobs{WorkerIdentity}, ledger::AgentLeaseLedger,
        delivery::AssignedDelivery, outcome::Protocol.AssignedResult)
    Protocol.validate(outcome)
    job = delivery.job
    outcome.fence == job.fence && outcome.result.job_id == job.request.job_id &&
        outcome.execution == job.execution &&
        outcome.result.operation == job.request.operation && outcome.result.input_hash == job.request.input_hash ||
        throw(AccessDenied(409, "Result does not match the submitted scientific request"))
    transport.identity.worker_id == outcome.fence.worker_id &&
        agent_lease_usable(ledger, outcome.fence) || throw(AccessDenied(409, "Result assignment is no longer usable"))
    reply = validate_assigned_delivery(transport, delivery)
    subject = Protocol.assigned_result_subject(outcome.fence, outcome.result.job_id)
    saved = assigned_result(transport, outcome.fence, outcome.result.job_id)
    if saved === nothing
        acknowledgement = job_request(transport, NATS.JetStream.PubAck, subject,
            (Protocol.encode_message(outcome), ["Nats-Msg-Id"=>subject]))
        saved = something(acknowledgement.duplicate, false) ?
            assigned_result(transport, outcome.fence, outcome.result.job_id) : outcome
    end
    saved !== nothing && saved.execution == job.execution && saved.result.input_hash == job.request.input_hash &&
        saved.result.operation == job.request.operation ||
        throw(AccessDenied(409, "Durable result does not match its original request"))
    job_request(transport, NATS.Msg, reply, "+ACK")
    return saved
end
