"""
    JobCoordinator(config, coordinator, science, events)

Own bounded durable job delivery and result reconciliation independently of
heartbeat/lease and preparation traffic. Construction is inert. SQLite retains
submission/cancellation intent; the exact JetStream result remains authoritative.
No scientific code is loaded or evaluated by this owner.
"""
mutable struct JobCoordinator{C<:LeaseCoordinator,S<:ScientificCoordinator}
    "Operator-approved endpoint and worker stream policies."
    config::ControlConfig
    "Existing lease authority, never a second lease registry."
    coordinator::C
    "Existing preparation and idempotent cancellation channel."
    science::S
    "Shared bounded and redacted diagnostics."
    events::ControlEvents
    "Independent job/result connection."
    connection::Union{Nothing,BrokerJobs{CoordinatorIdentity}}
    "At most four finite reconciliation tasks."
    flights::Dict{UUID,Task}
    "Monotonic next-poll times, bounded by pending admitted jobs."
    retry_at::Dict{UUID,Float64}
    "At most four concurrent authorized HTTP result reads; excess requests fail finitely."
    result_readers::Int
    "At most four private artifact reads; joined during shutdown."
    artifact_readers::Int
    "Owned scheduling task."
    task::Union{Nothing,Task}
    "Independent bounded connection attempt."
    connector::Union{Nothing,Task}
    "Next initial connection attempt in monotonic seconds."
    connect_at::Float64
    "Online, offline or stopped transport dimension."
    state::Symbol
    "Reject new work after shutdown begins."
    closed::Bool
    "Serialize brief admission and task ownership changes."
    lock::ReentrantLock
end

function JobCoordinator(config::ControlConfig,coordinator::LeaseCoordinator,science::ScientificCoordinator,events::ControlEvents)
    science.coordinator===coordinator || throw(ArgumentError("Jobs must share preparation's lease authority"))
    JobCoordinator(config,coordinator,science,events,nothing,Dict{UUID,Task}(),Dict{UUID,Float64}(),0,0,
        nothing,nothing,0.0,:offline,false,ReentrantLock())
end
job_store(service::JobCoordinator)=service.coordinator.assignments.inventory.store
job_clock(service::JobCoordinator)=lease_clock(service.coordinator)

function current_job_assignment(service::JobCoordinator,receipt::JobRecord)
    fence=receipt.job.fence
    principal=Principal(fence.owner)
    store=job_store(service)
    get_assignment(store,principal,UUID(fence.lease_id)).fence==fence &&
        get_run(store,principal,UUID(fence.run_id)).state in (:reserved,:starting,:running) &&
        assignment_usable(service.coordinator,principal,UUID(fence.lease_id))
end

"""
    submit_job!(service, principal, lease_id, operation, parameters; request_id=uuid4())

Authorize a registered operation against a live lease and fresh prepared target,
then save the exact server-authored request before any broker publication. An
identical retry returns its original receipt without requiring a new preparation
report or extending its deadline. The HTTP caller never waits for numerical work.
"""
function submit_job!(service::JobCoordinator,principal::Principal,id::UUID,operation::AbstractString,
        parameters::AbstractDict;request_id::UUID=uuid4())
    store=job_store(service)
    previous=prior_job_submission(store,principal,id,operation,parameters,request_id)
    previous===nothing || return previous
    lock(service.lock) do
        previous=prior_job_submission(store,principal,id,operation,parameters,request_id)
        previous===nothing || return previous
        service.closed && throw(AccessDenied(503,"Job service is stopped"))
        service.state==:online && service.connection!==nothing || throw(BrokerUnavailable())
        fence=remote_scientific_fence(service.science,principal,id)
        profile=service.config.profiles.definitions[fence.profile_id]
        operation in profile.operations || throw(AccessDenied(400,"Operation is not permitted by this profile"))
        ready=remote_scientific_status(service.science,principal,id)
        ready.preparation=="ready" || throw(AccessDenied(409,"Fresh explicit preparation is required"))
        target=Protocol.PreparedExecution(ready.executor_id,ready.executor_generation,ready.preparation_key)
        request=Protocol.new_job_request(operation,parameters;session_id=fence.run_id,
            timeout=Millisecond(ceil(Int,1000*profile.budget.job_seconds)))
        job=Protocol.validate(Protocol.AssignedJob("2.0",fence,request,target))
        receipt=lock(service.coordinator.lock) do
            assignment_usable(service.coordinator,principal,id) || throw(AccessDenied(409,"Job assignment changed"))
            reserve_job!(store,principal,job,request_id)
        end
        receipt.job.request.job_id==request.job_id && record_event!(service.events,:job_queued;job)
        return receipt
    end
end

function job_transition!(service::JobCoordinator,receipt::JobRecord,state::Symbol)
    changed=transition_job!(job_store(service),Principal(receipt.job.fence.owner),UUID(receipt.job.request.job_id),state)
    receipt.state==changed.state || record_event!(service.events,Symbol("job_",state);job=receipt.job)
    return changed
end

"""Persist exact owned cancellation intent; do not claim its operation has stopped."""
function cancel_job!(service::JobCoordinator,principal::Principal,id::UUID;request_id::UUID=uuid4())
    lock(service.lock) do
        receipt=get_job(job_store(service),principal,id)
        receipt.state in (:queued,:submitted) || return receipt
        service.closed && throw(AccessDenied(503,"Job service is stopped"))
        current_job_assignment(service,receipt) || throw(AccessDenied(409,"Job assignment is not usable"))
        previous=job_cancellation(job_store(service),principal,id)
        request_job_cancellation!(job_store(service),principal,id,request_id)
        previous===nothing && record_event!(service.events,:job_cancel_requested;job=receipt.job)
        return receipt
    end
end

function reconcile_job_cancellation!(service::JobCoordinator,receipt::JobRecord)
    job=receipt.job
    owner=Principal(job.fence.owner)
    id=UUID(job.request.job_id)
    intent=job_cancellation(job_store(service),owner,id)
    intent===nothing && return true
    intent.acknowledged && return true
    flight=request_scientific!(service.science,owner,UUID(job.fence.lease_id),"cancel_job";
        target_id=job.request.job_id,request_id=intent.request_id)
    command,report=lock(()->(flight.command,flight.report),service.science.lock)
    report===nothing && return false
    command.fence==job.fence && command.action=="cancel_job" && command.target_id==job.request.job_id &&
        command.request_id==string(intent.request_id) && report.accepted || return false
    acknowledge_job_cancellation!(job_store(service),owner,id,intent.request_id)
    return true
end

function reconcile_job!(service::JobCoordinator,id::UUID)
    # Scheduler identities originate only from the private, bounded receipt table.
    receipt=get_job(job_store(service),Principal("runtime-jobs";administrator=true),id)
    receipt.state in (:queued,:submitted) || return nothing
    current_job_assignment(service,receipt) || (job_transition!(service,receipt,:revoked); return nothing)
    service.closed && return nothing
    job=receipt.job
    deadline=Protocol.parse_utc_timestamp(job.request.deadline)
    if now(UTC)>=deadline+Second(5)
        job_transition!(service,receipt,:uncertain)
        return nothing
    end
    connection=lock(()->service.connection,service.lock)
    connection===nothing && throw(BrokerUnavailable())
    receipt.state==:queued && ensure_worker_streams!(connection,service.config.workers[job.fence.worker_id])
    saved=assigned_result(connection,job.fence,job.request.job_id)
    if saved!==nothing
        matching_result(job,saved) || throw(AccessDenied(409,"Stored result differs from its owned receipt"))
        failure=saved.result.failure
        state=failure===nothing ? :succeeded : failure.category=="canceled" ? :canceled :
            failure.category=="execution_uncertain" ? :uncertain : :failed
        job_transition!(service,receipt,state)
        return nothing
    end
    # Cancel-before-publication first records the worker tombstone, then delivers
    # the original job so the worker can persist its canceled terminal result.
    # Cancellation acknowledgement alone is never a fabricated job outcome.
    reconcile_job_cancellation!(service,receipt) || return nothing
    receipt.state==:queued || return nothing
    # Immutable publication retries stay strictly inside the stream's ten-minute
    # duplicate window. Afterwards only inspect for a result; never republish an
    # old job as a fresh message that could repeat an uncertain computation.
    now(UTC)<min(deadline,Protocol.parse_utc_timestamp(job.request.submitted_at)+Minute(2)) || return nothing
    service.closed && return nothing
    acknowledgment=publish_assigned_job!(connection,service.coordinator,Principal(job.fence.owner),job)
    acknowledgment.stream==job_stream(job.fence.worker_id) || throw(AccessDenied(409,"Job was stored in an unexpected stream"))
    job_transition!(service,receipt,:submitted)
    return nothing
end

"""
    owned_job_result(service, principal, job_id)

Read a result only through an owned receipt and exact request/target provenance.
Historical data may remain readable after lease loss; it does not regain current
execution authority. Return nothing only for an actual broker not-found response,
not when the transport is unavailable.
"""
function owned_job_result(service::JobCoordinator,principal::Principal,id::UUID)
    receipt=get_job(job_store(service),principal,id)
    connection=lock(service.lock) do
        !service.closed && service.connection!==nothing || throw(BrokerUnavailable())
        service.result_readers<4 || throw(AccessDenied(429,"Job result read capacity is occupied"))
        service.result_readers+=1
        service.connection
    end
    try
        result=assigned_result(connection,receipt.job.fence,string(id))
        result===nothing || matching_result(receipt.job,result) || throw(AccessDenied(409,"Job result provenance differs"))
        if result!==nothing && receipt.state==:uncertain
            failure=result.result.failure
            state=failure===nothing ? :succeeded : failure.category=="canceled" ? :canceled :
                failure.category=="execution_uncertain" ? :uncertain : :failed
            job_transition!(service,receipt,state)
        end
        return result
    finally
        lock(service.lock) do; service.result_readers-=1; end
    end
end

function tick_jobs!(service::JobCoordinator)
    store=job_store(service)
    ids=lock(store.lock) do
        [UUID(row.job_id) for row in sql_rows(store.db,
            "SELECT job_id FROM jobs WHERE state IN $PENDING_JOB_SQL ORDER BY created_at,job_id LIMIT ?",
            (service.coordinator.assignments.limits.total,))]
    end
    lock(service.lock) do
        service.closed && return nothing
        for (id,task) in collect(service.flights)
            istaskdone(task) && delete!(service.flights,id)
        end
        for id in collect(keys(service.retry_at))
            id in ids || delete!(service.retry_at,id)
        end
        connection=service.connection
        service.state=connection!==nothing && !connection.closed && NATS.status(connection.connection)==NATS.CONNECTED ? :online : :offline
        # Old jobs must not repeatedly occupy all slots during slow transport.
        # Least-recently attempted work goes first, with stable receipt ordering
        # for ties; the next delay starts after the finite attempt completes.
        sort!(ids;by=id->get(service.retry_at,id,-Inf),alg=MergeSort)
        for id in ids
            length(service.flights)<4 || break
            haskey(service.flights,id) && continue
            job_clock(service)>=get(service.retry_at,id,0.0) || continue
            service.retry_at[id]=job_clock(service)+1
            service.flights[id]=@async try
                reconcile_job!(service,id)
            catch
                # Keep the saved request intact. Neither a network failure nor
                # malformed stored result is permission to invent/replay work.
                service.state=:offline
            finally
                lock(service.lock) do
                    service.retry_at[id]=job_clock(service)+1
                end
            end
        end
    end
    return nothing
end

"""Retrieve a private artifact only through an owned receipt and exact durable result."""
function owned_job_artifact(service::JobCoordinator,principal::Principal,id::UUID)
    receipt=get_job(job_store(service),principal,id)
    location=lock(service.lock) do
        !service.closed && service.config.artifacts!==nothing || throw(ArtifactUnavailable())
        service.artifact_readers<4 || throw(AccessDenied(429,"Artifact read capacity is occupied"))
        service.artifact_readers+=1
        service.config.artifacts
    end
    try
        try
            outcome=owned_job_result(service,principal,id)
            outcome!==nothing && outcome.result.artifact!==nothing || throw(AccessDenied(404,"Job artifact not found"))
            bytes=read_job_artifact(location,receipt.job,outcome.result.artifact)
            bytes===nothing && throw(ArtifactUnavailable())
            return bytes
        catch error
            error isa AccessDenied && rethrow()
        end
        throw(ArtifactUnavailable())
    finally
        lock(service.lock) do;service.artifact_readers-=1;end
    end
end

function connect_jobs!(service::JobCoordinator)
    lock(service.lock) do
        service.closed && return
        service.connection===nothing || return
        service.connector!==nothing && !istaskdone(service.connector) && return
        job_clock(service)>=service.connect_at || return
        service.connect_at=job_clock(service)+2
        service.connector=@async begin
            connection=nothing
            try
                connection=BrokerJobs(service.config.endpoint,CoordinatorIdentity())
                # Provision/verify server-owned transport before announcing the
                # job connection. First-use stream/JSON compilation must not be
                # deferred into the first live lease's publication deadline.
                # This neither allocates workers nor imports numerical code.
                for trust in values(service.config.workers)
                    any(id->service.config.profiles.definitions[id].kind==:scientific,trust.profiles) || continue
                    service.closed && break
                    ensure_worker_streams!(connection,trust)
                end
                lock(service.lock) do
                    if !service.closed; service.connection=connection; connection=nothing; end
                end
            catch
                service.state=:offline
            finally
                connection===nothing || close(connection)
            end
        end
    end
end

"""Start the independent bounded job scheduler without waiting for broker availability."""
function start_jobs!(service::JobCoordinator)
    lock(service.lock) do
        service.closed && throw(ArgumentError("Job coordinator is closed"))
        service.task===nothing || return service
        service.task=@async while !service.closed
            try connect_jobs!(service); tick_jobs!(service) catch; service.state=:offline end
            sleep(0.05)
        end
    end
    return service
end

function Base.close(service::JobCoordinator)
    lock(service.lock) do; service.closed=true; end
    service.task===nothing || wait(service.task)
    service.connector===nothing || wait(service.connector)
    tasks=lock(()->collect(values(service.flights)),service.lock)
    drained()=all(istaskdone,tasks) && lock(()->service.result_readers==0 && service.artifact_readers==0,service.lock)
    timedwait(drained,30;pollint=0.025)==:ok || throw(ArgumentError("Job reconciliation cleanup remains unresolved"))
    foreach(wait,tasks)
    connection=lock(service.lock) do
        previous=service.connection; service.connection=nothing
        empty!(service.flights); empty!(service.retry_at); service.state=:stopped
        previous
    end
    connection===nothing || close(connection)
    return nothing
end
