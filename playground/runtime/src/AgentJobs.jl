"""Own one delivered input until its bounded terminal persistence/acknowledgement attempt ends."""
mutable struct AgentJobFlight
    "Exact broker delivery and prepared-execution target."
    delivery::AssignedDelivery
    "Owned execution/persistence task."
    task::Union{Nothing,Task}
    "Last bounded public stage; never arbitrary exception text."
    phase::Symbol
    "Completed result retained while retrying only persistence, not computation."
    outcome::Union{Nothing,Protocol.AssignedResult}
end

"""
    AgentJobService(endpoint, resources, ledger)

Consume only this agent's provisioned v2 job stream through an independent
connection. Construction is inert. The existing ScientificResources owner checks
lease, preparation, executor generation, deadlines and physical policy. Delivery
retries inspect a durable result first; uncertain computation is never rerun.
"""
mutable struct AgentJobService{R<:ScientificResources,L<:AgentLeaseLedger}
    "Operator-owned endpoint, never a browser broker credential."
    endpoint::BrokerEndpoint
    "The agent's existing scientific resource owner."
    resources::R
    "The same worker-incarnation lease authority as the control scheduler."
    ledger::L
    "Private result writer, outside every scientific child."
    artifacts::Union{Nothing,AbstractRuntimeArtifacts}
    "Separate job/result connection."
    connection::Union{Nothing,BrokerJobs{WorkerIdentity}}
    "At most capacity delivered jobs, including results awaiting persistence."
    flights::Dict{String,AgentJobFlight}
    "Owned bounded-pull scheduler."
    task::Union{Nothing,Task}
    "Independent finite initial connection attempt."
    connector::Union{Nothing,Task}
    "Next initial connection attempt in local monotonic seconds."
    retry_at::Float64
    "Online/offline/stopped transport dimension, not preparation readiness."
    state::Symbol
    "Reject new deliveries after shutdown begins."
    closed::Bool
    "Serialize short ownership updates, never numerical work."
    lock::ReentrantLock
end

function AgentJobService(endpoint::BrokerEndpoint,resources::ScientificResources,ledger::AgentLeaseLedger;
        artifacts::Union{Nothing,AbstractRuntimeArtifacts}=nothing)
    resources.ledger===ledger || throw(ArgumentError("job service must share the agent lease authority"))
    AgentJobService(endpoint,resources,ledger,artifacts,nothing,Dict{String,AgentJobFlight}(),
        nothing,nothing,0.0,:offline,false,ReentrantLock())
end

function job_outcome(job::Protocol.AssignedJob,started::String,output::ScientificOutput;artifact=nothing)
    result=Protocol.JobResult("1.0",job.request.job_id,job.request.operation,output.schema_version,
        job.request.input_hash,"environment-sha256:"*job.fence.fingerprint,job.fence.fingerprint,
        job.fence.worker_id,"bypass",started,Protocol.utc_timestamp(),artifact===nothing ? output.value : nothing,artifact,nothing,output.warnings)
    Protocol.validate(Protocol.AssignedResult("2.0",job.fence,result,job.execution))
end

function job_failure(job::Protocol.AssignedJob,started::String,code::String)
    code in ("not_admitted","operation_rejected","canceled","deadline","executor_failed",
        "execution_uncertain","result_payload_limit","artifact_unavailable") || throw(ArgumentError("Unknown scientific failure code"))
    failure=Protocol.FailureInfo(code,"Scientific job failed: $code",string(uuid4()),false)
    result=Protocol.JobResult("1.0",job.request.job_id,job.request.operation,"runtime.failure.v1",
        job.request.input_hash,"environment-sha256:"*job.fence.fingerprint,job.fence.fingerprint,
        job.fence.worker_id,"bypass",started,Protocol.utc_timestamp(),nothing,nothing,failure,String[])
    Protocol.validate(Protocol.AssignedResult("2.0",job.fence,result,job.execution))
end

function admit_delivered_job!(service::AgentJobService,flight::AgentJobFlight)
    job=flight.delivery.job
    until=time_ns()/1e9+5
    while true
        service.closed && throw(AccessDenied(409,"Job service is closed"))
        pending,admitted=lock(service.resources.lock) do
            check_job_cancellation(service.resources,job.fence,job.request.job_id)
            prepared_execution(service.resources,job.fence)==job.execution ||
                throw(AccessDenied(409,"Prepared executor changed before admission"))
            handle=service.resources.handles[job.fence.lease_id]
            # A read-only readiness probe may be in flight when its prior report
            # arrives at the client. It must not turn a valid job into a busy
            # rejection. Do not wait behind preparation or another calculation.
            if handle.request_kind==:inspect && handle.task!==nothing && !istaskdone(handle.task)
                return handle.task,false
            end
            return execute_assigned!(service.resources,job),true
        end
        admitted && return pending
        flight.phase=:admitting
        while !istaskdone(pending)
            !service.closed && agent_lease_usable(service.ledger,job.fence) &&
                time_ns()/1e9<until && now(UTC)<Protocol.parse_utc_timestamp(job.request.deadline) ||
                throw(AccessDenied(409,"Scientific inspection did not release job admission"))
            sleep(0.025)
        end
        time_ns()/1e9<until || throw(AccessDenied(409,"Scientific inspection admission timed out"))
        # Revalidate the prepared identity and reserve execution under one lock.
    end
end

function execute_delivered_job!(service::AgentJobService,flight::AgentJobFlight,connection)
    job=flight.delivery.job
    started=Protocol.utc_timestamp()
    numerical=nothing
    try
        saved=assigned_result(connection,job.fence,job.request.job_id)
        if saved!==nothing
            matching_result(job,saved) || throw(AccessDenied(409,"Stored job result differs from its request"))
            return saved
        end
        # Losing the process between computation and persistence leaves no proof
        # whether it ran. Do not turn broker redelivery into a new computation.
        delivery_count(connection,flight.delivery)>1 && return job_failure(job,started,"execution_uncertain")
        numerical=admit_delivered_job!(service,flight)
        flight.phase=:running
        output=fetch(numerical)
        bytes=collect(codeunits(JSON3.write(output.value)))
        artifact=nothing
        if length(bytes)>ASSIGNED_INLINE_BYTES && service.artifacts!==nothing
            flight.phase=:storing
            try
                artifact=store_job_artifact!(service.artifacts,job,bytes)
            catch
                return job_failure(job,started,"artifact_unavailable")
            end
        end
        outcome=job_outcome(job,started,output;artifact)
        ncodeunits(Protocol.encode_message(outcome))<=262144 || return job_failure(job,started,"result_payload_limit")
        return outcome
    catch error
        if numerical===nothing
            error isa ExecutionCore.OperationCanceled && return job_failure(job,started,"canceled")
            # A transport read failure cannot be mistaken for no saved result.
            error isa AccessDenied || rethrow()
            return job_failure(job,started,"not_admitted")
        end
        code=lock(service.resources.lock) do
            handle=get(service.resources.handles,job.fence.lease_id,nothing)
            handle!==nothing && handle.request_id==job.request.job_id ? something(handle.failure,"executor-failed") : "executor-failed"
        end
        code=replace(code,"-"=>"_")
        code in ("operation_rejected","canceled","deadline") || (code="executor_failed")
        return job_failure(job,started,code)
    end
end

function run_delivered_job!(service::AgentJobService,flight::AgentJobFlight,connection)
    delivery=flight.delivery
    watching=Ref(true)
    # Progress acks reset only the redelivery timer. The broker is never told
    # that computation succeeded until its exact terminal result is stored.
    heartbeat=@async begin
        next=0.0
        while watching[] && !service.closed && agent_lease_usable(service.ledger,delivery.job.fence)
            if service.ledger.clock()>=next
                try progress_assigned_delivery!(connection,service.ledger,delivery) catch; service.state=:offline end
                next=service.ledger.clock()+5
            end
            sleep(0.05)
        end
    end
    try
        flight.outcome=execute_delivered_job!(service,flight,connection)
        flight.phase=:persisting
        # Persistence may outlast the scientific deadline briefly, but never the
        # lease. Only the already completed outcome is retried in this loop.
        until=Protocol.parse_utc_timestamp(delivery.job.request.deadline)+Second(5)
        while !service.closed && agent_lease_usable(service.ledger,delivery.job.fence) && now(UTC)<until
            try
                persist_assigned_result!(connection,service.ledger,delivery,flight.outcome)
                flight.phase=:completed
                return nothing
            catch error
                flight.phase=:delivery_unavailable
                error isa AccessDenied && break
                sleep(0.25)
            end
        end
        flight.phase=:revoked
        # A stale input must not become authority after recovery. If transport
        # is down, the next bounded pull will reject its expired deadline/fence.
        try terminate_assigned_delivery!(connection,delivery) catch end
    catch
        flight.phase=:delivery_unavailable
        # Leave it unacknowledged. Redelivery either finds the durable result or
        # returns execution_uncertain; it never blindly repeats numerical work.
    finally
        watching[]=false
        wait(heartbeat)
    end
    return nothing
end

function accept_assigned_delivery!(service::AgentJobService,delivery::AssignedDelivery)
    lock(service.lock) do
        service.closed && throw(AccessDenied(409,"Job service is closed"))
        service.connection!==nothing || throw(BrokerUnavailable())
        validate_assigned_delivery(service.connection,delivery)
        job=delivery.job
        agent_lease_usable(service.ledger,job.fence) || throw(AccessDenied(409,"Job assignment is not usable"))
        if haskey(service.flights,job.request.job_id)
            previous=service.flights[job.request.job_id]
            previous.delivery.job==job || throw(AccessDenied(409,"Job identifier was reused with another request"))
            return previous
        end
        length(service.flights)<service.ledger.capacity || throw(AccessDenied(409,"Job delivery capacity is occupied"))
        flight=AgentJobFlight(delivery,nothing,:received,nothing)
        service.flights[job.request.job_id]=flight
        connection=service.connection
        flight.task=@async run_delivered_job!(service,flight,connection)
        return flight
    end
end

function tick_jobs!(service::AgentJobService)
    connection=lock(service.lock) do
        for (id,flight) in collect(service.flights)
            flight.task!==nothing && istaskdone(flight.task) && delete!(service.flights,id)
        end
        service.closed || length(service.flights)>=service.ledger.capacity ? nothing : service.connection
    end
    connection===nothing && return nothing
    service.state=!connection.closed && NATS.status(connection.connection)==NATS.CONNECTED ? :online : :offline
    service.state==:online || return nothing
    delivery=poll_assigned_job!(connection,service.ledger)
    delivery===nothing || accept_assigned_delivery!(service,delivery)
    return nothing
end

function connect_jobs!(service::AgentJobService)
    lock(service.lock) do
        service.closed && return
        service.connection===nothing || return
        service.connector!==nothing && !istaskdone(service.connector) && return
        service.ledger.clock()>=service.retry_at || return
        service.retry_at=service.ledger.clock()+2
        service.connector=@async begin
            connection=nothing
            try
                connection=BrokerJobs(service.endpoint,WorkerIdentity(service.ledger.worker_id))
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

"""Start bounded job polling only after the shared resource owner's recovery completes."""
function start_jobs!(service::AgentJobService)
    lock(service.lock) do
        service.closed && throw(ArgumentError("Job service is closed"))
        service.resources.recovered || throw(ArgumentError("Job resources are not recovered"))
        service.task===nothing || return service
        service.task=@async while !service.closed
            try connect_jobs!(service); tick_jobs!(service) catch; service.state=:offline end
            sleep(0.05)
        end
    end
    return service
end

function stop_jobs!(service::AgentJobService)
    lock(service.lock) do; service.closed=true; end
    service.task===nothing || wait(service.task;throw=false)
    service.connector===nothing || wait(service.connector;throw=false)
    deliveries=lock(()->[flight.delivery for flight in values(service.flights)],service.lock)
    for delivery in deliveries
        try
            cancel_assigned!(service.resources,delivery.job.fence,delivery.job.request.job_id)
        catch error
            # The root agent may already have revoked the lease. Its resource
            # guard/cleanup then owns termination; do not attempt new authority.
            error isa AccessDenied || rethrow()
        end
    end
    return nothing
end

function Base.close(service::AgentJobService)
    stop_jobs!(service)
    tasks=lock(()->[flight.task for flight in values(service.flights) if flight.task!==nothing],service.lock)
    timedwait(()->all(istaskdone,tasks),30;pollint=0.025)==:ok ||
        throw(ArgumentError("Job execution cleanup remains unresolved"))
    foreach(wait,tasks)
    connection=lock(service.lock) do
        prior=service.connection; service.connection=nothing; empty!(service.flights); prior
    end
    connection===nothing || close(connection)
    service.state=:stopped
    return nothing
end
