"""Retain command ordering and at most one child-inspection task per assignment."""
mutable struct AgentScientificLane
    "Exact assignment; never reused across a worker/run generation."
    fence::AssignmentFence
    "Last accepted command, including passive inputs for duplicate comparison."
    command::Protocol.ScientificCommand
    "Whether the last mutation/query was admitted; retained for duplicate replies."
    accepted::Bool
    "Owned asynchronous inspection/report task, if any."
    task::Union{Nothing,Task}
end

"""
    AgentScientificService(endpoint, resources, ledger)

Own transient preparation/status/cancel traffic separately from heartbeat and
lease control. Construction is inert. Every command checks the existing live
lease and enters ScientificResources; this service cannot bypass its required
physical driver or preparation orchestration. Durable jobs use another channel.
"""
mutable struct AgentScientificService{R<:ScientificResources,L<:AgentLeaseLedger}
    "Operator-owned broker endpoint."
    endpoint::BrokerEndpoint
    "Existing root scientific resource owner."
    resources::R
    "The same worker-incarnation lease authority as the agent."
    ledger::L
    "Separate bounded scientific-record connection."
    connection::Union{Nothing,BrokerControl{WorkerIdentity}}
    "Per-assignment ordering and inspection ownership."
    lanes::Dict{String,AgentScientificLane}
    "Owned finite-poll scheduler."
    task::Union{Nothing,Task}
    "Independent bounded connection attempt."
    connector::Union{Nothing,Task}
    "Next allowed initial connection attempt in local monotonic seconds."
    retry_at::Float64
    "Connection dimension, not preparation readiness."
    state::Symbol
    "Stop accepting new commands."
    closed::Bool
    "Serialize short command/connection state updates."
    lock::ReentrantLock
end

function AgentScientificService(endpoint::BrokerEndpoint,resources::ScientificResources,ledger::AgentLeaseLedger)
    resources.ledger === ledger || throw(ArgumentError("scientific channel must share agent lease authority"))
    AgentScientificService(endpoint,resources,ledger,nothing,Dict{String,AgentScientificLane}(),
        nothing,nothing,0.0,:offline,false,ReentrantLock())
end

function scientific_report(service::AgentScientificService,command;accepted=true,reason="accepted",fresh=false)
    # Keep the status identity and its remaining validity under the same lock;
    # a replacement executor cannot lend its TTL to a previous snapshot.
    lock(service.resources.lock) do
    status = scientific_status(service.resources,command.fence)
    handle = get(service.resources.handles,command.fence.lease_id,nothing)
    inspecting = handle !== nothing && handle.request_kind == :inspect
    phase = status.phase == :starting && handle !== nothing ?
        (inspecting ? "idle" : handle.request_kind == :operation ? "executing" : "starting") : String(status.phase)
    valid_ms = 0
    if fresh && accepted && status.preparation == :ready
        valid_ms = lock(service.resources.lock) do
            handle = get(service.resources.handles,command.fence.lease_id,nothing)
            handle === nothing && return 0
            lock(service.ledger.lock) do
                agent_lease_usable(service.ledger,command.fence) || return 0
                lease = service.ledger.leases[agent_key(command.fence)]
                remaining = min(5.0,handle.prepared_until-service.ledger.clock(),lease.expires_at-service.ledger.clock(),
                    service.ledger.probed_at+service.ledger.presence_seconds-service.ledger.clock())
                max(0,floor(Int,1000*remaining))
            end
        end
    end
    preparation = status.preparation == :ready ? (valid_ms>0 ? "ready" : "unknown") : String(status.preparation)
    failure = status.failure === nothing ? nothing : replace(status.failure,"-"=>"_")
    report = Protocol.ScientificReport("2.0",command.request_id,command.fence,command.revision,accepted,reason,
        phase,preparation,status.executor_id,status.generation,inspecting ? nothing : status.request_id,
        valid_ms>0 ? status.preparation_key : nothing,clamp(round(Int,status.progress*1000),0,1000),
        clamp(floor(Int,status.elapsed_seconds*1000),0,9_007_199_254_740_991),clamp(status.output_lines,0,1_000_000),failure,valid_ms)
    return Protocol.validate(report)
    end
end

function unavailable_scientific_report(command,reason)
    Protocol.ScientificReport("2.0",command.request_id,command.fence,command.revision,false,reason,
        "idle","unknown",nothing,0,nothing,nothing,0,0,0,nothing,0)
end

function send_scientific_report!(service::AgentScientificService,report)
    connection = lock(()->service.connection,service.lock)
    !service.closed && connection !== nothing || throw(BrokerUnavailable())
    send_control!(connection,report)
end

"""
    receive_scientific!(service, command)

Admit one ordered scientific control message. Preparation returns acceptance
without waiting for numerical work. Status queries retained child state only when
idle and prepared; it never performs implicit preparation. Duplicate mutations
do not execute twice, and old revisions cannot cancel or replace newer work.
"""
function receive_scientific!(service::AgentScientificService,command::Protocol.ScientificCommand;
        emit=report->send_scientific_report!(service,report))
    Protocol.validate(command)
    lock(service.lock) do
        service.closed && throw(AccessDenied(409,"Scientific channel is closed"))
        scientific_authority(service.resources,command.fence)
        lane = get(service.lanes,command.fence.lease_id,nothing)
        if lane === nothing
            length(service.lanes) < service.ledger.capacity || throw(AccessDenied(409,"Scientific channel capacity is occupied"))
            lane = AgentScientificLane(command.fence,command,true,nothing)
            service.lanes[command.fence.lease_id] = lane
        else
            lane.fence == command.fence || throw(AccessDenied(409,"Scientific channel fence differs"))
            if command.revision < lane.command.revision
                emit(unavailable_scientific_report(command,"stale_command")); return nothing
            elseif command.revision == lane.command.revision
                lane.command == command || throw(AccessDenied(409,"Scientific command identity was reused"))
                # Never repeat a preparation/cancel mutation or replay cached
                # ready evidence. A later status request can inspect the child.
                emit(scientific_report(service,command;accepted=lane.accepted,reason=lane.accepted ? "duplicate" : "not_admitted")); return nothing
            end
            lane.command = command
            lane.accepted=true
        end
        if command.action == "cancel_job"
            try
                running=cancel_job!(service.resources,command.fence,command.target_id)
                emit(scientific_report(service,command;reason=running ? "cancel_requested" : "cancel_recorded"))
            catch error
                error isa BrokerUnavailable && rethrow()
                lane.accepted=false
                emit(unavailable_scientific_report(command,"cancellation_not_admitted"))
            end
        elseif command.action == "cancel"
            canceled = cancel_assigned!(service.resources,command.fence,command.target_id)
            emit(scientific_report(service,command;reason=canceled ? "cancel_requested" : "not_pending"))
        elseif command.action == "prepare"
            try
                prepare_assigned!(service.resources,command.fence,command.parameters;request_id=command.request_id)
                emit(scientific_report(service,command))
            catch error
                error isa BrokerUnavailable && rethrow()
                lane.accepted=false
                emit(unavailable_scientific_report(command,"preparation_not_admitted"))
            end
        else
            status = scientific_status(service.resources,command.fence)
            if status.phase != :idle || status.preparation != :ready ||
                    (lane.task !== nothing && !istaskdone(lane.task))
                emit(scientific_report(service,command)); return nothing
            end
            lane.task = @async begin
                report = try
                    fetch(refresh_preparation!(service.resources,command.fence))
                    scientific_report(service,command;fresh=true)
                catch
                    unavailable_scientific_report(command,"inspection_failed")
                end
                # A canceled, superseded or revoked query cannot publish fresh
                # readiness after a later command or resource cleanup.
                lock(service.lock) do
                    if !service.closed && lane.command == command && agent_lease_usable(service.ledger,command.fence)
                        lane.accepted=report.accepted
                        try emit(report) catch; service.state=:offline end
                    end
                end
            end
        end
    end
    return nothing
end

function tick_science!(service::AgentScientificService)
    lock(service.lock) do
        for (id,lane) in collect(service.lanes)
            if !agent_lease_usable(service.ledger,lane.fence) && (lane.task===nothing || istaskdone(lane.task))
                delete!(service.lanes,id)
            end
        end
    end
    connection = lock(()->service.connection,service.lock)
    connection === nothing && return nothing
    service.state = NATS.status(connection.connection)==NATS.CONNECTED ? :online : :offline
    for envelope in poll_control!(connection;limit=16)
        try
            receive_scientific!(service,envelope.record)
        catch
            try send_scientific_report!(service,unavailable_scientific_report(envelope.record,"command_rejected")) catch end
        end
    end
    return nothing
end

function connect_science!(service::AgentScientificService)
    lock(service.lock) do
        service.closed && return
        service.connection === nothing || return
        service.connector !== nothing && !istaskdone(service.connector) && return
        service.ledger.clock() >= service.retry_at || return
        service.retry_at = service.ledger.clock()+2
        service.connector = @async begin
            connection = nothing
            try
                connection = BrokerControl(service.endpoint,WorkerIdentity(service.ledger.worker_id);traffic=:science)
                lock(service.lock) do
                    if !service.closed
                        service.connection=connection; connection=nothing
                    end
                end
            catch
                service.state=:offline
            finally
                connection===nothing || close(connection)
            end
        end
    end
end

function start_science!(service::AgentScientificService)
    lock(service.lock) do
        service.closed && throw(ArgumentError("scientific channel is closed"))
        service.resources.recovered || throw(ArgumentError("scientific resources are not recovered"))
        service.task===nothing || return service
        service.task = @async while !service.closed
            try connect_science!(service); tick_science!(service) catch; service.state=:offline end
            sleep(0.05)
        end
    end
    return service
end

function stop_science!(service::AgentScientificService)
    lock(service.lock) do; service.closed=true; end
    service.task===nothing || wait(service.task;throw=false)
    service.connector===nothing || wait(service.connector;throw=false)
    connection = lock(service.lock) do
        previous=service.connection; service.connection=nothing; previous
    end
    connection===nothing || close(connection)
    service.state=:stopped
    return nothing
end

function Base.close(service::AgentScientificService)
    stop_science!(service)
    tasks = lock(()->[lane.task for lane in values(service.lanes) if lane.task!==nothing],service.lock)
    timedwait(()->all(istaskdone,tasks),10;pollint=0.025) == :ok ||
        throw(ArgumentError("scientific inspection cleanup remains unresolved"))
    foreach(wait,tasks)
    empty!(service.lanes)
    return nothing
end
