"""Retain one bounded transient request and its locally expiring response."""
mutable struct ScientificFlight
    "Complete ordered command sent to one worker incarnation."
    command::Protocol.ScientificCommand
    "Local monotonic time before send, also used to subtract round-trip latency."
    sent_at::Float64
    "Finite response deadline in local monotonic seconds."
    deadline::Float64
    "Matching report, if accepted before deadline."
    report::Union{Nothing,Protocol.ScientificReport}
    "Local readiness expiry; never renewed merely by reading a report."
    ready_until::Float64
end

"""Bound command sequencing and mutation retry identities within one lease."""
mutable struct ScientificCommandLane
    "Exact assignment identity."
    fence::Protocol.AssignmentFence
    "Last allocated command revision."
    revision::Int
    "At most one outstanding command/report."
    flight::Union{Nothing,ScientificFlight}
    "Last correlated report while a read-only successor query is outstanding."
    observed::Union{Nothing,ScientificFlight}
    "Bounded mutation request-ID/digest tombstones; no old inputs retained."
    mutations::Dict{String,String}
end

"""
    ScientificCoordinator(endpoint, coordinator, worker_ids)

Own remote preparation/status/cancel requests over a separate bounded connection.
SQLite/LeaseCoordinator remain the sole assignment authority. Reports are volatile,
command-correlated and locally expiring; persisted rows or worker advertisements
never establish readiness. No scientific package is imported here.
"""
mutable struct ScientificCoordinator{C<:LeaseCoordinator}
    "Operator-owned connection settings."
    endpoint::BrokerEndpoint
    "Existing run/lease authorization owner."
    coordinator::C
    "Exact provisioned worker subscriptions."
    workers::Vector{String}
    "Independent scientific-record connection."
    connection::Union{Nothing,BrokerControl{CoordinatorIdentity}}
    "Current assignment command lanes, bounded by admitted capacity."
    lanes::Dict{String,ScientificCommandLane}
    "Owned response scheduler."
    task::Union{Nothing,Task}
    "Independent finite connection attempt."
    connector::Union{Nothing,Task}
    "Next initial retry time in local monotonic seconds."
    retry_at::Float64
    "Online, offline or stopped; not model readiness."
    state::Symbol
    "Reject new requests after shutdown begins."
    closed::Bool
    "Serialize short requests and matching report updates."
    lock::ReentrantLock
end

function ScientificCoordinator(endpoint::BrokerEndpoint,coordinator::LeaseCoordinator,worker_ids)
    ids = sort!(unique(Protocol.runtime_token.(collect(worker_ids))))
    1 <= length(ids) <= 128 || throw(ArgumentError("scientific channel requires 1:128 provisioned workers"))
    ScientificCoordinator(endpoint,coordinator,ids,nothing,Dict{String,ScientificCommandLane}(),
        nothing,nothing,0.0,:offline,false,ReentrantLock())
end
science_clock(service::ScientificCoordinator) = lease_clock(service.coordinator)

function remote_scientific_fence(service::ScientificCoordinator,principal::Principal,id::UUID)
    store = service.coordinator.assignments.inventory.store
    lease = get_assignment(store,principal,id)
    # get_assignment protects owner boundaries, including before availability
    # checks. Historical or merely reserved rows cannot authorize execution.
    assignment_usable(service.coordinator,principal,id) || throw(AccessDenied(409,"Scientific assignment is not usable"))
    profile = service.coordinator.assignments.inventory.profiles.definitions[lease.fence.profile_id]
    profile.kind == :scientific || throw(AccessDenied(400,"Scientific profile required"))
    lease.fence.worker_id in service.workers || throw(AccessDenied(409,"Scientific worker is not provisioned"))
    return lease.fence
end

"""
    request_scientific!(service, principal, lease_id, action;
        parameters=Dict(), target_id=nothing, request_id=uuid4())

Authorize and send one explicit transient request without waiting for scientific
work. Preparation/current-request mutations are never automatically replayed after an uncertain reply. Exact
HTTP retries do not execute them twice; changed inputs with a reused ID fail.
Status coalesces a pending query and never carries preparation inputs. Return
the current flight for a caller to describe/poll through the same owned API.
The idempotent `cancel_job` tombstone may be resent after a lost reply, using a
new command revision but the same request/target identity; it cannot run code.
"""
function request_scientific!(service::ScientificCoordinator,principal::Principal,id::UUID,action::String;
        parameters=Dict{String,Any}(),target_id=nothing,request_id=uuid4(),
        send=command->send_control!(service.connection,command))
    fence = remote_scientific_fence(service,principal,id)
    normalized = Protocol.normalize_wire(parameters)
    normalized isa Dict{String,Any} || throw(AccessDenied(400,"Expected scientific input object"))
    request = string(request_id)
    Protocol.runtime_uuid(request)
    # Validate even a coalesced status or replayed mutation before consulting
    # transient state. Coalescing must not make malformed inputs acceptable.
    Protocol.validate(Protocol.ScientificCommand("2.0",request,fence,1,action,normalized,target_id))
    digest = Protocol.input_hash("runtime.scientific_control",Dict("action"=>action,"parameters"=>normalized,"target_id"=>target_id))
    return lock(service.lock) do
        service.closed && throw(AccessDenied(503,"Scientific channel is stopped"))
        timestamp = science_clock(service)
        lane = get(service.lanes,fence.lease_id,nothing)
        if lane === nothing
            length(service.lanes) < service.coordinator.assignments.limits.total ||
                throw(AccessDenied(409,"Scientific channel capacity is occupied"))
            lane = ScientificCommandLane(fence,0,nothing,nothing,Dict{String,String}())
            service.lanes[fence.lease_id] = lane
        end
        lane.fence == fence || throw(AccessDenied(409,"Scientific assignment changed"))
        if haskey(lane.mutations,request)
            lane.mutations[request] == digest || throw(AccessDenied(409,"Scientific request identity was reused"))
            action=="cancel_job" || return lane.flight # preparation is never repeated automatically
            for prior in (lane.flight,lane.observed)
                prior===nothing && continue
                prior.command.request_id==request && prior.command.action=="cancel_job" || continue
                prior.report!==nothing && prior.report.accepted && return prior
                prior.report===nothing && timestamp<prior.deadline && return prior
            end
        end
        previous = lane.flight
        if previous !== nothing && previous.command.request_id == request
            command = previous.command
            (command.action,command.parameters,command.target_id) == (action,normalized,target_id) ||
                throw(AccessDenied(409,"Scientific request identity was reused"))
            action=="cancel_job" || return previous
        end
        if action == "status" && previous !== nothing && previous.report === nothing && timestamp < previous.deadline
            return previous
        end
        lane.revision < 9_007_199_254_740_991 || throw(AccessDenied(409,"Scientific command sequence exhausted"))
        command = Protocol.validate(Protocol.ScientificCommand("2.0",request,fence,lane.revision+1,action,normalized,target_id))
        if action != "status"
            haskey(lane.mutations,request) || length(lane.mutations) < 256 || throw(AccessDenied(409,"Scientific command history requires a new assignment"))
            lane.mutations[request] = digest
            lane.observed = nothing # mutations immediately invalidate prior evidence
        end
        lane.revision += 1
        flight = ScientificFlight(command,timestamp,timestamp+5,nothing,0.0)
        lane.flight = flight # correlate before send, including failed send
        try
            send(command)
        catch
            service.state=:offline
            # Keep the command identity but make the finite failure visible.
            # A status query or explicit new mutation may follow; no replay.
            flight.deadline=timestamp
            throw(BrokerUnavailable())
        end
        return flight
    end
end

function accept_scientific_report!(service::ScientificCoordinator,worker::String,report::Protocol.ScientificReport)
    Protocol.validate(report)
    worker == report.fence.worker_id && worker in service.workers || throw(AccessDenied(403,"Scientific report worker differs"))
    return lock(service.lock) do
        service.closed && return false
        lane = get(service.lanes,report.fence.lease_id,nothing)
        lane !== nothing && lane.fence == report.fence && lane.flight !== nothing || return false
        flight = lane.flight
        command = flight.command
        command.request_id == report.request_id && command.revision == report.revision || return false
        report.preparation == "ready" && command.action != "status" && return false
        # First correlated reply only. A duplicate cannot extend a deadline or
        # replace a prior answer with newly invented readiness.
        flight.report===nothing && science_clock(service)<flight.deadline || return false
        remote_scientific_fence(service,Principal(report.fence.owner),UUID(report.fence.lease_id)) == report.fence || return false
        flight.report=report
        flight.ready_until = report.preparation=="ready" ? min(flight.deadline,flight.sent_at+report.valid_for_ms/1000) : 0.0
        lane.observed=flight
        return true
    end
end

"""
    remote_scientific_status(service, principal, lease_id)

Read the current command/report without preparing, polling the child or extending
readiness. Return unknown when no fresh report exists. Authorization and live
lease checks apply on every read, independently of cached report content.
"""
function remote_scientific_status(service::ScientificCoordinator,principal::Principal,id::UUID)
    fence = remote_scientific_fence(service,principal,id)
    lock(service.lock) do
        lane = get(service.lanes,fence.lease_id,nothing)
        flight = lane===nothing ? nothing : lane.flight
        observed = lane===nothing ? nothing : lane.observed
        timestamp = science_clock(service)
        report = observed===nothing || timestamp>=observed.deadline ? nothing : observed.report
        remaining_ms = observed===nothing ? 0 : max(0,floor(Int,1000*(observed.ready_until-timestamp)))
        ready = !service.closed && service.state==:online && report!==nothing && report.preparation=="ready" && remaining_ms>0
        return (channel=String(service.state),phase=report===nothing ? "idle" : report.phase,
            preparation=ready ? "ready" : report===nothing || report.preparation=="ready" ? "unknown" : report.preparation,
            request_id=flight===nothing ? nothing : flight.command.request_id,
            revision=flight===nothing ? 0 : flight.command.revision,
            pending=flight!==nothing && flight.report===nothing && timestamp<flight.deadline,
            accepted=report===nothing ? nothing : report.accepted,
            reason=report===nothing ? "report_unavailable" : report.reason,
            executor_id=report===nothing ? nothing : report.executor_id,
            executor_generation=report===nothing ? 0 : report.executor_generation,
            current_request_id=report===nothing ? nothing : report.current_request_id,
            preparation_key=ready ? report.preparation_key : nothing,
            valid_for_ms=ready ? remaining_ms : 0,
            progress=report===nothing ? 0.0 : report.progress_milli/1000,
            elapsed_seconds=report===nothing ? 0.0 : report.elapsed_ms/1000,
            output_lines=report===nothing ? 0 : report.output_lines,
            failure=report===nothing ? nothing : report.failure)
    end
end

function tick_science!(service::ScientificCoordinator)
    lock(service.lock) do
        for (id,lane) in collect(service.lanes)
            usable = try remote_scientific_fence(service,Principal(lane.fence.owner),UUID(id))==lane.fence catch; false end
            usable || delete!(service.lanes,id)
        end
    end
    connection = lock(()->service.connection,service.lock)
    connection===nothing && return nothing
    service.state = NATS.status(connection.connection)==NATS.CONNECTED ? :online : :offline
    for envelope in poll_control!(connection)
        try accept_scientific_report!(service,envelope.worker_id,envelope.record) catch end
    end
    return nothing
end

function connect_science!(service::ScientificCoordinator)
    lock(service.lock) do
        service.closed && return
        service.connection===nothing || return
        service.connector!==nothing && !istaskdone(service.connector) && return
        science_clock(service)>=service.retry_at || return
        service.retry_at=science_clock(service)+2
        service.connector = @async begin
            connection=nothing
            try
                connection=BrokerControl(service.endpoint,CoordinatorIdentity();worker_ids=service.workers,traffic=:science)
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

function start_science!(service::ScientificCoordinator)
    lock(service.lock) do
        service.closed && throw(ArgumentError("scientific coordinator is closed"))
        service.task===nothing || return service
        service.task = @async while !service.closed
            try connect_science!(service); tick_science!(service) catch; service.state=:offline end
            sleep(0.05)
        end
    end
    return service
end

function Base.close(service::ScientificCoordinator)
    lock(service.lock) do; service.closed=true; end
    service.task===nothing || wait(service.task)
    service.connector===nothing || wait(service.connector)
    connection = lock(service.lock) do
        previous=service.connection; service.connection=nothing; empty!(service.lanes); previous
    end
    connection===nothing || close(connection)
    service.state=:stopped
    return nothing
end
