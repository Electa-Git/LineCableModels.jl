"""Retain one terminal request until its bounded reply or explicit uncertain outcome."""
mutable struct TerminalFlight
    "Original private command while sending; removed after completion or timeout."
    command::Union{Nothing,Protocol.TerminalCommand}
    "Request identity retained without the input payload."
    request_id::String
    "Exact monotonic assignment command revision."
    revision::Int
    "Action needed to validate its reply."
    action::String
    "Requested stream identity, absent for open."
    session_id::Union{Nothing,String}
    "Requested input sequence, or zero."
    input_sequence::Int
    "Requested output cursor, or zero."
    after::Int
    "Digest for an exact explicit retry, not raw input."
    digest::Vector{UInt8}
    "Finite reply deadline on the coordinator's monotonic clock."
    deadline::Float64
    "Next resend time for read-only queries and idempotent open only."
    resend_at::Float64
    "Owned finite send task, separate from heartbeat scheduling."
    sender::Union{Nothing,Task}
    "First matching bounded reply, never a lease or readiness grant."
    report::Union{Nothing,Protocol.TerminalReport}
end
Base.show(io::IO,::TerminalFlight)=print(io,"TerminalFlight(<private>)")

"""Keep one terminal assignment's command revision and latest private flight."""
mutable struct TerminalCommandLane
    "Full live assignment identity."
    fence::AssignmentFence
    "Latest allocated command revision."
    revision::Int
    "At most one pending request and bounded result."
    flight::Union{Nothing,TerminalFlight}
end
Base.show(io::IO,::TerminalCommandLane)=print(io,"TerminalCommandLane(<private>)")

"""
    TerminalCoordinator(endpoint, coordinator, worker_ids)

Own transient private terminal requests under the existing run/lease authority.
An administrator may manage inventory but cannot read or write another owner's
terminal. Commands and replies are bounded, fenced and request-correlated. Input,
restart and other mutations are never automatically replayed; explicit uncertain
retries preserve their original revision and digest. No terminal bytes persist
in SQLite, JetStream, scientific job records or ordinary diagnostics.
"""
mutable struct TerminalCoordinator{C<:LeaseCoordinator}
    "Server-owned broker settings."
    endpoint::BrokerEndpoint
    "Existing run/assignment authority, not a second lease store."
    coordinator::C
    "Provisioned worker identities."
    workers::Vector{String}
    "Separate exact-assignment terminal connection."
    connection::Union{Nothing,BrokerTerminal{CoordinatorIdentity}}
    "Bounded volatile request lanes."
    lanes::Dict{String,TerminalCommandLane}
    "At most one live browser attachment per private assignment."
    attachments::Dict{String,TerminalAttachment}
    "Independent scheduling task."
    task::Union{Nothing,Task}
    "Independent finite connection attempt."
    connector::Union{Nothing,Task}
    "Next allowed connection attempt in local monotonic seconds."
    retry_at::Float64
    "Online, offline or stopped; not REPL readiness."
    state::Symbol
    "Permanently reject new use after closure."
    closed::Bool
    "Short request, reply and connection ownership updates."
    lock::ReentrantLock
end
function TerminalCoordinator(endpoint::BrokerEndpoint,coordinator::LeaseCoordinator,worker_ids)
    workers=sort!(unique(Protocol.runtime_token.(collect(worker_ids))))
    1<=length(workers)<=128 || throw(ArgumentError("terminal channel requires 1:128 provisioned workers"))
    TerminalCoordinator(endpoint,coordinator,workers,nothing,Dict{String,TerminalCommandLane}(),Dict{String,TerminalAttachment}(),nothing,nothing,
        0.0,:offline,false,ReentrantLock())
end
Base.show(io::IO,::TerminalCoordinator)=print(io,"TerminalCoordinator(<private>)")
terminal_remote_clock(service)=lease_clock(service.coordinator)

function remote_terminal_fence(service::TerminalCoordinator,principal::Principal,id::UUID)
    store=service.coordinator.assignments.inventory.store
    lease=get_assignment(store,principal,id)
    lease.fence.owner==principal.id || throw(AccessDenied(404,"Terminal assignment not found"))
    assignment_usable(service.coordinator,principal,id) || throw(AccessDenied(409,"Terminal assignment is not usable"))
    profile=service.coordinator.assignments.inventory.profiles.definitions[lease.fence.profile_id]
    profile.kind==:terminal && profile.isolation==:container || throw(AccessDenied(400,"Terminal profile required"))
    lease.fence.worker_id in service.workers || throw(AccessDenied(409,"Terminal worker is not provisioned"))
    return lease.fence
end

function terminal_request_digest(command)
    SHA.sha256(JSON3.write((action=command.action,session_id=command.session_id,writer_id=command.writer_id,
        input_sequence=command.input_sequence,after=command.after,columns=command.columns,rows=command.rows,bytes=command.bytes)))
end

function send_terminal_flight!(service,lane,flight;send=nothing)
    command=flight.command
    command === nothing && return
    flight.sender !== nothing && !istaskdone(flight.sender) && return
    flight.resend_at=terminal_remote_clock(service)+0.25
    flight.sender=@async begin
        try
            lock(()->service.closed,service.lock) && throw(BrokerUnavailable())
            remote_terminal_fence(service,Principal(lane.fence.owner),UUID(lane.fence.lease_id))==lane.fence || throw(BrokerUnavailable())
            if send === nothing
                connection=lock(()->service.connection,service.lock)
                connection === nothing && throw(BrokerUnavailable())
                watch_terminal!(connection,lane.fence)
                send_terminal!(connection,command)
            else
                send(command) # finite test-only transport substitution
            end
        catch
            lock(service.lock) do
                flight.deadline=min(flight.deadline,terminal_remote_clock(service))
                service.state=:offline
            end
        end
        return nothing
    end
    return nothing
end

"""
    request_terminal!(service, principal, lease_id, action; session_id=nothing,
        writer_id=nothing, input_sequence=0, after=0, columns=0, rows=0,
        bytes=UInt8[], request_id=uuid4(), retry=false)

Authorize one private terminal action and return its bounded asynchronous flight.
Only an exact explicit retry can resend an uncertain mutation; it retains the
original revision. A later query may supersede an expired request, fencing its
late reply. Input sequence and stream checks remain mandatory at the agent.
"""
function request_terminal!(service::TerminalCoordinator,principal::Principal,id::UUID,action::AbstractString;
        session_id=nothing,writer_id=nothing,input_sequence=0,after=0,columns=0,rows=0,
        bytes::AbstractVector{UInt8}=UInt8[],request_id=uuid4(),retry::Bool=false,send=nothing)
    fence=remote_terminal_fence(service,principal,id)
    length(bytes)<=Protocol.MAX_TERMINAL_CHUNK_BYTES || throw(AccessDenied(400,"Terminal input exceeds its bound"))
    request=string(request_id)
    candidate=Protocol.validate(Protocol.TerminalCommand("2.0",request,fence,1,String(action),session_id,writer_id,
        input_sequence,after,columns,rows,collect(bytes)))
    digest=terminal_request_digest(candidate)
    return lock(service.lock) do
        service.closed && throw(AccessDenied(503,"Terminal channel is stopped"))
        timestamp=terminal_remote_clock(service)
        lane=get(service.lanes,fence.lease_id,nothing)
        if lane === nothing
            length(service.lanes)<min(256,service.coordinator.assignments.limits.total) || throw(CapacityUnavailable())
            lane=TerminalCommandLane(fence,0,nothing);service.lanes[fence.lease_id]=lane
        end
        lane.fence==fence || throw(AccessDenied(409,"Terminal assignment changed"))
        prior=lane.flight
        if prior !== nothing && prior.request_id==request
            prior.digest==digest || throw(AccessDenied(409,"Terminal retry differs"))
            if retry && prior.report === nothing && timestamp>=prior.deadline
                prior.sender === nothing || istaskdone(prior.sender) || throw(AccessDenied(409,"Terminal send is still pending"))
                prior.command=Protocol.TerminalCommand(candidate.protocol_version,request,fence,prior.revision,candidate.action,
                    session_id,writer_id,input_sequence,after,columns,rows,collect(bytes))
                prior.deadline=timestamp+(action=="restart" ? 30 : 5)
                send_terminal_flight!(service,lane,prior;send)
            end
            return prior
        end
        if prior !== nothing && (prior.sender !== nothing && !istaskdone(prior.sender) ||
                prior.report === nothing && timestamp<prior.deadline)
            throw(AccessDenied(409,"Terminal action is still pending"))
        end
        retry && throw(AccessDenied(409,"Terminal retry is no longer current"))
        lane.revision<9_007_199_254_740_991 || throw(AccessDenied(409,"Terminal command sequence exhausted"))
        command=Protocol.TerminalCommand(candidate.protocol_version,request,fence,lane.revision+1,candidate.action,
            session_id,writer_id,input_sequence,after,columns,rows,collect(bytes))
        flight=TerminalFlight(command,request,command.revision,command.action,session_id,input_sequence,after,digest,
            timestamp+(action=="restart" ? 30 : 5),timestamp,nothing,nothing)
        lane.revision=command.revision;lane.flight=flight
        send_terminal_flight!(service,lane,flight;send)
        return flight
    end
end

function accept_terminal_report!(service::TerminalCoordinator,worker::String,report::Protocol.TerminalReport)
    Protocol.validate(report)
    worker==report.fence.worker_id && worker in service.workers || throw(AccessDenied(403,"Terminal report worker differs"))
    return lock(service.lock) do
        service.closed && return false
        lane=get(service.lanes,report.fence.lease_id,nothing)
        lane !== nothing && lane.fence==report.fence && lane.flight !== nothing || return false
        flight=lane.flight
        report.request_id==flight.request_id && report.revision==flight.revision && flight.report === nothing &&
            terminal_remote_clock(service)<flight.deadline || return false
        remote_terminal_fence(service,Principal(lane.fence.owner),UUID(lane.fence.lease_id))==lane.fence || return false
        if report.accepted
            if flight.action in ("open","restart")
                report.session_id !== nothing || return false
                flight.action=="restart" && report.session_id==flight.session_id && return false
            else
                report.session_id==flight.session_id || return false
            end
            flight.action=="input" && report.input_sequence!=flight.input_sequence && return false
            if flight.action=="read"
                expected=flight.after+length(report.bytes)
                (report.gap ? report.cursor>=expected : report.cursor==expected) || return false
            else
                isempty(report.bytes) && report.cursor==0 && !report.gap || return false
            end
        end
        flight.report=report;flight.command=nothing
        return true
    end
end

"""Read only the currently authorized flight; inspection does not refresh its deadline."""
function terminal_flight(service::TerminalCoordinator,principal::Principal,id::UUID,request_id::UUID)
    fence=remote_terminal_fence(service,principal,id)
    lock(service.lock) do
        lane=get(service.lanes,fence.lease_id,nothing)
        lane !== nothing && lane.fence==fence && lane.flight !== nothing && lane.flight.request_id==string(request_id) ||
            throw(AccessDenied(409,"Terminal request is no longer current"))
        return lane.flight
    end
end

function connect_terminals!(service::TerminalCoordinator)
    lock(service.lock) do
        service.closed && return
        service.connection === nothing || return
        service.connector !== nothing && !istaskdone(service.connector) && return
        terminal_remote_clock(service)>=service.retry_at || return
        service.retry_at=terminal_remote_clock(service)+2
        service.connector=@async begin
            connection=nothing
            try
                connection=BrokerTerminal(service.endpoint,CoordinatorIdentity();worker_ids=service.workers,
                    capacity=min(256,service.coordinator.assignments.limits.total))
                lock(service.lock) do
                    if !service.closed;service.connection=connection;connection=nothing;end
                end
            catch
                service.state=:offline
            finally
                connection === nothing || close(connection)
            end
        end
    end
end

function tick_terminals!(service::TerminalCoordinator)
    connection=lock(()->service.connection,service.lock)
    if connection !== nothing
        service.state=!connection.closed && NATS.status(connection.connection)==NATS.CONNECTED ? :online : :offline
        if service.state==:offline
            # Retire old queued frames on loss, as at the agent. A fresh private
            # connection never replays buffered terminal commands after recovery.
            lock(service.lock) do
                service.connection===connection && (service.connection=nothing)
            end
            close(connection);connection=nothing
        else
            for envelope in poll_terminal!(connection)
                try accept_terminal_report!(service,envelope.worker_id,envelope.record) catch end
            end
        end
    end
    expired=AssignmentFence[]
    lock(service.lock) do
        timestamp=terminal_remote_clock(service)
        for (id,lane) in collect(service.lanes)
            flight=lane.flight
            usable=try remote_terminal_fence(service,Principal(lane.fence.owner),UUID(id))==lane.fence catch;false end
            if flight !== nothing && flight.sender !== nothing && istaskdone(flight.sender)
                flight.sender=nothing
            end
            if !usable && (flight === nothing || flight.sender === nothing)
                delete!(service.lanes,id);push!(expired,lane.fence);continue
            end
            flight === nothing && continue
            if timestamp>=flight.deadline || flight.report !== nothing
                flight.command=nothing
            elseif usable && service.state==:online && flight.action in ("open","status","read") &&
                    flight.sender === nothing && timestamp>=flight.resend_at
                send_terminal_flight!(service,lane,flight)
            end
        end
    end
    if connection !== nothing
        for fence in expired;unwatch_terminal!(connection,fence);end
    end
    return nothing
end

function start_terminals!(service::TerminalCoordinator)
    lock(service.lock) do
        service.closed && throw(ArgumentError("terminal coordinator is closed"))
        service.task === nothing || return service
        service.task=@async while !service.closed
            try connect_terminals!(service);tick_terminals!(service) catch;service.state=:offline end
            sleep(0.02)
        end
    end
    return service
end
function Base.close(service::TerminalCoordinator)
    lock(()->service.closed=true,service.lock)
    attachments=lock(()->collect(values(service.attachments)),service.lock)
    foreach(abort_terminal_attachment!,attachments)
    for attachment in attachments
        attachment.task===current_task() || wait(attachment.task)
    end
    service.task === nothing || wait(service.task)
    service.connector === nothing || wait(service.connector)
    for lane in values(service.lanes)
        lane.flight === nothing || lane.flight.sender === nothing || wait(lane.flight.sender)
    end
    connection=lock(service.lock) do
        previous=service.connection;service.connection=nothing;empty!(service.lanes);previous
    end
    connection === nothing || close(connection)
    service.state=:stopped
    return nothing
end
