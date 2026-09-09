"""Retain one ordered terminal request without retaining its input payload."""
mutable struct AgentTerminalLane
    "Full assignment identity, never reused for another incarnation."
    fence::AssignmentFence
    "Latest admitted command revision."
    revision::Int
    "Latest explicit command UUID."
    request_id::String
    "Digest for exact retry comparison; not raw input bytes."
    digest::Vector{UInt8}
    "At most one asynchronous action, independent of other assignments."
    task::Union{Nothing,Task}
    "Latest bounded reply, retained only for an exact retry."
    report::Union{Nothing,Protocol.TerminalReport}
end
Base.show(io::IO,::AgentTerminalLane)=print(io,"AgentTerminalLane(<private>)")

"""
    AgentTerminalService(endpoint, resources, ledger)

Own private terminal relay scheduling independently of heartbeat and scientific
traffic. Every action passes through TerminalResources and its exact live ledger.
Per-assignment command revisions and digests fence duplicate/late actions; no
input is submitted to scientific jobs or evaluated in the scheduler.
"""
mutable struct AgentTerminalService{R<:TerminalResources,L<:AgentLeaseLedger}
    "Operator-provisioned private broker endpoint."
    endpoint::BrokerEndpoint
    "Existing lease-bound terminal session owner."
    resources::R
    "The same authoritative ledger as the root agent."
    ledger::L
    "Separate bounded terminal connection."
    connection::Union{Nothing,BrokerTerminal{WorkerIdentity}}
    "Ordered action state bounded by occupied leases."
    lanes::Dict{String,AgentTerminalLane}
    "Independent finite initial-connection task."
    connector::Union{Nothing,Task}
    "Independent bounded polling task."
    task::Union{Nothing,Task}
    "Next allowed connect attempt on the agent's clock."
    retry_at::Float64
    "Offline, online or stopped; never REPL readiness."
    state::Symbol
    "Reject new commands once closing starts."
    closed::Bool
    "Short state and connection ownership updates."
    lock::ReentrantLock
end
function AgentTerminalService(endpoint::BrokerEndpoint,resources::TerminalResources,ledger::AgentLeaseLedger)
    resources.ledger === ledger || throw(ArgumentError("terminal channel must share agent lease authority"))
    AgentTerminalService(endpoint,resources,ledger,nothing,Dict{String,AgentTerminalLane}(),nothing,nothing,
        0.0,:offline,false,ReentrantLock())
end
Base.show(io::IO,::AgentTerminalService)=print(io,"AgentTerminalService(<private>)")

terminal_command_digest(command)=SHA.sha256(Protocol.encode_message(command))

function terminal_rejection(command,reason)
    Protocol.TerminalReport("2.0",command.request_id,command.fence,command.revision,false,reason,
        command.session_id,"unknown",false,0,0,0,false,UInt8[],nothing,false)
end

function terminal_action_report(resources,command,id;output=nothing)
    status=terminal_status(resources,command.fence,id)
    bytes=output === nothing ? UInt8[] : output.bytes
    cursor=output === nothing ? 0 : output.cursor
    sequence=output === nothing ? status.output_sequence : output.sequence
    Protocol.TerminalReport("2.0",command.request_id,command.fence,command.revision,true,"accepted",id,
        string(status.phase),status.writer_connected,status.input_sequence,cursor,sequence,
        output !== nothing && output.gap,bytes,status.failure === nothing ? nothing : string(status.failure),status.cleanup_pending)
end

function terminal_action!(::Val{:open},resources,command)
    id=open_terminal!(resources,command.fence,command.writer_id;columns=command.columns,rows=command.rows)
    terminal_action_report(resources,command,id)
end
terminal_action!(::Val{:status},resources,command)=terminal_action_report(resources,command,command.session_id)
function terminal_action!(::Val{:read},resources,command)
    output=read_terminal(resources,command.fence,command.session_id,command.after)
    terminal_action_report(resources,command,command.session_id;output)
end
function terminal_action!(::Val{:input},resources,command)
    write_terminal!(resources,command.fence,command.session_id,command.writer_id,command.input_sequence,command.bytes)
    terminal_action_report(resources,command,command.session_id)
end
function terminal_action!(::Val{:resize},resources,command)
    resize_terminal!(resources,command.fence,command.session_id,command.writer_id,command.columns,command.rows)
    terminal_action_report(resources,command,command.session_id)
end
function terminal_action!(::Val{:disconnect},resources,command)
    disconnect_terminal!(resources,command.fence,command.session_id,command.writer_id)
    terminal_action_report(resources,command,command.session_id)
end
function terminal_action!(::Val{:keepalive},resources,command)
    keepalive_terminal!(resources,command.fence,command.session_id,command.writer_id)
    terminal_action_report(resources,command,command.session_id)
end
function terminal_action!(::Val{:stop},resources,command)
    stop_terminal!(resources,command.fence,command.session_id,command.writer_id)
    terminal_action_report(resources,command,command.session_id)
end
function terminal_action!(::Val{:restart},resources,command)
    id=restart_terminal!(resources,command.fence,command.session_id,command.writer_id;columns=command.columns,rows=command.rows)
    terminal_action_report(resources,command,id)
end

function perform_terminal_action(service,command)
    try
        report=terminal_action!(Val(Symbol(command.action)),service.resources,command)
        return Protocol.validate(report)
    catch error
        reason=error isa TerminalFailure ? string(error.code) : error isa AccessDenied ? "authority_rejected" :
            error isa CapacityUnavailable ? "capacity_unavailable" : "action_rejected"
        return terminal_rejection(command,reason)
    end
end

"""Admit at most one action per exact lease; cached retries never reevaluate input."""
function receive_terminal_command!(service::AgentTerminalService,command::Protocol.TerminalCommand)
    Protocol.validate(command)
    return lock(service.lock) do
        service.closed && return terminal_rejection(command,"channel_closed")
        lock(service.resources.lock) do
            terminal_authority(service.resources,command.fence)
        end
        digest=terminal_command_digest(command)
        lane=get(service.lanes,command.fence.lease_id,nothing)
        if lane !== nothing
            lane.fence==command.fence || return terminal_rejection(command,"fence_changed")
            if command.revision==lane.revision
                (command.request_id==lane.request_id && digest==lane.digest) || return terminal_rejection(command,"retry_changed")
                return lane.report # nothing while the original task remains active
            end
            command.revision>lane.revision || return terminal_rejection(command,"stale_revision")
            command.request_id!=lane.request_id || return terminal_rejection(command,"request_reused")
            lane.task === nothing || return terminal_rejection(command,"action_pending")
        else
            length(service.lanes)<service.ledger.capacity || return terminal_rejection(command,"capacity_unavailable")
        end
        lane=AgentTerminalLane(command.fence,command.revision,command.request_id,digest,nothing,nothing)
        service.lanes[command.fence.lease_id]=lane
        lane.task=@async perform_terminal_action(service,command)
        return nothing
    end
end

function connect_terminals!(service::AgentTerminalService)
    lock(service.lock) do
        service.closed && return
        service.connection === nothing || return
        service.connector !== nothing && !istaskdone(service.connector) && return
        service.ledger.clock()>=service.retry_at || return
        service.retry_at=service.ledger.clock()+2
        service.connector=@async begin
            connection=nothing
            try
                connection=BrokerTerminal(service.endpoint,WorkerIdentity(service.ledger.worker_id);capacity=service.ledger.capacity)
                lock(service.lock) do
                    if !service.closed
                        service.connection=connection;connection=nothing
                    end
                end
            catch
                terminal_service_state!(service,:offline)
            finally
                connection === nothing || close(connection)
            end
        end
    end
end

function terminal_service_state!(service,state)
    previous=lock(service.lock) do
        previous=service.state;service.state=state;previous
    end
    if previous==:online && state!=:online
        # Losing this private channel starts grace even if control heartbeats
        # are still healthy. Output/status traffic cannot keep a writer alive.
        lock(service.resources.lock) do
            for session in values(service.resources.handles)
                session.disconnected_at === nothing && (session.disconnected_at=service.ledger.clock())
            end
        end
    end
    return nothing
end

function completed_terminal_reports!(service)
    reports=Protocol.TerminalReport[]
    lock(service.lock) do
        for lane in values(service.lanes)
            task=lane.task
            if task !== nothing && istaskdone(task)
                report=fetch(task)
                lane.task=nothing;lane.report=report
                push!(reports,report)
            end
        end
    end
    return reports
end

function terminal_active_fences(service)
    lock(service.ledger.lock) do
        [lease.fence for lease in values(service.ledger.leases)
            if haskey(service.resources.profiles.definitions,lease.fence.profile_id) &&
                agent_lease_usable(service.ledger,lease.fence)]
    end
end

function tick_terminals!(service::AgentTerminalService)
    connection=lock(()->service.connection,service.lock)
    connection === nothing && return
    connected=!connection.closed && NATS.status(connection.connection)==NATS.CONNECTED
    terminal_service_state!(service,connected ? :online : :offline)
    if !connected
        lock(service.lock) do
            service.connection === connection && (service.connection=nothing)
        end
        # Retire this connection's queued frames. A replacement channel requires
        # explicit writer reconnect, never replay of old incoming input.
        close(connection)
        return nothing
    end
    fences=terminal_active_fences(service)
    if connected
        for fence in fences;watch_terminal!(connection,fence);end
    end
    for entry in copy(connection.subscriptions)
        entry.fence in fences || unwatch_terminal!(connection,entry.fence)
    end
    # Retire completed action closures promptly: a task may have captured input,
    # whereas the retained retry state contains only its digest and bounded reply.
    for report in completed_terminal_reports!(service)
        if report.fence in fences && connected
            try send_terminal!(connection,report) catch error;error isa BrokerUnavailable || rethrow() end
        end
    end
    lock(service.lock) do
        for (id,lane) in collect(service.lanes)
            lane.fence in fences || lane.task !== nothing || delete!(service.lanes,id)
        end
    end
    for envelope in poll_terminal!(connection)
        report=try receive_terminal_command!(service,envelope.record) catch error
            error isa AccessDenied || rethrow()
            terminal_rejection(envelope.record,"authority_rejected")
        end
        report === nothing || send_terminal!(connection,report)
    end
    return nothing
end

function start_terminals!(service::AgentTerminalService)
    lock(service.lock) do
        service.closed && throw(ArgumentError("terminal service is closed"))
        service.task === nothing || return service
        service.task=@async begin
            while !service.closed
                try connect_terminals!(service);tick_terminals!(service) catch;terminal_service_state!(service,:offline) end
                sleep(0.02)
            end
        end
    end
    return service
end
function stop_terminals!(service::AgentTerminalService)
    lock(()->service.closed=true,service.lock)
    service.task === nothing || wait(service.task;throw=false)
    service.connector === nothing || wait(service.connector;throw=false)
    terminal_service_state!(service,:stopped)
    return nothing
end
function Base.close(service::AgentTerminalService)
    stop_terminals!(service)
    for lane in values(service.lanes)
        lane.task === nothing || wait(lane.task)
    end
    empty!(service.lanes)
    connection=lock(service.lock) do
        previous=service.connection;service.connection=nothing;previous
    end
    connection === nothing || close(connection)
    return nothing
end
