using .Protocol: AssignmentFence

"""Retain an exact terminal-assignment subscription and its original NATS handle."""
struct TerminalSubscription
    "Full assignment identity; the subject alone is insufficient authority."
    fence::AssignmentFence
    "Original bounded subscription, never a wildcard."
    subscription::NATS.Sub
end

"""
    BrokerTerminal(endpoint, identity; worker_ids=(), capacity=256)

Own a separate authenticated transient terminal connection. No subscriptions or
resources are acquired until explicitly watched. Each watched assignment has an
exact worker-boot/lease/generation subject and a bounded queue. Round-robin polling
does not execute callbacks, and disconnected input publication is never queued
for replay. The caller remains responsible for live lease and principal checks.
"""
mutable struct BrokerTerminal{I<:AbstractBrokerIdentity}
    "Server-provisioned broker identity, never a browser credential."
    identity::I
    "Separate private byte/control connection."
    connection::NATS.Connection
    "Operator-approved worker identities for a coordinator connection."
    worker_ids::Set{String}
    "Maximum simultaneously watched assignments."
    capacity::Int
    "Exact subscriptions, retaining the complete fence."
    subscriptions::Vector{TerminalSubscription}
    "Next queue in a fair bounded poll."
    cursor::Int
    "Saturating malformed-frame count, without raw payloads."
    rejected::Int
    "Permanently reject new use after close."
    closed::Bool
    "Serialize subscription ownership, bounded polling and teardown."
    lock::ReentrantLock
end
Base.show(io::IO,::BrokerTerminal)=print(io,"BrokerTerminal(<private>)")

function BrokerTerminal(endpoint::BrokerEndpoint,identity::AbstractBrokerIdentity;worker_ids=(),capacity=256)
    capacity isa Integer && !(capacity isa Bool) && 1<=capacity<=256 || throw(ArgumentError("invalid terminal channel capacity"))
    ids=Set{String}(Protocol.runtime_token.(collect(worker_ids)))
    if identity isa WorkerIdentity
        isempty(ids) || ids==Set([identity.worker_id]) || throw(ArgumentError("worker terminal channel cannot watch another worker"))
        push!(ids,identity.worker_id)
    else
        identity isa CoordinatorIdentity && 1<=length(ids)<=128 || throw(ArgumentError("terminal channel requires approved worker identities"))
    end
    connection=connect_broker(endpoint,identity)
    BrokerTerminal(identity,connection,ids,Int(capacity),TerminalSubscription[],1,0,false,ReentrantLock())
end

terminal_input_kind(::CoordinatorIdentity)=Protocol.TerminalReport
terminal_input_kind(::WorkerIdentity)=Protocol.TerminalCommand
terminal_input_direction(::CoordinatorIdentity)=:report
terminal_input_direction(::WorkerIdentity)=:command
terminal_output_direction(::CoordinatorIdentity,::Protocol.TerminalCommand)=:command
terminal_output_direction(::WorkerIdentity,::Protocol.TerminalReport)=:report
terminal_output_direction(::AbstractBrokerIdentity,::Protocol.RuntimeRecord)=throw(AccessDenied(403,"Terminal record is not permitted for this broker role"))

function terminal_subscription(channel,fence)
    findfirst(entry->entry.fence==fence,channel.subscriptions)
end

"""Watch one exact assignment without granting terminal or writer authority."""
function watch_terminal!(channel::BrokerTerminal,fence::AssignmentFence)
    Protocol.validate(fence)
    fence.worker_id in channel.worker_ids || throw(AccessDenied(403,"Terminal worker is not provisioned"))
    return lock(channel.lock) do
        !channel.closed && NATS.status(channel.connection)==NATS.CONNECTED || throw(BrokerUnavailable())
        terminal_subscription(channel,fence) === nothing || return nothing
        subject=Protocol.terminal_subject(fence,terminal_input_direction(channel.identity))
        # A reused subject with different owner/run/profile metadata is invalid.
        all(entry->entry.subscription.subject!=subject,channel.subscriptions) || throw(AccessDenied(409,"Terminal subscription fence differs"))
        length(channel.subscriptions)<channel.capacity || throw(CapacityUnavailable())
        subscription=nothing
        try
            Logging.with_logger(Logging.NullLogger()) do
                subscription=NATS.subscribe(channel.connection,subject;channel_size=4)
                NATS.ping(channel.connection;timeout=2,measure=false)
            end
            push!(channel.subscriptions,TerminalSubscription(fence,subscription))
        catch
            if subscription !== nothing
                Logging.with_logger(Logging.NullLogger()) do
                    try NATS.unsubscribe(channel.connection,subscription) catch end
                end
            end
            throw(BrokerUnavailable())
        end
        return nothing
    end
end

"""Remove only the original exact terminal subscription; this does not release a lease."""
function unwatch_terminal!(channel::BrokerTerminal,fence::AssignmentFence)
    lock(channel.lock) do
        index=terminal_subscription(channel,fence)
        index === nothing && return nothing
        entry=channel.subscriptions[index]
        Logging.with_logger(Logging.NullLogger()) do
            NATS.unsubscribe(channel.connection,entry.subscription)
        end
        deleteat!(channel.subscriptions,index)
    end
    return nothing
end

"""Publish one bounded role-appropriate terminal frame, with no disconnected replay."""
function send_terminal!(channel::BrokerTerminal,record::Union{Protocol.TerminalCommand,Protocol.TerminalReport})
    Protocol.validate(record)
    direction=terminal_output_direction(channel.identity,record)
    payload=Protocol.encode_message(record)
    ncodeunits(payload)<=Protocol.MAX_TERMINAL_FRAME_BYTES || throw(ArgumentError("terminal frame exceeds its bound"))
    lock(channel.lock) do
        !channel.closed && NATS.status(channel.connection)==NATS.CONNECTED || throw(BrokerUnavailable())
        terminal_subscription(channel,record.fence) !== nothing || throw(AccessDenied(403,"Terminal assignment is not watched"))
        try
            NATS.publish(channel.connection,Protocol.terminal_subject(record.fence,direction),payload)
        catch
            throw(BrokerUnavailable())
        end
    end
    return nothing
end

"""Poll bounded exact-assignment terminal frames fairly, without evaluating input."""
function poll_terminal!(channel::BrokerTerminal;limit::Integer=32)
    !(limit isa Bool) && 1<=limit<=128 || throw(ArgumentError("invalid terminal poll bound"))
    return lock(channel.lock) do
        result=ControlEnvelope[]
        (channel.closed || isempty(channel.subscriptions)) && return result
        queues=length(channel.subscriptions);empty_checks=0;consumed=0
        kind=terminal_input_kind(channel.identity)
        while empty_checks<queues && consumed<limit
            index=mod1(channel.cursor,queues);channel.cursor=mod1(index+1,queues)
            entry=channel.subscriptions[index]
            message=NATS.next(channel.connection,entry.subscription;no_wait=true,no_throw=true)
            if message === nothing
                empty_checks+=1;continue
            end
            empty_checks=0;consumed+=1
            try
                message.subject==entry.subscription.subject && message.reply_to === nothing &&
                    length(message.payload)<=Protocol.MAX_TERMINAL_FRAME_BYTES || throw(ArgumentError("invalid terminal frame"))
                record=Protocol.decode_terminal_message(kind,NATS.payload(message))
                record.fence==entry.fence || throw(ArgumentError("terminal assignment fence differs"))
                push!(result,ControlEnvelope(entry.fence.worker_id,record))
            catch
                channel.rejected==typemax(Int) || (channel.rejected+=1)
            end
        end
        return result
    end
end

function Base.close(channel::BrokerTerminal)
    lock(channel.lock) do
        channel.closed && return nothing
        channel.closed=true
        Logging.with_logger(Logging.NullLogger()) do
            for entry in channel.subscriptions
                try NATS.unsubscribe(channel.connection,entry.subscription) catch end
            end
            NATS.drain(channel.connection)
        end
        empty!(channel.subscriptions)
    end
    return nothing
end
