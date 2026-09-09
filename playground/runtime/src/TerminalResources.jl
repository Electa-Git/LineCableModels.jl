"""Configure finite startup, disconnect and idle deadlines for private sessions."""
struct TerminalSessionLimits
    "Maximum guarded REPL startup time in seconds."
    startup_seconds::Float64
    "Maximum disconnected writer grace in seconds."
    disconnect_seconds::Float64
    "Maximum interval without fresh writer presence, in seconds."
    presence_seconds::Float64
    "Maximum interval without explicit writer activity in seconds."
    idle_seconds::Float64
    "Maximum wait per session-task cleanup attempt in seconds."
    cleanup_seconds::Float64
    function TerminalSessionLimits(;startup_seconds=120,disconnect_seconds=30,presence_seconds=15,idle_seconds=1800,cleanup_seconds=15)
        all(v->v isa Real && !(v isa Bool) && isfinite(v) && 0<v<=3600,
            (startup_seconds,disconnect_seconds,presence_seconds,idle_seconds,cleanup_seconds)) ||
            throw(ArgumentError("terminal session deadlines must be in (0, 3600] seconds"))
        new(startup_seconds,disconnect_seconds,presence_seconds,idle_seconds,cleanup_seconds)
    end
end

"""Retain one private stream's identity and lifecycle, including failed cleanup."""
mutable struct TerminalSession
    "Exact live assignment; never changed during a session."
    fence::AssignmentFence
    "Fresh stream identity, not a PID or lease identifier."
    id::String
    "Sole writer connection identity; never reassigned to another connection."
    writer::String
    "Original driver-owned PTY, absent before acquisition."
    process::Union{Nothing,TerminalProcess}
    "Independent startup, deadline and physical-retirement task."
    task::Union{Nothing,Task}
    "Starting, ready, exited, failed or closing."
    phase::Symbol
    "Latest explicit writer activity on the agent's monotonic clock."
    activity_at::Float64
    "Latest fresh writer presence; heartbeat alone does not extend idle life."
    writer_seen_at::Float64
    "Original disconnect time; repeated disconnects do not extend grace."
    disconnected_at::Union{Nothing,Float64}
    "Latest locally admitted input sequence, not an evaluation acknowledgement."
    input_sequence::Int
    "Digest of only the last admitted chunk; raw input is not retained here."
    input_digest::Union{Nothing,Vector{UInt8}}
    "Private startup prefix length omitted from the exposed output stream."
    output_origin::Int
    "Non-sensitive failure code, never raw Julia output or exception context."
    failure::Union{Nothing,Symbol}
    "Permanently reject new input or implicit restart."
    closing::Bool
    "Whether exact physical retirement has been confirmed."
    cleanup_complete::Bool
    "Serialize retries of physical retirement for this session only."
    cleanup_lock::ReentrantLock
end
Base.show(io::IO,session::TerminalSession)=print(io,"TerminalSession(",session.id,", <private>)")

"""
    TerminalResources(driver; limits=TerminalSessionLimits())

Own private terminal sessions under the agent's existing live lease ledger.
Construction and binding are passive. Explicit open starts an independently
guarded PTY; no input is accepted until its acquisition-bound startup marker and
post-start physical checks pass. Each session has one writer, sequenced bounded
input and bounded output with explicit gaps. A dead session never restarts itself.

This owner does not expose HTTP/NATS routes or replace gateway authentication.
Callers must authorize the principal and run before supplying the exact fence.
"""
mutable struct TerminalResources{D<:AbstractTerminalDriver} <: AbstractAgentResources
    "Approved physical terminal adapter."
    driver::D
    "Verified container-only terminal definitions."
    profiles::ProfileRegistry
    "The same authority used for every kind of agent resource."
    ledger::Union{Nothing,AgentLeaseLedger}
    "Finite session policy independent of scientific preparation."
    limits::TerminalSessionLimits
    "Sessions keyed by lease, including ended sessions awaiting lease release."
    handles::Dict{String,TerminalSession}
    "Whether original physical ownership has been recovered."
    recovered::Bool
    "Whether new actions are permanently forbidden."
    closed::Bool
    "Short authority, sequence and status updates; no physical operations."
    lock::ReentrantLock
    "Serialize recovery and whole-owner closure."
    lifecycle_lock::ReentrantLock
end

function TerminalResources(driver::AbstractTerminalDriver;limits=TerminalSessionLimits())
    RequiredInterfaces.check_interface_implemented(AbstractTerminalDriver,typeof(driver)) === true ||
        throw(ArgumentError("terminal driver does not implement required ownership hooks"))
    profiles=installed_profiles(driver)
    profiles isa ProfileRegistry && all(p->p.kind==:terminal && p.isolation==:container,
        values(profiles.definitions)) || throw(ArgumentError("terminal driver must return verified terminal profiles"))
    TerminalResources(driver,profiles,nothing,limits,Dict{String,TerminalSession}(),false,false,ReentrantLock(),ReentrantLock())
end
installed_profiles(resources::TerminalResources)=resources.profiles
Base.show(io::IO,::TerminalResources)=print(io,"TerminalResources(<private>)")

function bind_agent!(resources::TerminalResources,ledger::AgentLeaseLedger)
    lock(resources.lock) do
        !resources.closed && (resources.ledger === nothing || resources.ledger === ledger) ||
            throw(ArgumentError("terminal resources cannot change agent incarnation"))
        all(get(ledger.profiles.definitions,id,nothing) === profile for (id,profile) in resources.profiles.definitions) ||
            throw(ArgumentError("terminal definitions do not match lease authority"))
        resources.ledger=ledger
    end
    return nothing
end

function recover_owned!(resources::TerminalResources)
    lock(resources.lifecycle_lock) do
        !resources.closed && resources.ledger !== nothing || throw(ArgumentError("terminal resources lack live ownership"))
        resources.recovered && return nothing
        recover_owned!(resources.driver) === nothing || throw(ArgumentError("terminal recovery failed"))
        resources.recovered=true
    end
    return nothing
end

function terminal_authority(resources,fence)
    !resources.closed && resources.recovered && resources.ledger !== nothing &&
        agent_lease_usable(resources.ledger,fence) || throw(AccessDenied(409,"Terminal assignment is not usable"))
    profile=get(resources.profiles.definitions,fence.profile_id,nothing)
    profile !== nothing && (string(profile.version),profile.fingerprint)==(fence.profile_version,fence.fingerprint) ||
        throw(AccessDenied(409,"Terminal profile does not match its assignment"))
    return profile
end

function terminal_session(resources,fence,id)
    terminal_authority(resources,fence)
    session=get(resources.handles,fence.lease_id,nothing)
    session !== nothing && session.fence==fence && session.id==id ||
        throw(AccessDenied(409,"Terminal stream identity is not current"))
    return session
end

function terminal_deadline(resources,session)
    now=resources.ledger.clock()
    0<=now-session.activity_at<resources.limits.idle_seconds || return :idle_timeout
    if session.disconnected_at === nothing && now-session.writer_seen_at>=resources.limits.presence_seconds
        session.disconnected_at=session.writer_seen_at+resources.limits.presence_seconds
    end
    disconnected=session.disconnected_at
    disconnected === nothing || 0<=now-disconnected<resources.limits.disconnect_seconds || return :disconnect_timeout
    return nothing
end

function terminal_writer(resources,session,writer;ready=true)
    session.writer==writer || throw(AccessDenied(409,"Terminal already has a different writer"))
    !session.closing && terminal_deadline(resources,session) === nothing || throw(TerminalFailure(:session_closed))
    !ready || (session.phase==:ready && session.disconnected_at === nothing) || throw(TerminalFailure(:not_ready))
    return nothing
end

function terminal_start_authority(resources,session)
    terminal_authority(resources,session.fence)
    session.closing && throw(TerminalFailure(:session_closed))
    reason=terminal_deadline(resources,session)
    reason === nothing || throw(TerminalFailure(reason))
    return nothing
end

"""
    open_terminal!(resources, fence, writer; columns=100, rows=30) -> String

Explicitly admit one terminal and return its fresh stream UUID. Retrying with the
same writer UUID returns the same stream; it never replaces an exited process.
The same writer may reconnect within grace. A different writer is rejected even
while disconnected. Poll terminal_status for startup completion.
"""
function open_terminal!(resources::TerminalResources,fence::AssignmentFence,writer::AbstractString;columns=100,rows=30)
    Protocol.runtime_uuid(writer);terminal_size(columns,rows)
    return lock(resources.lock) do
        profile=terminal_authority(resources,fence)
        found=get(resources.handles,fence.lease_id,nothing)
        if found !== nothing
            found.fence==fence || throw(AccessDenied(409,"Terminal lease identity differs"))
            terminal_writer(resources,found,writer;ready=false)
            found.disconnected_at=nothing;found.activity_at=found.writer_seen_at=resources.ledger.clock()
            return found.id
        end
        length(resources.handles)<resources.ledger.capacity || throw(CapacityUnavailable())
        session=TerminalSession(fence,string(uuid4()),String(writer),nothing,nothing,:starting,
            resources.ledger.clock(),resources.ledger.clock(),nothing,0,nothing,0,nothing,false,false,ReentrantLock())
        resources.handles[fence.lease_id]=session
        session.task=@async run_terminal_session!(resources,session,profile,Int(columns),Int(rows))
        return session.id
    end
end

# Startup output has a small independent bound. No user input is admitted yet,
# so an acquisition marker cannot be forged by submitted REPL expressions.
function terminal_startup_offset(process,marker)
    1<=length(marker)<=128 || throw(TerminalFailure(:startup_marker_invalid))
    return lock(process.lock) do
        process.failure === nothing && !process.closing || throw(TerminalFailure(:startup_failed))
        process.output.first==0 && process.output.sequence<=8192 || throw(TerminalFailure(:startup_output_limit))
        output=read_terminal(process.output,0,8192).bytes
        for i in 1:length(output)-length(marker)+1
            @views output[i:i+length(marker)-1]==marker && return i+length(marker)-1
        end
        return nothing
    end
end

function retire_terminal_session!(resources,session)
    lock(session.cleanup_lock) do
        session.cleanup_complete && return true
        # Revoke byte transport even if physical deletion is temporarily denied.
        # Still attempt container cleanup when local PTY retirement is unresolved.
        local_complete=session.process === nothing || try close(session.process);true catch;false end
        physical_complete=try release_owned!(resources.driver,session.fence) === true catch;false end
        complete=local_complete && physical_complete
        lock(resources.lock) do
            session.cleanup_complete=complete
        end
        return complete
    end
end

function run_terminal_session!(resources,session,profile,columns,rows)
    deadline=terminal_clock()+resources.limits.startup_seconds
    try
        lock(resources.lock) do
            terminal_start_authority(resources,session)
        end
        process=terminal_for!(resources.driver,profile,session.fence)
        lock(resources.lock) do
            session.process=process
            terminal_start_authority(resources,session)
        end
        terminal_clock()<deadline || throw(TerminalFailure(:startup_timeout))
        start_owned_terminal!(resources.driver,profile,session.fence,columns,rows) === process ||
            throw(TerminalFailure(:process_identity_changed))
        marker=terminal_ready_marker(resources.driver,profile,session.fence)
        while true
            lock(resources.lock) do
                terminal_start_authority(resources,session)
            end
            terminal_clock()<deadline || throw(TerminalFailure(:startup_timeout))
            offset=terminal_startup_offset(process,marker)
            if offset !== nothing
                verify_terminal!(resources.driver,profile,session.fence)
                lock(resources.lock) do
                    terminal_start_authority(resources,session)
                    lock(process.lock) do
                        !process.closing && process.failure === nothing && process.process !== nothing &&
                            process_running(process.process) || throw(TerminalFailure(:startup_failed))
                    end
                    session.output_origin=offset;session.phase=:ready
                end
                break
            end
            sleep(0.02)
        end
        while true
            lock(resources.lock) do
                terminal_authority(resources,session.fence)
                reason=terminal_deadline(resources,session)
                reason === nothing || throw(TerminalFailure(reason))
            end
            lock(()->session.closing,resources.lock) && break
            ended,failure=lock(()->(process.closing,process.failure),process.lock)
            failure === nothing || throw(TerminalFailure(failure))
            if ended
                lock(()->session.phase=:exited,resources.lock)
                break
            end
            sleep(0.05)
        end
    catch error
        lock(resources.lock) do
            if !session.closing
                session.phase=:failed
                session.failure=error isa TerminalFailure ? error.code : error isa AccessDenied ? :lease_lost : :startup_failed
            end
        end
    finally
        lock(()->session.closing=true,resources.lock)
        retire_terminal_session!(resources,session)
    end
    return nothing
end

"""Return bounded session diagnostics without bytes, writer IDs or credentials."""
function terminal_status(resources::TerminalResources,fence::AssignmentFence,id::AbstractString)
    return lock(resources.lock) do
        session=terminal_session(resources,fence,id)
        sequence=session.process === nothing ? 0 : lock(()->session.process.output.sequence,session.process.lock)
        (;session_id=session.id,phase=session.phase,writer_connected=session.disconnected_at === nothing,
            input_sequence=session.input_sequence,output_sequence=max(0,sequence-session.output_origin),
            failure=session.failure,cleanup_pending=session.closing && !session.cleanup_complete)
    end
end

"""
    write_terminal!(resources, fence, id, writer, sequence, bytes) -> Int

Admit one bounded input chunk from the sole connected writer. A matching retry of
the last sequence is acknowledged without replay; altered, older or skipped
sequences fail. The acknowledgement means queued bytes, not completed Julia code.
"""
function write_terminal!(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,
        writer::AbstractString,sequence::Integer,bytes::AbstractVector{UInt8})
    !(sequence isa Bool) && 1<=sequence<=9007199254740991 || throw(ArgumentError("invalid terminal input sequence"))
    return lock(resources.lock) do
        session=terminal_session(resources,fence,id)
        terminal_writer(resources,session,writer)
        0<length(bytes)<=session.process.limits.chunk_bytes || throw(ArgumentError("terminal input exceeds its bound"))
        payload=Vector{UInt8}(bytes);digest=SHA.sha256(payload)
        if sequence==session.input_sequence
            digest==session.input_digest || throw(AccessDenied(409,"Terminal input retry differs"))
            return session.input_sequence
        end
        sequence==session.input_sequence+1 || throw(AccessDenied(409,"Terminal input sequence is not next"))
        # Admission is nonblocking. Keep lease revocation atomic with queuing;
        # no other thread may revoke the lease between this check and the write.
        return lock(resources.ledger.lock) do
            terminal_authority(resources,fence)
            write_terminal!(session.process,payload)
            session.input_sequence=Int(sequence);session.input_digest=digest
            session.activity_at=session.writer_seen_at=resources.ledger.clock()
            session.input_sequence
        end
    end
end

"""Read bounded private output with relative byte cursors and explicit gaps."""
function read_terminal(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,after::Integer)
    return lock(resources.lock) do
        session=terminal_session(resources,fence,id)
        session.output_origin>0 && session.process !== nothing || throw(TerminalFailure(:not_ready))
        !(after isa Bool) && 0<=after<=9007199254740991-session.output_origin || throw(ArgumentError("invalid terminal output cursor"))
        output=read_terminal(session.process,Int(after)+session.output_origin)
        TerminalRead(output.cursor-session.output_origin,output.sequence-session.output_origin,output.gap,output.bytes)
    end
end

"""Resize a current writer's PTY in terminal cells, without changing its identity."""
function resize_terminal!(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,
        writer::AbstractString,columns,rows)
    terminal_size(columns,rows)
    lock(resources.lock) do
        session=terminal_session(resources,fence,id);terminal_writer(resources,session,writer)
        resize_terminal!(session.process,columns,rows)
        session.activity_at=session.writer_seen_at=resources.ledger.clock()
    end
    return nothing
end

"""Refresh this connected writer's presence without renewing idle or lease authority."""
function keepalive_terminal!(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,writer::AbstractString)
    lock(resources.lock) do
        session=terminal_session(resources,fence,id)
        terminal_writer(resources,session,writer;ready=false)
        session.disconnected_at === nothing || throw(TerminalFailure(:not_connected))
        session.writer_seen_at=resources.ledger.clock()
    end
    return nothing
end

"""Begin finite reconnect grace once; repeated disconnects cannot extend it."""
function disconnect_terminal!(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,writer::AbstractString)
    lock(resources.lock) do
        session=terminal_session(resources,fence,id)
        session.writer==writer || throw(AccessDenied(409,"Terminal writer differs"))
        session.disconnected_at === nothing && (session.disconnected_at=resources.ledger.clock())
    end
    return nothing
end

"""Stop this writer's exact stream without freeing or implicitly replacing its identity."""
function stop_terminal!(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,writer::AbstractString)
    return lock(resources.lock) do
        session=terminal_session(resources,fence,id)
        session.writer==writer || throw(AccessDenied(409,"Terminal writer differs"))
        session.closing=true
        session.phase in (:failed,:exited) || (session.phase=:closing)
        return session
    end
end

"""
    restart_terminal!(resources, fence, id, writer; columns=100, rows=30) -> String

Explicitly retire the writer's exact stream and open a fresh one only after
confirmed physical cleanup. Lost variables are not restored. Failed retirement
retains the old identity and prevents replacement. The relay must deduplicate
the enclosing request; it may not automatically replay a restart command.
"""
function restart_terminal!(resources::TerminalResources,fence::AssignmentFence,id::AbstractString,
        writer::AbstractString;columns=100,rows=30)
    terminal_size(columns,rows)
    session=stop_terminal!(resources,fence,id,writer)
    task=session.task
    task === nothing || timedwait(()->istaskdone(task),resources.limits.cleanup_seconds;pollint=0.02)==:ok ||
        throw(TerminalFailure(:cleanup_unresolved))
    retire_terminal_session!(resources,session) || throw(TerminalFailure(:cleanup_unresolved))
    return lock(resources.lock) do
        terminal_session(resources,fence,id) === session || throw(TerminalFailure(:session_closed))
        delete!(resources.handles,fence.lease_id)
        open_terminal!(resources,fence,writer;columns,rows)
    end
end

function release_owned!(resources::TerminalResources,fence::AssignmentFence)
    session=lock(resources.lock) do
        found=get(resources.handles,fence.lease_id,nothing)
        if found !== nothing
            found.fence==fence || throw(ArgumentError("terminal cleanup fence differs"))
            found.closing=true
            found.phase in (:exited,:failed) || (found.phase=:closing)
        end
        found
    end
    session === nothing && return release_owned!(resources.driver,fence)
    task=session.task
    task === nothing || timedwait(()->istaskdone(task),resources.limits.cleanup_seconds;pollint=0.02)==:ok || return false
    retire_terminal_session!(resources,session) || return false
    lock(resources.lock) do
        get(resources.handles,fence.lease_id,nothing) === session && delete!(resources.handles,fence.lease_id)
    end
    return true
end

function Base.close(resources::TerminalResources)
    lock(resources.lifecycle_lock) do
        fences=lock(resources.lock) do
            resources.closed=true
            for session in values(resources.handles);session.closing=true;end
            [session.fence for session in values(resources.handles)]
        end
        completed=true
        for fence in fences
            completed &= try release_owned!(resources,fence) catch;false end
        end
        completed || throw(ArgumentError("terminal cleanup remains unresolved"))
        close(resources.driver)
    end
    return nothing
end
