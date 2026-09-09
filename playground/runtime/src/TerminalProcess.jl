"""Report a fixed terminal transport failure without commands, input or output."""
struct TerminalFailure <: Exception
    "Non-sensitive transport failure code."
    code::Symbol
end
Base.showerror(io::IO, error::TerminalFailure) = print(io, "Terminal transport failed: ", error.code)

"""
    TerminalIOLimits(; chunk_bytes=8192, input_bytes=65536,
        output_bytes=1048576, output_bytes_per_second=1048576,
        write_seconds=2, lifetime_seconds=3600, cleanup_seconds=2)

Bound one terminal's local byte transport and attached command lifetime. Output
uses a fixed-size ring and a token-bucket rate limit with a one-second burst.
These bounds do not replace container CPU, memory, PID or scratch requirements.
"""
struct TerminalIOLimits
    "Maximum accepted input or returned output chunk in bytes."
    chunk_bytes::Int
    "Maximum pending input bytes, including a partially written chunk."
    input_bytes::Int
    "Maximum retained output bytes."
    output_bytes::Int
    "Maximum sustained output bytes per second and initial burst bytes."
    output_bytes_per_second::Int
    "Maximum time to write one queued chunk in seconds."
    write_seconds::Float64
    "Maximum attached command lifetime in seconds."
    lifetime_seconds::Float64
    "Each signal or task-join deadline in seconds."
    cleanup_seconds::Float64
end

function TerminalIOLimits(; chunk_bytes=8192, input_bytes=65536, output_bytes=1024^2,
        output_bytes_per_second=1024^2, write_seconds=2, lifetime_seconds=3600, cleanup_seconds=2)
    all(v->v isa Integer && !(v isa Bool), (chunk_bytes,input_bytes,output_bytes,output_bytes_per_second)) ||
        throw(ArgumentError("terminal byte limits require integers"))
    1 <= chunk_bytes <= 16384 && chunk_bytes <= input_bytes <= 262144 &&
        chunk_bytes <= output_bytes <= 4*1024^2 && chunk_bytes <= output_bytes_per_second <= 4*1024^2 ||
        throw(ArgumentError("terminal byte limits are outside their bounds"))
    all(v->v isa Real && !(v isa Bool) && isfinite(v) && 0 < v <= 30, (write_seconds,cleanup_seconds)) &&
        lifetime_seconds isa Real && !(lifetime_seconds isa Bool) && isfinite(lifetime_seconds) &&
        0 < lifetime_seconds <= 86400 || throw(ArgumentError("terminal deadlines are outside their bounds"))
    return TerminalIOLimits(chunk_bytes,input_bytes,output_bytes,output_bytes_per_second,
        Float64(write_seconds),Float64(lifetime_seconds),Float64(cleanup_seconds))
end

"""Retain terminal output in constant storage, addressed by byte position."""
mutable struct TerminalBuffer
    "Fixed-capacity private byte ring."
    data::Vector{UInt8}
    "Total appended byte count; also the exclusive upper cursor."
    sequence::Int
    "Earliest retained byte cursor."
    first::Int
end
TerminalBuffer(capacity::Int) = TerminalBuffer(zeros(UInt8,capacity),0,0)
Base.show(io::IO, buffer::TerminalBuffer) = print(io,"TerminalBuffer(",length(buffer.data)," bytes, <private>)")

function append_terminal!(buffer::TerminalBuffer, bytes::AbstractVector{UInt8})
    length(bytes) <= 9007199254740991-buffer.sequence || throw(TerminalFailure(:sequence_limit))
    capacity=length(buffer.data)
    for i in max(1,length(bytes)-capacity+1):length(bytes)
        buffer.data[mod(buffer.sequence+i-1,capacity)+1]=bytes[i]
    end
    buffer.sequence+=length(bytes)
    buffer.first=max(buffer.first,buffer.sequence-capacity)
    return nothing
end

"""Copy one bounded private output chunk with its cursor and explicit gap flag."""
struct TerminalRead
    "Cursor immediately after the returned bytes."
    cursor::Int
    "Latest produced byte cursor; more data may remain."
    sequence::Int
    "Whether the requested earlier bytes are no longer retained."
    gap::Bool
    "Raw output bytes; UTF-8 or escape sequences may span chunks."
    bytes::Vector{UInt8}
end
Base.show(io::IO, value::TerminalRead) = print(io,"TerminalRead(",length(value.bytes)," bytes, gap=",value.gap,", <private>)")

function read_terminal(buffer::TerminalBuffer, after::Integer, maximum::Int)
    !(after isa Bool) && 0 <= after <= buffer.sequence && maximum > 0 ||
        throw(ArgumentError("invalid terminal output cursor or bound"))
    cursor=max(Int(after),buffer.first)
    count=min(maximum,buffer.sequence-cursor)
    bytes=Vector{UInt8}(undef,count)
    for i in 1:count
        bytes[i]=buffer.data[mod(cursor+i-1,length(buffer.data))+1]
    end
    return TerminalRead(cursor+count,buffer.sequence,after<buffer.first,bytes)
end

"""
    TerminalProcess(; limits=TerminalIOLimits())

Reserve passive transport state for one operator-built command. Construction
does not open a PTY or start a process. No shell, Julia parser, broker, container permission or
lease is supplied here. Production callers must use the approved container
driver; this transport alone is not a sandbox and does not establish REPL readiness.

Input is admitted without waiting for the child. A nonblocking descriptor and
bounded I/O pump keep floods or a non-reading command from blocking the scheduler.
Natural exit, transport failure and the lifetime deadline retire the original
process handle and descriptor. Container deletion remains the physical owner's
responsibility, even when its attached CLI has already exited.
"""
mutable struct TerminalProcess
    "Fresh local identity, not a reusable PID."
    id::UUID
    "Finite transport requirements."
    limits::TerminalIOLimits
    "Original attached command process; absent before spawn."
    process::Union{Nothing,Base.Process}
    "Owned nonblocking PTY master descriptor; -1 after confirmed closure."
    master::Cint
    "Owned PTY pump task."
    pump::Union{Nothing,Task}
    "Independent lifetime/cleanup task."
    lifetime::Union{Nothing,Task}
    "Bounded pending input chunks, never logged or persisted."
    input::Vector{Vector{UInt8}}
    "Written byte count within the first queued chunk."
    input_offset::Int
    "Total queued bytes not yet written."
    pending_bytes::Int
    "Monotonic start of the current input write, or zero."
    writing_since::Float64
    "Bounded output ring."
    output::TerminalBuffer
    "Token-bucket balance for output admission."
    output_tokens::Float64
    "Monotonic time of the previous output-token update."
    token_time::Float64
    "Monotonic start of the original command."
    started_at::Float64
    "Fixed failure code, without raw exception context."
    failure::Union{Nothing,Symbol}
    "Whether input and new I/O have been permanently closed."
    closing::Bool
    "Whether process, PTY and pump retirement have all completed."
    cleanup_complete::Bool
    "Short state and byte-ring updates only."
    lock::ReentrantLock
    "Serialize physical transport retirement and retry."
    cleanup_lock::ReentrantLock
end
Base.show(io::IO, terminal::TerminalProcess) = print(io,"TerminalProcess(",terminal.id,", <private>)")
Base.show(io::IO, ::MIME"text/plain", terminal::TerminalProcess) = show(io,terminal)
terminal_clock() = time_ns()/1e9

function terminal_size(columns,rows)
    all(v->v isa Integer && !(v isa Bool) && 1 <= v <= 1000,(columns,rows)) ||
        throw(ArgumentError("terminal dimensions must be integers in 1:1000"))
    return UInt16[rows,columns,0,0]
end

function terminal_pty()
    Sys.islinux() || throw(TerminalFailure(:linux_required))
    master=Cint(-1); slave=Cint(-1); failure=nothing
    try
        master=ccall(:posix_openpt,Cint,(Cint,),Base.JL_O_RDWR|Base.JL_O_NOCTTY|Base.JL_O_NONBLOCK|Base.JL_O_CLOEXEC)
        master >= 0 || throw(TerminalFailure(:pty_unavailable))
        ccall(:grantpt,Cint,(Cint,),master)==0 && ccall(:unlockpt,Cint,(Cint,),master)==0 ||
            throw(TerminalFailure(:pty_unavailable))
        name=zeros(UInt8,256)
        ccall(:ptsname_r,Cint,(Cint,Ptr{UInt8},Csize_t),master,name,length(name))==0 ||
            throw(TerminalFailure(:pty_unavailable))
        slave=ccall(:open,Cint,(Ptr{UInt8},Cint),name,Base.JL_O_RDWR|Base.JL_O_NOCTTY|Base.JL_O_CLOEXEC)
        slave >= 0 || throw(TerminalFailure(:pty_unavailable))
    catch
        failure=:pty_unavailable
    end
    if failure !== nothing
        master < 0 || ccall(:close,Cint,(Cint,),master)
        slave < 0 || ccall(:close,Cint,(Cint,),slave)
        throw(TerminalFailure(failure))
    end
    return master,slave
end

function TerminalProcess(;limits::TerminalIOLimits=TerminalIOLimits())
    now=terminal_clock()
    return TerminalProcess(uuid4(),limits,nothing,-1,nothing,nothing,Vector{UInt8}[],0,0,0,
        TerminalBuffer(limits.output_bytes),Float64(limits.output_bytes_per_second),now,0,
        nothing,false,false,ReentrantLock(),ReentrantLock())
end

"""
    start_terminal!(terminal, command; columns=100, rows=30) -> TerminalProcess

Start the exact operator-built command on a new Linux PTY, using only its explicit
environment. The caller must retain `terminal` before this acquisition. Failed
partial cleanup therefore remains owned and retryable through `close`. A handle
is single-use; restart requires a new handle and fresh higher-level authority.
"""
function start_terminal!(terminal::TerminalProcess,command::Cmd;columns=100,rows=30)
    size=terminal_size(columns,rows)
    command.env === nothing && throw(ArgumentError("terminal command requires an explicit environment"))
    lock(terminal.cleanup_lock) do
        lock(()->terminal.closing,terminal.lock) && throw(TerminalFailure(:closed))
        terminal.started_at==0 || throw(TerminalFailure(:already_started))
        terminal.started_at=terminal_clock()
        slave=Cint(-1); failure=nothing
        try
            terminal.master,slave=terminal_pty()
            ccall(:ioctl,Cint,(Cint,Culong,Ptr{UInt16}),terminal.master,0x5414,size)==0 ||
                throw(TerminalFailure(:resize_failed)) # Linux TIOCSWINSZ.
            fd=RawFD(slave)
            terminal.process=run(pipeline(ignorestatus(Cmd(command;detach=true));stdin=fd,stdout=fd,stderr=fd);wait=false)
            terminal.pump=@async pump_terminal!(terminal)
            terminal.lifetime=@async terminal_lifetime!(terminal)
        catch error
            failure=error isa TerminalFailure ? error.code : :spawn_failed
        finally
            slave < 0 || ccall(:close,Cint,(Cint,),slave)
        end
        if failure !== nothing
            lock(()->terminal.closing=true,terminal.lock)
            retired=try retire_terminal!(terminal) catch; false end
            retired || (failure=:cleanup_unresolved)
            terminal.failure=failure
            throw(TerminalFailure(failure))
        end
    end
    return terminal
end

function fail_terminal!(terminal,code)
    lock(terminal.lock) do
        terminal.failure === nothing && (terminal.failure=code)
    end
    return nothing
end

function pump_terminal!(terminal)
    block=Vector{UInt8}(undef,terminal.limits.chunk_bytes)
    try
        while !lock(()->terminal.closing,terminal.lock)
            count=ccall(:read,Int,(Cint,Ptr{UInt8},Csize_t),terminal.master,block,length(block))
            read_errno=count<0 ? Base.Libc.errno() : 0
            if count>0
                lock(terminal.lock) do
                    now=terminal_clock(); rate=terminal.limits.output_bytes_per_second
                    terminal.output_tokens=min(Float64(rate),terminal.output_tokens+(now-terminal.token_time)*rate)
                    terminal.token_time=now
                    count <= terminal.output_tokens || throw(TerminalFailure(:output_rate_limit))
                    terminal.output_tokens-=count
                    append_terminal!(terminal.output,@view block[1:count])
                end
            elseif count==0 || read_errno==Base.Libc.EIO
                break # Linux PTY reports EIO when its final slave closes.
            elseif !(read_errno in (Base.Libc.EAGAIN,Base.Libc.EINTR))
                throw(TerminalFailure(:read_failed))
            end
            bytes,offset=lock(terminal.lock) do
                isempty(terminal.input) && return (nothing,0)
                terminal.writing_since==0 && (terminal.writing_since=terminal_clock())
                (first(terminal.input),terminal.input_offset)
            end
            wrote=false
            if bytes !== nothing
                written=GC.@preserve bytes ccall(:write,Int,(Cint,Ptr{UInt8},Csize_t),
                    terminal.master,pointer(bytes,offset+1),length(bytes)-offset)
                write_errno=written<0 ? Base.Libc.errno() : 0
                if written>0
                    wrote=true
                    lock(terminal.lock) do
                        terminal.input_offset+=written; terminal.pending_bytes-=written
                        if terminal.input_offset==length(bytes)
                            fill!(popfirst!(terminal.input),0)
                            terminal.input_offset=0; terminal.writing_since=0
                        end
                    end
                elseif written<0 && !(write_errno in (Base.Libc.EAGAIN,Base.Libc.EINTR))
                    throw(TerminalFailure(:write_failed))
                end
            end
            if count<=0 && !wrote
                FileWatching.poll_fd(RawFD(terminal.master),0.02;readable=true)
            end
            yield() # Continuous output must give control/lifetime tasks a turn.
        end
    catch error
        lock(()->terminal.closing,terminal.lock) || fail_terminal!(terminal,error isa TerminalFailure ? error.code : :stream_failed)
    end
    return nothing
end

function retire_terminal!(terminal)
    lock(terminal.cleanup_lock) do
        terminal.cleanup_complete && return true
        lock(()->terminal.closing=true,terminal.lock)
        process=terminal.process
        if process !== nothing
            for signal in (Base.SIGTERM,Base.SIGKILL)
                process_running(process) || break
                try kill(process,signal) catch end
                timedwait(()->!process_running(process),terminal.limits.cleanup_seconds;pollint=0.01)
            end
            process_running(process) && return false
            wait(process,false)
        end
        pump=terminal.pump
        if pump !== nothing
            timedwait(()->istaskdone(pump),terminal.limits.cleanup_seconds;pollint=0.01)
            istaskdone(pump) || return false
        end
        # poll_fd's watcher is joined before the descriptor can be reused.
        if terminal.master>=0
            fd=terminal.master; terminal.master=-1
            ccall(:close,Cint,(Cint,),fd) # Linux releases fd even on EINTR; never retry that number.
        end
        lock(terminal.lock) do
            foreach(bytes->fill!(bytes,0),terminal.input);empty!(terminal.input)
            terminal.pending_bytes=0;terminal.input_offset=0;terminal.writing_since=0
            terminal.cleanup_complete=true
        end
        return true
    end
end

function terminal_lifetime!(terminal)
    try
        while true
            closing,failure,writing_since=lock(()->(terminal.closing,terminal.failure,terminal.writing_since),terminal.lock)
            closing && break
            if !process_running(terminal.process)
                # Drain the finite trailing output before retiring a normal exit.
                timedwait(()->istaskdone(terminal.pump),terminal.limits.cleanup_seconds;pollint=0.01)
                break
            end
            failure === nothing || break
            now=terminal_clock()
            if now-terminal.started_at >= terminal.limits.lifetime_seconds
                fail_terminal!(terminal,:lifetime_limit);break
            end
            if writing_since>0 && now-writing_since>=terminal.limits.write_seconds
                fail_terminal!(terminal,:input_stalled);break
            end
            if istaskdone(terminal.pump)
                # Process-exit notification can arrive just after PTY EIO.
                timedwait(()->!process_running(terminal.process),0.1;pollint=0.01)
                process_running(terminal.process) && fail_terminal!(terminal,:stream_closed)
                break
            end
            sleep(0.01)
        end
    catch
        fail_terminal!(terminal,:lifetime_failed)
    end
    retired=try retire_terminal!(terminal) catch; false end
    retired || fail_terminal!(terminal,:cleanup_unresolved)
    return nothing
end

"""
    write_terminal!(terminal, bytes) -> Nothing

Admit a bounded copy of input without waiting for the process. Full admission
fails explicitly and does not append a partial chunk. This is local queue
acceptance, not an acknowledgement that Julia evaluated the input. The higher
session owner must enforce the sole writer and input sequence; never auto-replay.
"""
function write_terminal!(terminal::TerminalProcess,bytes::AbstractVector{UInt8})
    0 < length(bytes) <= terminal.limits.chunk_bytes || throw(ArgumentError("terminal input chunk exceeds its bound"))
    lock(terminal.lock) do
        (terminal.process === nothing || terminal.closing || terminal.failure !== nothing) && throw(TerminalFailure(:closed))
        terminal.pending_bytes+length(bytes) <= terminal.limits.input_bytes || throw(TerminalFailure(:input_full))
        push!(terminal.input,Vector{UInt8}(bytes));terminal.pending_bytes+=length(bytes)
    end
    return nothing
end

"""
    read_terminal(terminal, after) -> TerminalRead

Return a bounded copy of available output after the supplied byte cursor.
A lagging cursor receives an explicit gap; this read does not consume history.
"""
read_terminal(terminal::TerminalProcess,after::Integer) = lock(()->
    read_terminal(terminal.output,after,terminal.limits.chunk_bytes),terminal.lock)

"""
    resize_terminal!(terminal, columns, rows) -> Nothing

Resize the owned PTY and notify only its original attached command process.
Dimensions are terminal cells, each in 1:1000; no container command is constructed.
"""
function resize_terminal!(terminal::TerminalProcess,columns,rows)
    size=terminal_size(columns,rows)
    lock(terminal.lock) do
        (terminal.process === nothing || terminal.closing) && throw(TerminalFailure(:closed))
        ccall(:ioctl,Cint,(Cint,Culong,Ptr{UInt16}),terminal.master,0x5414,size)==0 || throw(TerminalFailure(:resize_failed))
        process_running(terminal.process) || throw(TerminalFailure(:closed))
        signaled=try kill(terminal.process,28);true catch;false end # Linux SIGWINCH.
        signaled || throw(TerminalFailure(:resize_failed))
    end
    return nothing
end

"""
    close(terminal::TerminalProcess) -> Nothing

Close admission and join exact command/PTY ownership. Repeated close retries
unresolved retirement; a failure retains the handle and throws `TerminalFailure`.
This does not remove a container or release a higher-level application lease.
"""
function Base.close(terminal::TerminalProcess)
    lock(()->terminal.closing=true,terminal.lock)
    if terminal.lifetime !== nothing
        timedwait(()->istaskdone(terminal.lifetime),4*terminal.limits.cleanup_seconds+1;pollint=0.01)==:ok ||
            throw(TerminalFailure(:cleanup_unresolved))
    end
    retired=try retire_terminal!(terminal) catch; false end
    retired || throw(TerminalFailure(:cleanup_unresolved))
    return nothing
end
