"""
    CommandFailure

Report a fixed failure code for an owned host command. Executable arguments,
environment, output and original exception contexts are deliberately omitted.
"""
struct CommandFailure <: Exception
    "Finite reason, such as deadline, output_limit or cleanup_unresolved."
    code::Symbol
end
Base.showerror(io::IO, error::CommandFailure) = print(io, "Owned host command failed: ", error.code)

"""Hold bounded private command output; ordinary display never prints it."""
struct CommandResult
    "Operating-system exit code."
    exitcode::Int
    "Bounded stdout, for internal inspection only."
    output::String
    "Bounded stderr, never copied into public status."
    diagnostic::String
end
Base.show(io::IO, result::CommandResult) = print(io, "CommandResult(exitcode=", result.exitcode, ", <private>)")
Base.show(io::IO, ::MIME"text/plain", result::CommandResult) = show(io, result)

"""Retain one original host command and every pipe task until teardown completes."""
mutable struct OwnedCommand
    "Fresh local identity, not a stored PID."
    id::UUID
    "Original process handle; absent before successful spawn."
    process::Union{Nothing,Base.Process}
    "Owned stdout and stderr pipes."
    pipes::Tuple{Pipe,Pipe}
    "Bounded reader tasks."
    readers::Vector{Task}
    "Combined byte count across both output streams."
    bytes::Int
    "Fixed reader failure, or nothing."
    failure::Union{Nothing,Symbol}
    "Serialize repeated cleanup of this command only."
    cleanup_lock::ReentrantLock
end
Base.show(io::IO, command::OwnedCommand) = print(io, "OwnedCommand(", command.id, ", <private>)")

"""
    CommandRunner(; capacity=8, timeout_seconds=10, maximum_bytes=1048576,
                    cleanup_seconds=2)

Bound concurrent operator-owned host commands and their combined output. The
constructor is passive. A full runner rejects work immediately; it never queues.
Closed pipes and observed process exit are required before releasing a slot.
An unresolved retirement retains its original handle for a later close retry.

This owns CLI processes, not the containers those commands may create. Container
receipts and exact-ID reconciliation remain the physical driver's responsibility.
"""
mutable struct CommandRunner
    "Maximum commands including unresolved cleanup."
    capacity::Int
    "Default wall-clock deadline in seconds."
    timeout_seconds::Float64
    "Maximum combined stdout/stderr bytes per command."
    maximum_bytes::Int
    "Each bounded signal/reader join interval in seconds."
    cleanup_seconds::Float64
    "Original occupied command handles."
    active::Dict{UUID,OwnedCommand}
    "Serialize short admission and closed-state changes."
    lock::ReentrantLock
    "Serialize whole-runner teardown."
    cleanup_lock::ReentrantLock
    "Whether further commands are forbidden."
    closed::Bool
end

function CommandRunner(; capacity=8, timeout_seconds=10, maximum_bytes=1024^2, cleanup_seconds=2)
    capacity isa Integer && !(capacity isa Bool) && 1 <= capacity <= 32 ||
        throw(ArgumentError("host command capacity must be in 1:32"))
    all(v -> v isa Real && !(v isa Bool) && isfinite(v) && 0 < v <= 120,
        (timeout_seconds, cleanup_seconds)) || throw(ArgumentError("host command time bounds must be in (0, 120] seconds"))
    maximum_bytes isa Integer && !(maximum_bytes isa Bool) && 1 <= maximum_bytes <= 4 * 1024^2 ||
        throw(ArgumentError("host command output limit must be in 1:4194304 bytes"))
    return CommandRunner(capacity, timeout_seconds, maximum_bytes, cleanup_seconds,
        Dict{UUID,OwnedCommand}(), ReentrantLock(), ReentrantLock(), false)
end
Base.show(io::IO, runner::CommandRunner) = print(io, "CommandRunner(capacity=", runner.capacity, ", <owned>)")

function command_reader!(runner::CommandRunner, handle::OwnedCommand, pipe::Pipe)
    output = UInt8[]
    sizehint!(output, min(runner.maximum_bytes, 8192))
    block = Vector{UInt8}(undef, min(runner.maximum_bytes + 1, 8192))
    try
        while !eof(pipe)
            count = readbytes!(pipe, block, length(block))
            handle.bytes += count
            if handle.bytes > runner.maximum_bytes
                handle.failure = :output_limit
                break
            end
            append!(output, @view block[1:count])
            # A continuously readable command must not monopolize control tasks.
            yield()
        end
    catch
        handle.failure === nothing && (handle.failure = :pipe_failed)
    end
    return String(output)
end

function retire_command!(runner::CommandRunner, handle::OwnedCommand)
    return lock(handle.cleanup_lock) do
        process = handle.process
        if process !== nothing
            for signal in (Base.SIGTERM, Base.SIGKILL)
                Base.process_running(process) || break
                try
                    kill(process, signal)
                catch
                    # Actual exit, not signal success, establishes release.
                end
                timedwait(() -> !Base.process_running(process), runner.cleanup_seconds; pollint=0.01)
            end
            Base.process_running(process) && return false
            wait(process, false)
        end
        for pipe in handle.pipes
            try
                close(pipe)
            catch
                return false
            end
        end
        timedwait(() -> all(istaskdone, handle.readers), runner.cleanup_seconds; pollint=0.01)
        all(istaskdone, handle.readers) || return false
        lock(runner.lock) do
            get(runner.active, handle.id, nothing) === handle && delete!(runner.active, handle.id)
        end
        return true
    end
end

"""
    run_owned_command!(runner, command; timeout_seconds=runner.timeout_seconds,
                       token=nothing) -> CommandResult

Execute one operator-built Cmd without a shell or inherited stdin. Bound both
output streams together, process lifetime and reader retirement. A nonzero exit
returns a result; cancellation, deadline, overflow and cleanup failure throw a
fixed CommandFailure outside the original exception context. An optional execution
core CancellationToken interrupts only this command.
"""
function run_owned_command!(runner::CommandRunner, command::Cmd;
        timeout_seconds=runner.timeout_seconds, token=nothing)
    timeout_seconds isa Real && !(timeout_seconds isa Bool) && isfinite(timeout_seconds) &&
        0 < timeout_seconds <= 120 || throw(ArgumentError("host command timeout must be in (0, 120] seconds"))
    token === nothing || token isa ExecutionCore.CancellationToken ||
        throw(ArgumentError("host command cancellation requires an owned token"))
    handle = OwnedCommand(uuid4(), nothing, (Pipe(), Pipe()), Task[], 0, nothing, ReentrantLock())
    failure = nothing
    result = nothing
    admitted = false
    try
        lock(runner.lock) do
            runner.closed && throw(CommandFailure(:closed))
            length(runner.active) < runner.capacity || throw(CommandFailure(:busy))
            token === nothing || !ExecutionCore.iscanceled(token) || throw(CommandFailure(:canceled))
            runner.active[handle.id] = handle
            admitted = true
            handle.process = run(pipeline(ignorestatus(command);
                stdin=devnull, stdout=handle.pipes[1], stderr=handle.pipes[2]); wait=false)
            foreach(pipe -> close(pipe.in), handle.pipes)
            for pipe in handle.pipes
                push!(handle.readers, @async command_reader!(runner, handle, pipe))
            end
        end
        started = time_ns()
        while true
            runner.closed && throw(CommandFailure(:closed))
            token === nothing || !ExecutionCore.iscanceled(token) || throw(CommandFailure(:canceled))
            handle.failure === nothing || throw(CommandFailure(handle.failure))
            if !Base.process_running(handle.process) && all(istaskdone, handle.readers)
                result = CommandResult(Int(handle.process.exitcode),
                    fetch(handle.readers[1]), fetch(handle.readers[2]))
                break
            end
            (time_ns() - started) / 1e9 < timeout_seconds || throw(CommandFailure(:deadline))
            sleep(0.01)
        end
    catch error
        failure = error isa CommandFailure ? error.code : :command_failed
    finally
        if admitted
            retired = try
                retire_command!(runner, handle)
            catch
                false
            end
            retired || (failure = :cleanup_unresolved)
        else
            foreach(close, handle.pipes)
        end
    end
    failure === nothing || throw(CommandFailure(failure))
    return result::CommandResult
end

"""
    close(runner::CommandRunner)

Close admission, join all occupied command retirements and retain unresolved
handles. Repeating close retries those same resources without reopening admission.
"""
function Base.close(runner::CommandRunner)
    lock(runner.cleanup_lock) do
        handles = lock(runner.lock) do
            runner.closed = true
            collect(values(runner.active))
        end
        tasks = [@async(try retire_command!(runner, handle) catch; false end) for handle in handles]
        results = [try fetch(task) catch; false end for task in tasks]
        all(results) || throw(CommandFailure(:cleanup_unresolved))
    end
    return nothing
end

"""
    container_command_environment()

Return the host CLI's minimal environment. Home/runtime-directory references are
for the operator's local container configuration, never mounted or forwarded to
user code. Broker, storage, proxy and remote-engine override variables are absent.
"""
function container_command_environment()
    result = Dict(name => ENV[name] for name in
        ("PATH", "HOME", "USER", "LOGNAME", "XDG_RUNTIME_DIR", "DBUS_SESSION_BUS_ADDRESS")
        if haskey(ENV, name))
    result["LC_ALL"] = "C"
    return result
end

function container_probe(runner::CommandRunner, arguments::Vector{String})
    result = run_owned_command!(runner, setenv(Cmd(arguments), container_command_environment()))
    # Runtime warnings are private; discovery only needs stdout plus the shim marker.
    shim = length(arguments) == 2 && arguments[2] == "version" &&
        occursin("podman", lowercase(result.diagnostic)) ? "\nPodman" : ""
    return result.exitcode == 0, result.output * shim
end
