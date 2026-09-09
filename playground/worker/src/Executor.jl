mutable struct CancellationToken
    requested::Threads.Atomic{Bool}
end

CancellationToken() = CancellationToken(Threads.Atomic{Bool}(false))
cancel!(token::CancellationToken) = (token.requested[] = true)
iscanceled(token::CancellationToken) = token.requested[]

mutable struct CancellationRegistry
    pending::Dict{String,Float64}
    active::Dict{String,CancellationToken}
    ttl_seconds::Float64
    lock::ReentrantLock
end

CancellationRegistry(; ttl_seconds::Real=86_400) = CancellationRegistry(
    Dict{String,Float64}(),
    Dict{String,CancellationToken}(),
    Float64(ttl_seconds),
    ReentrantLock()
)

function prune_cancellations!(registry::CancellationRegistry)
    cutoff = time() - registry.ttl_seconds
    filter!(pair -> last(pair) >= cutoff, registry.pending)
    return registry
end

function request_cancellation!(registry::CancellationRegistry, job_id::AbstractString)
    lock(registry.lock) do
        prune_cancellations!(registry)
        id = string(job_id)
        token = get(registry.active, id, nothing)
        if isnothing(token)
            registry.pending[id] = time()
        else
            cancel!(token)
        end
    end
    return registry
end

function activate_cancellation!(registry::CancellationRegistry, job_id::AbstractString)
    return lock(registry.lock) do
        prune_cancellations!(registry)
        id = string(job_id)
        token = CancellationToken()
        if pop!(registry.pending, id, nothing) !== nothing
            cancel!(token)
        end
        registry.active[id] = token
        token
    end
end

function finish_cancellation!(registry::CancellationRegistry, job_id::AbstractString)
    lock(registry.lock) do
        pop!(registry.active, string(job_id), nothing)
        pop!(registry.pending, string(job_id), nothing)
    end
    return registry
end

struct OperationCanceled <: Exception end
Base.showerror(io::IO, ::OperationCanceled) = print(io, "operation canceled")

struct OperationDeadlineExpired <: Exception end
Base.showerror(io::IO, ::OperationDeadlineExpired) = print(
    io,
    "operation deadline expired"
)

mutable struct ExecutionContext
    job_id::String
    token::CancellationToken
    emit_progress::Function
    emit_log::Function
    warnings::Vector{String}
    prepared_cache::PreparedResourceCache
    deadline::Union{Nothing,DateTime}
end

function check_canceled(context::ExecutionContext)
    iscanceled(context.token) && throw(OperationCanceled())
    !isnothing(context.deadline) && Dates.now(Dates.UTC) >= context.deadline &&
        throw(OperationDeadlineExpired())
    return nothing
end

function progress!(
        context::ExecutionContext,
        fraction::Real,
        stage::AbstractString;
        message::Union{Nothing,AbstractString}=nothing
    )
    check_canceled(context)
    0 <= fraction <= 1 || throw(ArgumentError("progress must be between zero and one"))
    context.emit_progress(
        Float64(fraction),
        string(stage),
        isnothing(message) ? nothing : string(message)
    )
    return nothing
end

log!(context::ExecutionContext, message::AbstractString) =
    context.emit_log(string(message))

function bounded_warning_text(message; maximum_bytes::Integer=4 * 1024)
    text = string(message)
    ncodeunits(text) <= maximum_bytes && return text
    io = IOBuffer()
    used = 0
    for character in text
        bytes = ncodeunits(string(character))
        used + bytes > maximum_bytes - 3 && break
        write(io, character)
        used += bytes
    end
    return String(take!(io)) * "..."
end

function warning!(context::ExecutionContext, message)
    text = bounded_warning_text(message)
    if length(context.warnings) < 128 && !(text in context.warnings)
        push!(context.warnings, text)
    end
    context.emit_log("Warning: $text")
    return nothing
end

function execute_operation(
        spec::OperationSpec,
        context::ExecutionContext,
        parameters::Dict{String,Any}
    )
    validated = Base.invokelatest(spec.validator, parameters)
    check_canceled(context)
    result = Base.invokelatest(spec.executor, context, validated)
    check_canceled(context)
    return normalize_wire(result)
end

mutable struct ExecutorSupervisor
    project::String
    process::Union{Nothing,Base.Process}
    lock::ReentrantLock
    generation::Int
    startup_timeout_seconds::Float64
    "Optional operator-owned profile command; never browser-supplied."
    command::Union{Nothing,Cmd}
    "Discard private adapter stderr without adding another process or pipe task."
    discard_stderr::Bool
    "Serialize process retirement without waiting for the operation lock."
    stop_lock::ReentrantLock
    "Owned bounded pipe tasks, joined before forgetting a process."
    io_tasks::Set{Task}
end

const EXECUTOR_FRAME_PREFIX = "@LCM_EXECUTOR_FRAME@"
const EXECUTOR_MAX_LINE_BYTES = 1024 * 1024

"""
    ExecutorSupervisor(project=dirname(Base.active_project()); startup_timeout_seconds=30,
                       command=nothing, discard_stderr=false)

Supervise one scientific process with a separately bounded bootstrap. A bootstrap
acknowledgement means that the command reader is available, not that scientific
packages or a model have been prepared. Operation deadlines remain independent.
An optional `command` is an operator-owned local launch command for a profile;
it must never be constructed from a browser-supplied executable or arguments.
Set `discard_stderr=true` for a container adapter whose engine diagnostics must
not enter ordinary logs. This redirects only stderr to the null device; arbitrary
pipelines remain unsupported because this owner supervises exactly one process.
"""
function ExecutorSupervisor(project::AbstractString=dirname(Base.active_project()); startup_timeout_seconds=30,
        command::Union{Nothing,Cmd}=nothing,discard_stderr::Bool=false)
    startup_timeout_seconds isa Real && !(startup_timeout_seconds isa Bool) &&
        isfinite(startup_timeout_seconds) && 0 < startup_timeout_seconds <= 120 ||
        throw(ArgumentError("executor startup timeout must be in (0, 120] seconds"))
    return ExecutorSupervisor(abspath(project), nothing, ReentrantLock(), 0,
        Float64(startup_timeout_seconds), command, discard_stderr, ReentrantLock(), Set{Task}())
end

function executor_command(supervisor::ExecutorSupervisor)
    supervisor.command === nothing || return supervisor.command
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    command = `$julia --startup-file=no --compiled-modules=existing --project=$(supervisor.project) -e "using LineCableModelsWorker; LineCableModelsWorker.executor_main()"`
    return addenv(command, "LCM_EXECUTOR_PROCESS" => "1")
end

function start_executor!(supervisor::ExecutorSupervisor, context::Union{Nothing,ExecutionContext}=nothing)
    return lock(() -> start_executor_locked!(supervisor, context), supervisor.lock)
end

function start_executor_locked!(supervisor::ExecutorSupervisor, context::Union{Nothing,ExecutionContext})
    context === nothing || check_canceled(context)
    process = supervisor.process
    if !isnothing(process) && process_running(process)
        return process
    end
    # Reap/close a naturally exited predecessor before replacing its handle.
    process === nothing || stop_executor!(supervisor)
    process = lock(supervisor.stop_lock) do
        command = executor_command(supervisor)
        owned = supervisor.discard_stderr ? pipeline(command;stderr=devnull) : command
        supervisor.process = open(owned, "r+")
        supervisor.generation += 1
        supervisor.process
    end
    started = time_ns()
    try
        while true
            response_task = executor_response_task(supervisor, process)
            while !istaskdone(response_task)
                context === nothing || check_canceled(context)
                (time_ns() - started) / 1e9 <= supervisor.startup_timeout_seconds ||
                    throw(RetryableOperationError("Scientific executor startup timed out"))
                sleep(0.025)
            end
            response = fetch_executor_io!(supervisor, response_task)
            response.type == "bootstrapped" && break
            response.type == "engine_output" ||
                throw(RetryableOperationError("Scientific executor bootstrap protocol failed"))
            context === nothing || context.emit_log(bounded_warning_text(response.message))
            (time_ns() - started) / 1e9 <= supervisor.startup_timeout_seconds ||
                throw(RetryableOperationError("Scientific executor startup timed out"))
        end
        context === nothing || check_canceled(context)
    catch error
        stop_executor!(supervisor)
        error isa Union{OperationCanceled,OperationDeadlineExpired,RetryableOperationError} && rethrow()
        throw(RetryableOperationError("Scientific executor stopped during startup"))
    end
    return process
end

"""
    stop_executor!(supervisor; grace_seconds=2, kill_seconds=2)

Retire the owned process and join its pipe tasks. Allow a finite shutdown-command
grace, then terminate, then kill if necessary. Each signal wait and pipe-task join
is bounded. Return nothing only after physical process exit and pipe cleanup.
An unresolved stop throws RetryableOperationError and retains the process handle
and pending tasks, allowing its owner to retry without releasing capacity.

This stops one trusted scientific process, not an arbitrary subprocess tree.
Container/process-group ownership belongs to the agent's isolation adapter.
"""
function stop_executor!(supervisor::ExecutorSupervisor; grace_seconds=2, kill_seconds=2)
    all(v -> v isa Real && !(v isa Bool) && isfinite(v) && 0 < v <= 10,
        (grace_seconds, kill_seconds)) || throw(ArgumentError("executor stop bounds must be in (0, 10] seconds"))
    return lock(supervisor.stop_lock) do
        process = supervisor.process
        isnothing(process) && return nothing
        if process_running(process)
            # A blocked stdin writer must never block retirement itself.
            executor_io_task(supervisor, process) do
                write(process, JSON3.write(Dict("command" => "shutdown")), '\n')
                flush(process)
            end
            timedwait(() -> !process_running(process), grace_seconds; pollint=0.025)
        end
        for (signal, seconds) in ((Base.SIGTERM, grace_seconds), (Base.SIGKILL, kill_seconds))
            process_running(process) || break
            try
                kill(process, signal)
            catch
                # The exit check below, not a successful signal call, decides cleanup.
            end
            timedwait(() -> !process_running(process), seconds; pollint=0.025)
        end
        process_running(process) && throw(RetryableOperationError("Scientific executor cleanup remains unresolved"))
        # Process exit is observed before wait; do not join arbitrary command-
        # supplied sync tasks. All scientific pipe tasks are explicitly owned here.
        wait(process, false)
        close(process)
        timedwait(() -> all(istaskdone, supervisor.io_tasks), kill_seconds; pollint=0.025)
        all(istaskdone, supervisor.io_tasks) ||
            throw(RetryableOperationError("Scientific executor pipe cleanup remains unresolved"))
        empty!(supervisor.io_tasks)
        supervisor.process = nothing
        return nothing
    end
end

function decode_executor_line(line::AbstractString)
    if startswith(line, EXECUTOR_FRAME_PREFIX)
        payload = line[(ncodeunits(EXECUTOR_FRAME_PREFIX) + 1):end]
        return JSON3.read(payload)
    end
    return (
        type="engine_output",
        message=string(line),
    )
end

function read_executor_line(io; maximum_bytes::Int=EXECUTOR_MAX_LINE_BYTES)
    maximum_bytes > 0 || throw(ArgumentError("executor line limit must be positive"))
    bytes = UInt8[]
    sizehint!(bytes, min(maximum_bytes, 4096))
    while !eof(io)
        byte = read(io, UInt8)
        byte == 0x0a && return String(bytes)
        length(bytes) < maximum_bytes || throw(ArgumentError("executor line exceeds its byte limit"))
        push!(bytes, byte)
    end
    isempty(bytes) && throw(EOFError())
    return String(bytes)
end

function executor_io_task(action::Function, supervisor::ExecutorSupervisor, process::Base.Process)
    return lock(supervisor.stop_lock) do
        supervisor.process === process || throw(RetryableOperationError("Scientific executor was retired"))
        filter!(task -> !istaskdone(task), supervisor.io_tasks)
        length(supervisor.io_tasks) < 4 || throw(RetryableOperationError("Scientific executor pipe task limit exceeded"))
        task = Task(action)
        push!(supervisor.io_tasks, task)
        schedule(task)
        task
    end
end

function fetch_executor_io!(supervisor::ExecutorSupervisor, task::Task)
    try
        return fetch(task)
    finally
        lock(supervisor.stop_lock) do
            delete!(supervisor.io_tasks, task)
        end
    end
end

function executor_response_task(supervisor::ExecutorSupervisor, process::Base.Process)
    return executor_io_task(supervisor, process) do
        line = read_executor_line(process)
        return decode_executor_line(line)
    end
end

function check_executor_deadline(context::ExecutionContext, started::UInt64, seconds::Real)
    check_canceled(context)
    (time_ns() - started) / 1e9 <= seconds || throw(PermanentOperationError(
        "timeout", "Calculation exceeded its $seconds second limit"))
    return nothing
end

"""
    execute_supervised!(supervisor, spec, context, parameters;
        control=:operation, preparation_key=nothing, on_result_schema=nothing)

Execute one bounded request in the owned child and return its decoded result.
With `preparation_key`, require the existing prepared process; do not start a
replacement implicitly. Cancellation, deadline and transport failures stop the
child. When supplied, `on_result_schema` receives the validated schema version
from the child's operation registry before the result is returned. The default
preserves the result-only interface used by legacy callers.
"""
function execute_supervised!(
        supervisor::ExecutorSupervisor,
        spec::OperationSpec,
        context::ExecutionContext,
        parameters::Dict{String,Any};
        control::Symbol=:operation,
        preparation_key::Union{Nothing,String}=nothing,
        on_result_schema=nothing
    )
    return lock(supervisor.lock) do
        control in (:operation, :prepare, :preparation_status) || throw(ArgumentError("unsupported executor control"))
        preparation_key === nothing || occursin(r"^[a-f0-9]{64}$", preparation_key) ||
            throw(ArgumentError("invalid preparation identity"))
        command = JSON3.write(Dict(
            (control == :operation ? "operation" : "command") => (control == :operation ? spec.name : String(control)),
            "parameters" => parameters,
            "preparation_key" => preparation_key,
        ))
        ncodeunits(command) <= EXECUTOR_MAX_LINE_BYTES || throw(ArgumentError("executor command exceeds its byte limit"))
        process = if control == :preparation_status || preparation_key !== nothing
            lock(supervisor.stop_lock) do
                current = supervisor.process
                current !== nothing && process_running(current) ||
                    throw(RetryableOperationError("Prepared scientific executor is not running"))
                check_canceled(context)
                current
            end
        else
            start_executor!(supervisor, context)
        end
        started = time_ns()
        try
            writer = executor_io_task(supervisor, process) do
                write(process, command, '\n')
                flush(process)
            end
            while !istaskdone(writer)
                check_executor_deadline(context, started, spec.timeout_seconds)
                sleep(0.025)
            end
            fetch_executor_io!(supervisor, writer)
            check_executor_deadline(context, started, spec.timeout_seconds)
        catch error
            stop_executor!(supervisor)
            error isa Union{OperationCanceled,OperationDeadlineExpired,PermanentOperationError} && rethrow()
            throw(RetryableOperationError(
                "Could not send work to the scientific executor"
            ))
        end

        while true
            response_task = executor_response_task(supervisor, process)
            while !istaskdone(response_task)
                if iscanceled(context.token)
                    stop_executor!(supervisor)
                    throw(OperationCanceled())
                elseif !isnothing(context.deadline) &&
                        Dates.now(Dates.UTC) >= context.deadline
                    stop_executor!(supervisor)
                    throw(OperationDeadlineExpired())
                elseif (time_ns() - started) / 1e9 > spec.timeout_seconds
                    stop_executor!(supervisor)
                    throw(PermanentOperationError(
                        "timeout",
                        "Calculation exceeded its $(spec.timeout_seconds) second limit"
                    ))
                end
                sleep(0.025)
            end
            response = try
                fetch_executor_io!(supervisor, response_task)
            catch error
                stop_executor!(supervisor)
                error isa RetryableOperationError && rethrow()
                throw(RetryableOperationError(
                    "scientific executor stopped: $(sprint(showerror, error))"
                ))
            end
            try
                check_canceled(context)
                (time_ns() - started) / 1e9 <= spec.timeout_seconds || throw(PermanentOperationError(
                    "timeout", "Calculation exceeded its $(spec.timeout_seconds) second limit"))
            catch
                stop_executor!(supervisor)
                rethrow()
            end
            kind = try
                string(response.type)
            catch
                stop_executor!(supervisor)
                throw(RetryableOperationError("Scientific executor response is missing its type"))
            end
            try
                if kind == "progress"
                    context.emit_progress(
                        Float64(response.progress),
                        string(response.stage),
                        haskey(response, :message) && !isnothing(response.message) ?
                            string(response.message) : nothing
                    )
                elseif kind == "log"
                    context.emit_log(string(response.message))
                elseif kind == "engine_output"
                    context.emit_log(string(response.message))
                elseif kind == "warning"
                    warning!(context, response.message)
                elseif kind == "result"
                    if on_result_schema !== nothing
                        schema = get(response,:schema_version,nothing)
                        schema isa String && ncodeunits(schema)<=32 && occursin(r"^[0-9]+(?:\.[0-9]+){1,2}$",schema) ||
                            throw(RetryableOperationError("Scientific result schema is missing or invalid"))
                        on_result_schema(schema)
                    end
                    return normalize_wire(response.result)
                elseif kind == "permanent_error"
                    throw(PermanentOperationError(
                        string(response.category),
                        string(response.message)
                    ))
                elseif kind == "retryable_error"
                    stop_executor!(supervisor)
                    throw(RetryableOperationError(string(response.message)))
                else
                    stop_executor!(supervisor)
                    throw(RetryableOperationError(
                        "scientific executor returned an unknown response"
                    ))
                end
            catch error
                # A failed observer or malformed response must not leave old
                # output queued to be consumed as the next operation's result.
                if !(kind == "permanent_error" && error isa PermanentOperationError)
                    stop_executor!(supervisor)
                end
                rethrow()
            end
        end
    end
end

"""
    prepare_supervised!(supervisor, context, parameters; timeout_seconds=300)

Explicitly prepare the current profile process through its existing framed
channel. Return representative-workload evidence only after that child completes
its hook. Cancellation, timeout or process loss terminates the child; a subsequent
start has a new supervisor generation and must be prepared again.
"""
function prepare_supervised!(supervisor::ExecutorSupervisor, context::ExecutionContext,
        parameters::Dict{String,Any}; timeout_seconds::Real=300)
    spec = OperationSpec("runtime.prepare", identity, (_, _) -> nothing;
        timeout_seconds, execution_mode=:supervised, cache_policy=:none)
    return execute_supervised!(supervisor, spec, context, parameters; control=:prepare)
end

"""
    inspect_preparation!(supervisor, context, key; timeout_seconds=5)

Query preparation evidence in the current child without building a model or
extending its cache lifetime. The caller must already own a live process; this
function never starts a cold replacement merely to ask whether it is ready.
"""
function inspect_preparation!(supervisor::ExecutorSupervisor, context::ExecutionContext,
        key::String; timeout_seconds::Real=5)
    occursin(r"^[a-f0-9]{64}$", key) || throw(ArgumentError("invalid preparation identity"))
    return lock(supervisor.lock) do
        process = supervisor.process
        process !== nothing && process_running(process) ||
            throw(RetryableOperationError("Prepared scientific executor is not running"))
        spec = OperationSpec("runtime.preparation_status", identity, (_,_) -> nothing;
            timeout_seconds, execution_mode=:supervised, cache_policy=:none)
        execute_supervised!(supervisor, spec, context, Dict{String,Any}("key"=>key); control=:preparation_status)
    end
end

function executor_write(message)
    encoded = JSON3.write(message)
    ncodeunits(encoded) + ncodeunits(EXECUTOR_FRAME_PREFIX) <= EXECUTOR_MAX_LINE_BYTES ||
        throw(ArgumentError("executor response exceeds its byte limit"))
    println(stdout, EXECUTOR_FRAME_PREFIX, encoded)
    flush(stdout)
    return nothing
end

struct ExecutorWireLogger <: Logging.AbstractLogger end

Logging.min_enabled_level(::ExecutorWireLogger) = Logging.Warn
Logging.shouldlog(::ExecutorWireLogger, level, _module, group, id) =
    level >= Logging.Warn
Logging.catch_exceptions(::ExecutorWireLogger) = true

function Logging.handle_message(
        ::ExecutorWireLogger,
        level,
        message,
        _module,
        group,
        id,
        file,
        line;
        kwargs...
    )
    executor_write(Dict(
        "type" => (level >= Logging.Error ? "engine_output" : "warning"),
        "message" => bounded_warning_text(message),
    ))
    return nothing
end

function executor_main()
    return executor_main(default_registry(), operation ->
        startswith(operation, "system.") ? nothing : load_scientific_packages!())
end

function executor_main(registry::OperationRegistry, load_operation!;
        prepared_cache=PreparedResourceCache(), prepare=nothing, preparation_status=nothing, cleanup=() -> nothing)
    executor_write(Dict("type" => "bootstrapped"))
    try
        while !eof(stdin)
            line = read_executor_line(stdin)
            isempty(strip(line)) && continue
            try
                command = JSON3.read(line)
                if haskey(command, :command) && command.command == "shutdown"
                    break
                end
                preparing = haskey(command, :command) && command.command == "prepare"
                inspecting = haskey(command, :command) && command.command == "preparation_status"
                if haskey(command, :command) && !preparing && !inspecting
                    throw(PermanentOperationError("invalid_command", "Unknown executor command"))
                end
                preparing && prepare === nothing && throw(PermanentOperationError(
                    "preparation_unavailable", "This executor has no preparation contract"))
                inspecting && preparation_status === nothing && throw(PermanentOperationError(
                    "preparation_unavailable", "This executor has no preparation status contract"))
                spec = if preparing || inspecting
                    nothing
                else
                    operation = string(command.operation)
                    registered = registered_operation(registry, operation)
                    registered.execution_mode == :supervised || throw(PermanentOperationError(
                        "execution_boundary", "Operation $operation is not allowed in the scientific executor"))
                    load_operation!(operation)
                    registered
                end
                preparation_key = get(command, :preparation_key, nothing)
                if !preparing && !inspecting && preparation_key !== nothing
                    preparation_status === nothing && throw(PermanentOperationError(
                        "preparation_unavailable", "This executor cannot verify preparation"))
                    status = Base.invokelatest(preparation_status, String(preparation_key))
                    status["ready"] === true || throw(PermanentOperationError(
                        "preparation_expired", "Prepare this executor before executing the operation"))
                end
                token = CancellationToken()
                context = ExecutionContext(
                    "executor",
                    token,
                    (progress, stage, message) -> executor_write(Dict(
                        "type" => "progress",
                        "progress" => progress,
                        "stage" => stage,
                        "message" => message,
                    )),
                    message -> executor_write(Dict(
                        "type" => "log",
                        "message" => message,
                    )),
                    String[],
                    prepared_cache,
                    nothing
                )
                result = Logging.with_logger(ExecutorWireLogger()) do
                    if preparing
                        Base.invokelatest(prepare, context, normalize_wire(command.parameters))
                    elseif inspecting
                        Base.invokelatest(preparation_status, String(command.parameters.key))
                    else
                        execute_operation(spec, context, normalize_wire(command.parameters))
                    end
                end
                executor_write(Dict("type" => "result", "result" => result,
                    "schema_version" => spec === nothing ? "1.0" : spec.schema_version))
            catch error
                if error isa PermanentOperationError
                    executor_write(Dict(
                        "type" => "permanent_error",
                        "category" => error.category,
                        "message" => error.message,
                    ))
                else
                    diagnostic_id = string(uuid4())
                    @error "Scientific executor operation failed" diagnostic_id exception=(
                        error,
                        catch_backtrace()
                    )
                    executor_write(Dict(
                        "type" => "retryable_error",
                        "message" => "Scientific executor failed; diagnostic $diagnostic_id",
                    ))
                end
            end
        end
    finally
        cleanup()
    end
    return nothing
end
