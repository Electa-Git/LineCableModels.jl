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
end

const EXECUTOR_FRAME_PREFIX = "@LCM_EXECUTOR_FRAME@"

function ExecutorSupervisor(project::AbstractString=dirname(Base.active_project()))
    return ExecutorSupervisor(abspath(project), nothing, ReentrantLock(), 0)
end

function executor_command(supervisor::ExecutorSupervisor)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    command = `$julia --startup-file=no --project=$(supervisor.project) -e "using LineCableModelsWorker; LineCableModelsWorker.executor_main()"`
    return addenv(command, "LCM_EXECUTOR_PROCESS" => "1")
end

function start_executor!(supervisor::ExecutorSupervisor)
    process = supervisor.process
    if !isnothing(process) && process_running(process)
        return process
    end
    supervisor.process = open(executor_command(supervisor), "r+")
    supervisor.generation += 1
    return supervisor.process
end

function stop_executor!(supervisor::ExecutorSupervisor)
    process = supervisor.process
    isnothing(process) && return nothing
    if process_running(process)
        try
            write(process, JSON3.write(Dict("command" => "shutdown")), '\n')
            flush(process)
        catch
        end
        timedwait(() -> !process_running(process), 2.0; pollint=0.025)
    end
    if process_running(process)
        try
            kill(process, Base.SIGKILL)
        catch
        end
    end
    try
        wait(process)
    catch
    end
    supervisor.process = nothing
    return nothing
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

function executor_response_task(process)
    return @async begin
        line = readline(process)
        return decode_executor_line(line)
    end
end

function execute_supervised!(
        supervisor::ExecutorSupervisor,
        spec::OperationSpec,
        context::ExecutionContext,
        parameters::Dict{String,Any}
    )
    return lock(supervisor.lock) do
        process = start_executor!(supervisor)
        command = JSON3.write(Dict(
            "operation" => spec.name,
            "parameters" => parameters,
        ))
        try
            write(process, command, '\n')
            flush(process)
        catch error
            stop_executor!(supervisor)
            throw(RetryableOperationError(
                "could not send work to the scientific executor: $(sprint(showerror, error))"
            ))
        end

        started = time()
        while true
            response_task = executor_response_task(process)
            while !istaskdone(response_task)
                if iscanceled(context.token)
                    stop_executor!(supervisor)
                    throw(OperationCanceled())
                elseif !isnothing(context.deadline) &&
                        Dates.now(Dates.UTC) >= context.deadline
                    stop_executor!(supervisor)
                    throw(OperationDeadlineExpired())
                elseif time() - started > spec.timeout_seconds
                    stop_executor!(supervisor)
                    throw(PermanentOperationError(
                        "timeout",
                        "Calculation exceeded its $(spec.timeout_seconds) second limit"
                    ))
                end
                sleep(0.025)
            end
            response = try
                fetch(response_task)
            catch error
                stop_executor!(supervisor)
                error isa RetryableOperationError && rethrow()
                throw(RetryableOperationError(
                    "scientific executor stopped: $(sprint(showerror, error))"
                ))
            end
            kind = string(response.type)
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
        end
    end
end

function executor_write(message)
    println(stdout, EXECUTOR_FRAME_PREFIX, JSON3.write(message))
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
    registry = default_registry()
    prepared_cache = PreparedResourceCache()
    for line in eachline(stdin)
        isempty(strip(line)) && continue
        try
            command = JSON3.read(line)
            if haskey(command, :command) && command.command == "shutdown"
                break
            end
            operation = string(command.operation)
            spec = registered_operation(registry, operation)
            spec.execution_mode == :supervised || throw(PermanentOperationError(
                "execution_boundary",
                "Operation $operation is not allowed in the scientific executor"
            ))
            startswith(operation, "system.") || load_scientific_packages!()
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
                execute_operation(
                    spec,
                    context,
                    normalize_wire(command.parameters)
                )
            end
            executor_write(Dict("type" => "result", "result" => result))
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
    return nothing
end
