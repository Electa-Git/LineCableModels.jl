module ShowcaseREPLWorker

using Serialization

const PROTOCOL_OUTPUT = stdout
const PROTOCOL_LOCK = ReentrantLock()
const MAX_OUTPUT_BYTES = 1_048_576

function emit(event)
    lock(PROTOCOL_LOCK) do
        serialize(PROTOCOL_OUTPUT, event)
        flush(PROTOCOL_OUTPUT)
    end
    return nothing
end

function forward_stream(
        pipe::Pipe,
        kind::Symbol,
        emitted::Base.RefValue{Int},
        truncated::Base.RefValue{Bool}
)
    try
        while true
            bytes = readavailable(pipe)
            if !isempty(bytes)
                remaining = max(0, MAX_OUTPUT_BYTES - emitted[])
                count = min(length(bytes), remaining)
                if count > 0
                    emit((; kind, data = bytes[1:count]))
                    emitted[] += count
                end
                if count < length(bytes) && !truncated[]
                    truncated[] = true
                    emit((;
                        kind = :system,
                        data = "\n[output truncated after 1 MiB]\n"
                    ))
                end
            end
            eof(pipe) && break
            yield()
        end
    catch error
        error isa EOFError || emit((;
            kind = :system,
            data = "\n[output forwarding failed: $(sprint(showerror, error))]\n"
        ))
    end
    return nothing
end

function displayed_result(value)
    io = IOBuffer()
    context = IOContext(io, :limit => true, :compact => false, :displaysize => (24, 100))
    show(context, MIME("text/plain"), value)
    return String(take!(io))
end

function displayed_error(error, backtrace)
    io = IOBuffer()
    context = IOContext(io, :limit => true, :displaysize => (24, 100))
    showerror(context, error, backtrace)
    return String(take!(io))
end

function evaluate!(workspace::Module, source::String)
    output_pipe = Pipe()
    error_pipe = Pipe()
    emitted = Ref(0)
    truncated = Ref(false)
    output_task = @async forward_stream(output_pipe, :stdout, emitted, truncated)
    error_task = @async forward_stream(error_pipe, :stderr, emitted, truncated)

    result = nothing
    failure = nothing
    backtrace = nothing
    try
        result = redirect_stdout(output_pipe) do
            redirect_stderr(error_pipe) do
                expression = Meta.parseall(source; filename = "showcase-repl")
                Core.eval(workspace, expression)
            end
        end
        Core.eval(workspace, Expr(:(=), :ans, QuoteNode(result)))
    catch error
        failure = error
        backtrace = catch_backtrace()
    finally
        close(output_pipe.in)
        close(error_pipe.in)
        wait(output_task)
        wait(error_task)
        close(output_pipe.out)
        close(error_pipe.out)
    end

    if isnothing(failure)
        suppressed = endswith(rstrip(source), ";")
        if !suppressed && !isnothing(result)
            try
                emit((; kind = :result, data = displayed_result(result)))
            catch error
                emit((;
                    kind = :error,
                    data = "Unable to display result: $(sprint(showerror, error))"
                ))
            end
        end
    elseif failure isa InterruptException
        emit((; kind = :interrupted, data = "Interrupted"))
    else
        emit((; kind = :error, data = displayed_error(failure, backtrace)))
    end
    emit((; kind = :done))
    return nothing
end

function workspace_module()
    workspace = Module(gensym(:ShowcaseWorkspace))
    Core.eval(workspace, :(using Base))
    return workspace
end

function serve()
    Base.exit_on_sigint(false)
    workspace = workspace_module()
    emit((; kind = :ready, pid = getpid(), version = string(VERSION)))
    while true
        request = try
            deserialize(stdin)
        catch error
            error isa EOFError && return
            rethrow()
        end
        request.kind === :evaluate && evaluate!(workspace, String(request.source))
        request.kind === :shutdown && return
    end
end

serve()

end
