module LiveJuliaPage

using Bonito
using Markdown
using Observables: off
using Serialization
using ..PageAuthoring

const SHOWCASE_DIRECTORY = dirname(@__DIR__)
const REPOSITORY_DIRECTORY = normpath(joinpath(SHOWCASE_DIRECTORY, "..", ".."))
const WORKER_SCRIPT = joinpath(SHOWCASE_DIRECTORY, "repl_worker.jl")
const INITIAL_SOURCE = """using LineCableModels

# This namespace remains alive between evaluations.
exported_names = names(LineCableModels; all = false)
first(exported_names, min(12, length(exported_names)))
"""

mutable struct JuliaWorker
    process::Base.Process
    pid::Int
    version::String
end

mutable struct WorkerController
    worker::Union{Nothing, JuliaWorker}
    epoch::Int
    busy::Bool
    closed::Bool
    evaluation::Int
end

function worker_command()
    project = dirname(Base.active_project())
    arguments = [Base.julia_cmd().exec;
                 "--project=$project";
                 "--startup-file=no";
                 "--history-file=no";
                 "--color=yes";
                 WORKER_SCRIPT]
    return Cmd(Cmd(arguments); dir = REPOSITORY_DIRECTORY, ignorestatus = true)
end

function start_worker()
    process = open(worker_command(), "r+")
    event = try
        deserialize(process)
    catch
        close(process)
        rethrow()
    end
    event.kind === :ready || begin
        close(process)
        error("Julia worker did not send its ready event")
    end
    return JuliaWorker(process, event.pid, event.version)
end

worker_running(worker::JuliaWorker) = Base.process_running(worker.process)

function stop_worker!(worker::JuliaWorker; force::Bool = false)
    if worker_running(worker)
        try
            if force
                Base.kill(worker.process, Base.SIGKILL)
            else
                serialize(worker.process, (; kind = :shutdown))
                flush(worker.process)
            end
        catch error
            if error isa Base.IOError
                worker_running(worker) && Base.kill(worker.process, Base.SIGKILL)
            else
                rethrow()
            end
        end
    end
    try
        wait(worker.process)
    catch error
        error isa Base.ProcessFailedException || rethrow()
    end
    close(worker.process)
    return nothing
end

function evaluate!(callback, worker::JuliaWorker, source::AbstractString)
    worker_running(worker) || error("Julia worker is not running")
    serialize(worker.process, (; kind = :evaluate, source = String(source)))
    flush(worker.process)
    while true
        event = deserialize(worker.process)
        callback(event)
        event.kind === :done && return nothing
    end
end

function prompt(source::AbstractString, evaluation::Integer)
    lines = split(replace(String(source), "\r\n" => "\n"), '\n'; keepempty = true)
    isempty(lines) && return "\e[1;36mjulia>\e[0m\n"
    rendered = String["\e[1;36mjulia>\e[0m $(first(lines))"]
    append!(rendered, ("       $line" for line in Iterators.drop(lines, 1)))
    return "\e[2m# evaluation $evaluation\e[0m\n$(join(rendered, '\n'))\n"
end

function append_event!(terminal::TerminalOutput, event)
    if event.kind === :stdout || event.kind === :stderr
        push!(terminal, event.data)
    elseif event.kind === :result
        push!(terminal, "\n\e[32m$(event.data)\e[0m\n")
    elseif event.kind === :error
        push!(terminal, "\n\e[1;31mERROR:\e[0m $(event.data)\n")
    elseif event.kind === :interrupted
        push!(terminal, "\n\e[1;33m$(event.data)\e[0m\n")
    elseif event.kind === :system
        push!(terminal, "\e[33m$(event.data)\e[0m")
    end
    return nothing
end

function worker_banner(worker::JuliaWorker)
    return "\e[1;36mJulia $(worker.version)\e[0m · isolated worker PID $(worker.pid)\n" *
           "State persists between runs. Ctrl+Enter evaluates the editor.\n\n"
end

function build(session::Session)
    editor = CodeEditor(
        "julia";
        initial_source = INITIAL_SOURCE,
        theme = "tomorrow_night",
        height = 260,
        fontSize = 16,
        minLines = 8,
        maxLines = 20,
        showLineNumbers = true,
        showGutter = true,
        highlightActiveLine = true,
        displayIndentGuides = true,
        showPrintMargin = false,
        style = Styles("width" => "100%")
    )
    terminal_style = Styles(
        "width" => "100%",
        "height" => "auto",
        "min-height" => "100%",
        "padding" => "0.85rem 1rem",
        "overflow" => "visible",
        "background" => "#171b1b",
        "color" => "#dbdee0",
        "font-family" => "JuliaMono, monospace",
        "font-size" => "0.88rem",
        "line-height" => "1.45",
        "white-space" => "pre-wrap",
        "overflow-wrap" => "anywhere"
    )
    terminal = TerminalOutput(; style = terminal_style)
    run_button = Bonito.Button(
        "Run";
        id = "live-julia-run",
        class = "lc-repl-button is-primary",
        style = nothing,
        title = "Run the editor contents (Ctrl+Enter)",
        ariaLabel = "Run Julia code"
    )
    interrupt_button = Bonito.Button(
        "Interrupt";
        id = "live-julia-interrupt",
        class = "lc-repl-button",
        style = nothing,
        title = "Stop the current evaluation and restart its worker",
        ariaLabel = "Interrupt Julia evaluation"
    )
    restart_button = Bonito.Button(
        "Restart worker";
        id = "live-julia-restart",
        class = "lc-repl-button",
        style = nothing,
        title = "Discard all REPL state and start a fresh Julia worker",
        ariaLabel = "Restart Julia worker"
    )
    clear_button = Bonito.Button(
        "Clear console";
        id = "live-julia-clear",
        class = "lc-repl-button",
        style = nothing,
        title = "Clear the console transcript",
        ariaLabel = "Clear Julia console"
    )
    status = Observable("Starting isolated Julia worker…")
    controller = WorkerController(nothing, 1, false, false, 0)

    try
        controller.worker = start_worker()
        status[] = "Worker PID $(controller.worker.pid) · ready"
        push!(terminal, worker_banner(controller.worker))
    catch error
        status[] = "Worker failed to start"
        push!(terminal, "\e[1;31mWorker startup failed:\e[0m $(sprint(showerror, error))\n")
    end

    toolbar = Grid(
        run_button,
        interrupt_button,
        restart_button,
        clear_button,
        status_line(status; class = "lc-repl-status");
        columns = "repeat(4, max-content) minmax(0, 1fr)",
        gap = "0.5rem",
        align_items = "center",
        class = "lc-repl-toolbar"
    )
    editor_panel = webpart(
        toolbar,
        editor;
        kind = :panel,
        id = "live-julia-editor",
        class = "lc-repl-editor"
    )
    terminal_host = Grid(
        terminal;
        rows = "max-content",
        gap = "0",
        id = "live-julia-terminal",
        class = "lc-repl-terminal",
        style = Styles("min-height" => "0", "overflow" => "auto")
    )
    console_panel = webpart(
        terminal_host;
        kind = :panel,
        title = "Julia console",
        meta = "process-isolated · stateful",
        id = "live-julia-console",
        class = "lc-repl-console"
    )
    canvas = webgrid(
        reshape([:editor, :console], 2, 1);
        rows = "minmax(12rem, 2fr) minmax(11rem, 3fr)",
        height = "100%",
        gap = "0.75rem",
        stack = false,
        class = "lc-repl-grid",
        editor = editor_panel,
        console = console_panel
    )

    function current(worker, epoch)
        return !controller.closed && controller.epoch == epoch &&
               controller.worker === worker
    end

    function replace_worker!(message::AbstractString; force::Bool = controller.busy)
        controller.closed && return nothing
        controller.epoch += 1
        epoch = controller.epoch
        previous = controller.worker
        controller.worker = nothing
        controller.busy = false
        status[] = message
        @async begin
            isnothing(previous) || stop_worker!(previous; force)
            replacement = try
                start_worker()
            catch error
                controller.epoch == epoch || return
                status[] = "Worker restart failed"
                push!(terminal,
                    "\n\e[1;31mWorker restart failed:\e[0m $(sprint(showerror, error))\n")
                return
            end
            if controller.closed || controller.epoch != epoch
                stop_worker!(replacement)
                return
            end
            controller.worker = replacement
            status[] = "Worker PID $(replacement.pid) · ready"
            push!(terminal, "\n$(worker_banner(replacement))")
        end
        return nothing
    end

    run_observer = on(run_button.value) do clicked
        clicked || return nothing
        controller.closed && return nothing
        if controller.busy
            status[] = "Worker is already evaluating · use Interrupt if needed"
            return nothing
        end
        worker = controller.worker
        if isnothing(worker) || !worker_running(worker)
            status[] = "Worker unavailable · restart it before running code"
            return nothing
        end
        source = editor.onchange[]
        isempty(strip(source)) && return nothing
        epoch = controller.epoch
        controller.busy = true
        controller.evaluation += 1
        evaluation = controller.evaluation
        status[] = "Evaluation $evaluation running…"
        push!(terminal, prompt(source, evaluation))
        @async begin
            try
                evaluate!(worker, source) do event
                    current(worker, epoch) && append_event!(terminal, event)
                end
                current(worker, epoch) && (status[] = "Worker PID $(worker.pid) · ready")
            catch error
                if current(worker, epoch)
                    controller.worker = nothing
                    status[] = "Worker stopped unexpectedly · restart required"
                    push!(terminal,
                        "\n\e[1;31mWorker connection closed:\e[0m $(sprint(showerror, error))\n")
                end
            finally
                current(worker, epoch) && (controller.busy = false)
            end
        end
        return nothing
    end

    interrupt_observer = on(interrupt_button.value) do clicked
        clicked || return nothing
        worker = controller.worker
        if controller.busy && !isnothing(worker) && worker_running(worker)
            push!(terminal,
                "\n\e[1;33mInterrupt requested; restarting the worker and clearing its state.\e[0m\n")
            replace_worker!("Interrupting evaluation and restarting worker…"; force = true)
        else
            status[] = "Nothing is currently running"
        end
        return nothing
    end

    restart_observer = on(restart_button.value) do clicked
        clicked || return nothing
        return replace_worker!("Restarting isolated Julia worker…")
    end

    clear_observer = on(clear_button.value) do clicked
        clicked || return nothing
        empty!(terminal)
        worker = controller.worker
        isnothing(worker) || push!(terminal, worker_banner(worker))
        return nothing
    end

    Bonito.onload(
        session,
        editor_panel,
        js"""
(root) => {
    root.addEventListener("keydown", (event) => {
        if ((event.ctrlKey || event.metaKey) && event.key === "Enter") {
            event.preventDefault();
            root.querySelector("#live-julia-run")?.click();
        }
    });
}
"""
    )
    Bonito.onjs(
        session,
        terminal.content,
        js"""
() => {
    const terminalHost = $(terminal_host);
    requestAnimationFrame(() => {
        terminalHost.scrollTop = terminalHost.scrollHeight;
    });
}
"""
    )

    on(session.on_close) do _
        controller.closed = true
        controller.epoch += 1
        worker = controller.worker
        force = controller.busy
        controller.worker = nothing
        isnothing(worker) || stop_worker!(worker; force)
        off(run_observer)
        off(interrupt_observer)
        off(restart_observer)
        off(clear_observer)
        return nothing
    end

    return (;
        body = slide(
            "Live Julia workspace",
            canvas;
            lede = md"""
            Evaluate Julia in a session-owned worker process. Definitions persist until you restart the worker; the Bonito server remains outside the execution process.
            """
        ),
        editor,
        terminal,
        controller,
        run_button,
        interrupt_button,
        restart_button,
        clear_button,
        status
    )
end

const PAGE = page_descriptor(
    id = "live-julia-workspace",
    group = "Developer tools",
    title = "Live Julia workspace",
    order = 40,
    render = true,
    class = "lc-repl-slide",
    build = build
)

end

LiveJuliaPage.PAGE
