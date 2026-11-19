"""
    PlaygroundWorker

Evaluate the cells of one presentation and serve each as a Bonito app at `/cells/<n>`. Run it with
`julia --threads=2,1 -m PlaygroundWorker <cells.json> --prefix <path>`. The file `cells.json` is the list that
`deck.lua` writes beside the rendered presentation. The browser loads the pages of the cells and their assets
under the path `--prefix`. Once the server listens, the worker prints `ready <port>` and logs the version of
Bonito it loaded. It then evaluates the cells marked `autorun` in order and returns at the end of its input.

Each browser session that asks for a cell gets its own copy of the cell, with the source of the presentation. The
value of that copy goes to that session alone. A cell marked `autorun` runs again for the session, and the worker
serves its page when that run ends. So a WebGL figure belongs to one session, and a reloaded frame draws its own
figure.

The worker answers 403 to the requests whose header `Playground-Secret` does not hold the value of the
environment variable `PLAYGROUND_SECRET`.

The cells run one at a time in `Main`, with the soft scope of the REPL, on the default thread pool. The page of a
cell shows its source in an editor when `echo` is true, and a button that runs it again. Below them come its
printed output, capped at 1 MiB, and its value or the message and trace of its error. The page fills its frame,
and the value takes the height that the editor, the button and the output leave. A run keeps the previous value
until it ends.
"""
module PlaygroundWorker

using Bonito: Bonito, App, Asset, Button, CodeEditor, DOM, Observable, Session, on, route!
using HTTP: HTTP
using JSON: JSON
using REPL: REPL

struct Cell
    source::Observable{String}
    autorun::Bool
    echo::Bool
    output::Observable{String}
    value::Observable{Any}
    function Cell(source::String, autorun::Bool, echo::Bool)
        return new(
            Observable(source), autorun, echo, Observable(""), Observable{Any}(DOM.span()))
    end
end

function Bonito.jsrender(session::Session, cell::Cell)
    assets = normpath(@__DIR__, "..", "..", "assets")
    height = 24 * (count('\n', cell.source[]) + 2)
    editor = CodeEditor("julia"; initial_source = cell.source[], height)
    button = Button("Run")
    on(button.value) do _
        cell.source[] = editor.onchange[]
    end
    shell = DOM.div(
        Asset(joinpath(assets, "brand.css")),
        Asset(joinpath(assets, "code-theme.css")),
        Asset(joinpath(assets, "control-contract.css")),
        DOM.script(read(joinpath(@__DIR__, "lifecycle.js"), String)),
        DOM.style("body { margin: 0; }"),
        cell.echo ? editor : DOM.span(),
        DOM.div(button),
        DOM.pre(cell.output; style = "margin: 0"),
        DOM.div(cell.value;
            style = "display: grid; grid-template: minmax(0, 1fr) / minmax(0, 1fr); flex: 1 1 0; min-height: 0");
        class = "lc-widget-shell",
        style = "display: flex; flex-direction: column; gap: 8px; height: 100vh; box-sizing: border-box; " *
                "padding: 8px; background: var(--lc-bg); color: var(--lc-text)"
    )
    return Bonito.jsrender(session, shell)
end

function (@main)(arguments)
    prefix = ""
    position = 2
    while position <= length(arguments)
        option = arguments[position]
        if option == "--prefix"
            prefix = arguments[position + 1]
            position += 2
        else
            throw(ArgumentError("unknown option $option; the option is --prefix"))
        end
    end
    cells = [Cell(entry["source"], entry["autorun"], entry["echo"])
             for entry in JSON.parsefile(arguments[1])]
    log = Base.Filesystem.File(Libc.dup(RawFD(2)))
    Base.CoreLogging.global_logger(Base.CoreLogging.ConsoleLogger(log))
    secret = ENV["PLAYGROUND_SECRET"]
    server = Bonito.Server("127.0.0.1", 0; proxy_url = prefix,
        handler = context -> HTTP.header(context.request, "Playground-Secret") == secret ?
                             Bonito.HTTPServer.delegate(context) : HTTP.Response(403))
    queue = Channel{Pair{Int, Cell}}(Inf)
    for (index, cell) in enumerate(cells)
        on(cell.source) do _
            put!(queue, index => cell)
        end
        route!(server,
            "/cells/$index" => App(; title = "Julia cell $index") do
                shown = Cell(cell.source[], cell.autorun, cell.echo)
                on(_ -> put!(queue, index => shown), shown.source)
                if shown.autorun
                    ran = Base.Event()
                    on(_ -> notify(ran), shown.value)
                    notify(shown.source)
                    wait(ran)
                end
                return shown
            end)
    end
    @info "Bonito $(pkgversion(Bonito))"
    println("ready ", server.port)
    flush(stdout)
    Threads.@spawn :default for (index, cell) in queue
        cell.output[] = ""
        pipe = Pipe()
        forwarding = @async begin
            printed = 0
            while true
                bytes = readavailable(pipe)
                count = min(length(bytes), 1_048_576 - printed)
                count > 0 && (cell.output[] *= String(bytes[1:count]))
                count < length(bytes) && printed <= 1_048_576 &&
                    (cell.output[] *= "\n[output truncated after 1 MiB]\n")
                printed += length(bytes)
                eof(pipe) && break
            end
        end
        try
            result = redirect_stdout(pipe) do
                redirect_stderr(pipe) do
                    include_string(REPL.softscope, Main, cell.source[], "cell $index")
                end
            end
            cell.value[] = endswith(rstrip(cell.source[]), ";") ? DOM.span() : result
        catch error
            trace = stacktrace(catch_backtrace())
            last = findlast(frame -> String(frame.file) == "cell $index", trace)
            text = sprint(
                showerror, error, last === nothing ? empty(trace) : trace[1:last];
                context = :limit => true)
            cell.value[] = DOM.pre(text; style = "color: var(--lc-danger)")
        finally
            close(pipe.in)
            wait(forwarding)
            close(pipe.out)
        end
    end
    for cell in cells
        cell.autorun && notify(cell.source)
    end
    read(stdin)
    close(server)
    return 0
end

end
