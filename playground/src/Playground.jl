"""
    Playground

Render the Quarto project of this folder and serve its site: the landing page at `/` and each presentation at
`/presentations/<name>/`. Run it with `julia --project=playground -m Playground`, with the options `--host`
(default `127.0.0.1`), `--port` (default `8080`) and `--no-render`, which serves the site of the last render.

Each presentation's cells run in a Julia process of their own, a [`Presentation`](@ref). `GET /api/presentations`
answers the state of each presentation's process, `stopped`, `starting` or `running`, and how many times it
started. `POST /api/presentations/<name>/restart` ends the process and starts it again. The HTTP requests and
the WebSocket connections under `/presentations/<name>/live/` go to the presentation's process with that prefix
removed. The first of them starts the process. While the process starts, or after it stopped, they get a page that
says so and reloads once the process runs.
"""
module Playground

using HTTP: HTTP
using JSON: JSON
using Random: RandomDevice
using quarto_jll: quarto

"""
    Presentation(name)

The Julia process that evaluates the cells of the presentation `presentations/<name>/` and serves them on
`127.0.0.1:<port>`. The process runs `PlaygroundWorker` on the presentation's `cells.json`, with the
presentation's folder as its project when the folder holds a `Project.toml`, and the worker's own project
otherwise. The constructor returns once the process listens. The process ends when its input closes:
[`close`](@ref) closes it, and so does the end of this process. Each line that the process writes to its error
stream appears on this one, prefixed by the name. When the process ends before it listens, the constructor fails
with the last 20 lines.

The process gets a secret, drawn for it, in the environment variable `PLAYGROUND_SECRET`. It refuses the requests
whose header `Playground-Secret` does not hold the secret.
"""
struct Presentation
    name::String
    process::Base.Process
    port::Int
    secret::String
end

function Presentation(name::String)
    root = pkgdir(@__MODULE__)
    folder = joinpath(root, "presentations", name)
    worker = joinpath(root, "worker")
    project = isfile(joinpath(folder, "Project.toml")) ? folder : worker
    cells = joinpath(root, "_site", "presentations", name, "cells.json")
    command = `$(Base.julia_cmd()) --startup-file=no --threads=2,1 --project=$project
               -m PlaygroundWorker $cells --prefix /presentations/$name/live`
    secret = bytes2hex(rand(RandomDevice(), UInt8, 32))
    environment = addenv(
        command, "JULIA_LOAD_PATH" => "@:$worker:@stdlib", "PLAYGROUND_SECRET" => secret)
    errors = Pipe()
    process = open(pipeline(environment; stderr = errors), "r+")
    close(errors.in)
    lines = String[]
    forwarding = @async for line in eachline(errors)
        println(stderr, name, ": ", line)
        push!(lines, line)
        length(lines) > 20 && popfirst!(lines)
    end
    line = readline(process)
    if !startswith(line, "ready ")
        wait(process)
        wait(forwarding)
        error("the process of $name ended before it was ready:\n", join(lines, '\n'))
    end
    return Presentation(name, process, parse(Int, line[7:end]), secret)
end

"""
    close(presentation::Presentation)

Close the input of the presentation's process and wait for the process to end. After 10 s, kill it.
"""
function Base.close(presentation::Presentation)
    close(presentation.process.in)
    timedwait(() -> process_exited(presentation.process), 10)
    process_exited(presentation.process) || kill(presentation.process, Base.SIGKILL)
    wait(presentation.process)
    return
end

"""
    isopen(presentation::Presentation)

Whether the presentation's process runs.
"""
Base.isopen(presentation::Presentation) = process_running(presentation.process)

function (@main)(arguments)
    host = "127.0.0.1"
    port = 8080
    render = true
    position = 1
    while position <= length(arguments)
        option = arguments[position]
        if option == "--host"
            host = arguments[position + 1]
            position += 2
        elseif option == "--port"
            port = parse(Int, arguments[position + 1])
            position += 2
        elseif option == "--no-render"
            render = false
            position += 1
        else
            throw(ArgumentError("unknown option $option; the options are --host, --port and --no-render"))
        end
    end
    root = pkgdir(@__MODULE__)
    names = filter(readdir(joinpath(root, "presentations"))) do name
        isfile(joinpath(root, "presentations", name, "index.qmd"))
    end
    for name in names
        occursin(r"^[a-z0-9][a-z0-9-]*$", name) && !occursin("http", name) || error(
            "the presentation '$name' needs a name of lowercase letters, digits and hyphens that starts ",
            "with a letter or a digit and does not contain 'http': Bonito turns the first 'http' of the ",
            "address of its WebSocket into 'ws'")
    end
    render && run(`$(quarto()) render $root -M engine:markdown`)

    started = Dict{String, Task}()
    starts = Dict(name => 0 for name in names)
    guard = ReentrantLock()
    router = HTTP.Router(HTTP.streamhandler(HTTP.fileserver(joinpath(root, "_site"))))
    HTTP.register!(router,
        "GET",
        "/api/presentations",
        HTTP.streamhandler() do _
            states = lock(guard) do
                Dict(map(names) do name
                    if !haskey(started, name)
                        name => Dict("state" => "stopped", "starts" => starts[name])
                    elseif !istaskdone(started[name])
                        name => Dict("state" => "starting", "starts" => starts[name])
                    elseif istaskfailed(started[name])
                        name => Dict("state" => "stopped", "starts" => starts[name],
                            "error" => sprint(showerror, started[name].exception))
                    else
                        name => Dict(
                            "state" => isopen(fetch(started[name])) ? "running" : "stopped",
                            "starts" => starts[name])
                    end
                end)
            end
            return HTTP.Response(200, ["Content-Type" => "application/json"], JSON.json(states))
        end)
    HTTP.register!(router,
        "POST",
        "/api/presentations/{name}/restart",
        HTTP.streamhandler() do request
            name = HTTP.getparam(request, "name")
            name in names || return HTTP.Response(404)
            lock(guard) do
                previous = haskey(started, name) ? [started[name]] : Task[]
                starts[name] += 1
                started[name] = Threads.@spawn begin
                    waitall(previous; failfast = false, throw = false)
                    foreach(task -> istaskfailed(task) || close(fetch(task)), previous)
                    Presentation(name)
                end
            end
            return HTTP.Response(202)
        end)
    HTTP.register!(router, "/presentations/{name}/live/**") do stream
        name = HTTP.getparam(stream.message, "name")
        name in names || return HTTP.streamhandler(_ -> HTTP.Response(404))(stream)
        task = lock(guard) do
            get!(started, name) do
                starts[name] += 1
                Threads.@spawn Presentation(name)
            end
        end
        if !istaskdone(task) || istaskfailed(task) || !isopen(fetch(task))
            # The page joins the deck as a live view does: the lifecycle script sends the deck lcm:child-ready and the
            # class lc-widget-shell marks the page. It polls the state and reloads once the process runs.
            state = istaskdone(task) ? "stopped" : "starting"
            page = """
                <!doctype html>
                <html lang="en">
                <head>
                <meta charset="utf-8">
                <title>Julia $state</title>
                <style>
                $(read(joinpath(root, "assets", "brand.css"), String))
                body { margin: 0; background: var(--lc-bg); color: var(--lc-muted);
                       font: 14px/1.4 Lato, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; }
                .lc-widget-shell { display: grid; place-content: center; min-height: 100vh; padding: 0 1em; }
                pre { color: var(--lc-danger); white-space: pre-wrap; }
                </style>
                </head>
                <body>
                <div class="lc-widget-shell"><p>Julia $state</p><pre></pre></div>
                <script>$(read(joinpath(root, "worker", "src", "lifecycle.js"), String))</script>
                <script>
                setInterval(async () => {
                  try {
                    const presentation = (await (await fetch('/api/presentations')).json())['$name'];
                    if (presentation.state === 'running') location.reload();
                    document.querySelector('p').textContent = 'Julia ' + presentation.state;
                    document.querySelector('pre').textContent = presentation.error || '';
                  } catch (_) { /* The next poll tries again. */ }
                }, 1000);
                </script>
                </body>
                </html>
                """
            return HTTP.streamhandler(_ -> HTTP.Response(
                503, ["Content-Type" => "text/html; charset=utf-8"], page))(stream)
        end
        presentation = fetch(task)
        worker = "127.0.0.1:$(presentation.port)" *
                 chopprefix(stream.message.target, "/presentations/$name/live")
        secret = "Playground-Secret" => presentation.secret
        HTTP.WebSockets.isupgrade(stream.message) ||
            return HTTP.streamhandler(request -> HTTP.request(
                request.method, "http://$worker", [request.headers; secret], request.body;
                decompress = false, redirect = false, retry = false, status_exception = false))(stream)
        HTTP.WebSockets.upgrade(stream; maxframesize = typemax(Int) ÷ 4,
            maxfragmentation = typemax(Int)) do browser
            HTTP.WebSockets.open("ws://$worker"; headers = [secret], redirect = false,
                maxframesize = typemax(Int) ÷ 4, maxfragmentation = typemax(Int)) do process
                directions = map(((browser, process), (process, browser))) do (from, to)
                    Threads.@spawn try
                        for message in from
                            HTTP.WebSockets.send(to, message)
                        end
                        close(to, from.closebody)
                    catch error
                        error isa HTTP.WebSockets.WebSocketError || rethrow()
                        close(to, error.message)
                    end
                end
                waitall(directions; failfast = false, throw = false)
            end
        end
    end
    try
        HTTP.listen(router, host, port)
    finally
        tasks = lock(() -> collect(values(started)), guard)
        waitall(tasks; failfast = false, throw = false)
        foreach(task -> istaskfailed(task) || close(fetch(task)), tasks)
    end
    return 0
end

end
