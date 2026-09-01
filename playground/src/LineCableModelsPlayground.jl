module LineCableModelsPlayground

using Bonito
using BonitoWidgets
using Quarto

const PLAYGROUND_ROOT = normpath(joinpath(@__DIR__, ".."))
const SITE_DIR = joinpath(PLAYGROUND_ROOT, "_site")
const DEFAULT_HOST = "127.0.0.1"
const DEFAULT_PORT = 8080

include("toolbar.jl")
include("console.jl")
include("widgets.jl")

export ConsoleEntry,
    ConsoleView,
    Toolbar,
    ToolbarBinding,
    ToolbarButton,
    ToolbarDropdown,
    ToolbarEvent,
    ToolbarSeparator,
    append_console!,
    clear_console!,
    set_console_status!,
    toolbar,
    toolbar_event_name,
    toolbar_icon

struct StaticFileHandler
    path::String
end

function Bonito.HTTPServer.apply_handler(file::StaticFileHandler, context)
    headers = [
        "Access-Control-Allow-Origin" => "*",
        "Content-Type" => Bonito.file_mimetype(file.path),
    ]
    return Bonito.HTTP.Response(200, headers; body=read(file.path))
end

function register_static_site_routes!(server, directory=SITE_DIR)
    isdir(directory) || throw(ArgumentError("static site directory not found: $directory"))

    for (root, _, files) in walkdir(directory)
        for name in files
            path = joinpath(root, name)
            relative = replace(relpath(path, directory), '\\' => '/')
            route = "/$relative"
            handler = StaticFileHandler(path)
            Bonito.route!(server, route => handler)

            if name == "index.html"
                directory_route = dirname(route)
                directory_route == "/" || Bonito.route!(server, directory_route => handler)
                Bonito.route!(
                    server,
                    (directory_route == "/" ? "/" : "$directory_route/") => handler
                )
            end
        end
    end
    return server
end

function usage(io::IO = stdout)
    println(io, "LineCableModels playground")
    println(io)
    println(io, "Usage:")
    println(io, "  linecablemodels playground build [--quiet]")
    println(io, "  linecablemodels playground start [options]")
    println(io)
    println(io, "Start options:")
    println(io, "  --host HOST         Bind address (default: HOST or 127.0.0.1)")
    println(io, "  --port PORT         Listening port (default: PORT or 8080)")
    println(io, "  --proxy-url URL     Public Bonito URL (default: PROXY_URL or auto-detected)")
    println(io, "  --no-render         Serve the existing Quarto build")
    println(io, "  --no-open           Do not open the default browser")
    println(io, "  -h, --help          Show this help")
    return nothing
end

function inferred_proxy_url(port)
    codespace = get(ENV, "CODESPACE_NAME", "")
    forwarding_domain = get(ENV, "GITHUB_CODESPACES_PORT_FORWARDING_DOMAIN", "")
    if !isempty(codespace) && !isempty(forwarding_domain)
        return "https://$codespace-$port.$forwarding_domain"
    end
    return ""
end

function option_value(arguments, index, name)
    argument = arguments[index]
    prefix = "$name="
    startswith(argument, prefix) && return argument[(length(prefix) + 1):end], index
    argument == name || return nothing, index
    index == length(arguments) && throw(ArgumentError("$name requires a value"))
    return arguments[index + 1], index + 1
end

function parse_start_options(arguments)
    host = get(ENV, "HOST", DEFAULT_HOST)
    port = tryparse(Int, get(ENV, "PORT", string(DEFAULT_PORT)))
    isnothing(port) && throw(ArgumentError("PORT must be an integer"))
    proxy_url = get(ENV, "PROXY_URL", nothing)
    open_browser = true
    render_before_start = true

    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument in ("-h", "--help") && return nothing
        if argument == "--no-open"
            open_browser = false
        elseif argument == "--no-render"
            render_before_start = false
        elseif argument == "--host" || startswith(argument, "--host=")
            host, index = option_value(arguments, index, "--host")
        elseif argument == "--port" || startswith(argument, "--port=")
            value, index = option_value(arguments, index, "--port")
            port = tryparse(Int, value)
            isnothing(port) && throw(ArgumentError("--port must be an integer"))
        elseif argument == "--proxy-url" || startswith(argument, "--proxy-url=")
            proxy_url, index = option_value(arguments, index, "--proxy-url")
        else
            throw(ArgumentError("unknown option: $argument"))
        end
        index += 1
    end

    0 < port <= 65535 || throw(ArgumentError("port must be between 1 and 65535"))
    if isnothing(proxy_url) || proxy_url in (".", "./")
        proxy_url = inferred_proxy_url(port)
    end
    return (; host, port, proxy_url, open_browser, render_before_start)
end

function require_quarto()
    executable = Quarto.path()
    if isnothing(executable)
        data_home = get(
            ENV,
            "XDG_DATA_HOME",
            joinpath(homedir(), ".local", "share")
        )
        managed_executable = joinpath(
            data_home,
            "linecablemodels-playground",
            "quarto",
            "current",
            "bin",
            "quarto"
        )
        if isfile(managed_executable)
            ENV["QUARTO_PATH"] = managed_executable
            executable = managed_executable
        end
    end
    isnothing(executable) && throw(ArgumentError(
        "Quarto CLI not found. Run `./playground/bootstrap.sh` or set QUARTO_PATH."
    ))
    return executable
end

function render_site(; quiet = false)
    executable = require_quarto()
    println("Rendering playground with Quarto ($executable)")
    cd(PLAYGROUND_ROOT) do
        Quarto.render(PLAYGROUND_ROOT; execute = false, quiet)
    end
    index = joinpath(SITE_DIR, "index.html")
    isfile(index) || error("Quarto did not produce $index")
    println("Rendered $index")
    return index
end

function browser_command(url)
    if Sys.islinux()
        opener = Sys.which("xdg-open")
        return isnothing(opener) ? nothing : `$opener $url`
    elseif Sys.isapple()
        opener = Sys.which("open")
        return isnothing(opener) ? nothing : `$opener $url`
    elseif Sys.iswindows()
        return `cmd /c start $url`
    end
    return nothing
end

function open_default_browser(url)
    command = browser_command(url)
    if isnothing(command)
        @warn "No default-browser launcher was found" url
        return nothing
    end
    try
        return run(command; wait = false)
    catch error
        @warn "Could not open the default browser" url exception = (
            error,
            catch_backtrace()
        )
        return nothing
    end
end

function server_shutdown(server)
    closed = Ref(false)
    return function ()
        closed[] && return nothing
        closed[] = true
        try
            close(server)
        catch error
            @debug "Error while closing playground server" exception = (
                error,
                catch_backtrace()
            )
        end
        return nothing
    end
end

function shutdown_on_stdin(shutdown)
    return errormonitor(@async begin
        try
            read(stdin, UInt8)
        catch error
            error isa EOFError || rethrow()
        end
        shutdown()
    end)
end

function start_server(;
        host,
        port,
        proxy_url,
        open_browser,
        render_before_start
    )
    if render_before_start
        render_site()
    elseif !isfile(joinpath(SITE_DIR, "index.html"))
        throw(ArgumentError(
            "No rendered site exists. Run `linecablemodels playground build` first."
        ))
    end

    server = Bonito.Server(host, port; proxy_url)
    register_widget_routes!(server)
    register_static_site_routes!(server)
    url = Bonito.online_url(server, "/")
    println("LineCableModels playground listening at $url")
    open_browser && open_default_browser(url)

    shutdown = server_shutdown(server)
    atexit(shutdown)

    supervised = get(ENV, "LINECABLEMODELS_SUPERVISED", "") == "1"
    if supervised
        shutdown_on_stdin(shutdown)
    else
        # Direct invocation remains fail-safe: terminate instead of allowing a
        # Bonito background task to swallow a catchable InterruptException.
        Base.exit_on_sigint(true)
    end
    try
        wait(server)
    finally
        shutdown()
    end
    return nothing
end

function run_cli(arguments)
    arguments in (["-h"], ["--help"]) && return usage()
    length(arguments) >= 2 || throw(ArgumentError(
        "expected `linecablemodels playground build` or `playground start`"
    ))
    arguments[1] == "playground" || throw(ArgumentError(
        "expected `linecablemodels playground build` or `playground start`"
    ))

    command = arguments[2]
    if command == "build"
        options = arguments[3:end]
        any(option -> option in ("-h", "--help"), options) && return usage()
        all(==("--quiet"), options) || throw(ArgumentError(
            "unknown build option; only --quiet is supported"
        ))
        render_site(; quiet = "--quiet" in options)
        return nothing
    elseif command == "start"
        options = parse_start_options(arguments[3:end])
        isnothing(options) && return usage()
        return start_server(; options...)
    end

    throw(ArgumentError("unknown playground command: $command"))
end

function (@main)(arguments)
    try
        return run_cli(arguments)
    catch error
        error isa InterruptException && return nothing
        println(stderr, "error: ", sprint(showerror, error))
        println(stderr, "Run `linecablemodels --help` for usage.")
        return 2
    end
end

end
