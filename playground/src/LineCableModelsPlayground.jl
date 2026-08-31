module LineCableModelsPlayground

using Bonito
using Quarto

const PLAYGROUND_ROOT = normpath(joinpath(@__DIR__, ".."))
const SITE_DIR = joinpath(PLAYGROUND_ROOT, "_site")
const DEFAULT_HOST = "127.0.0.1"
const DEFAULT_PORT = 8080
const DEFAULT_PROXY_URL = "."

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
    println(io, "  --proxy-url URL     Bonito proxy URL (default: PROXY_URL or .)")
    println(io, "  --no-render         Serve the existing Quarto build")
    println(io, "  --no-open           Do not open the default browser")
    println(io, "  -h, --help          Show this help")
    return nothing
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
    proxy_url = get(ENV, "PROXY_URL", DEFAULT_PROXY_URL)
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
    Bonito.route!(server, r".*" => Bonito.FolderServer(SITE_DIR))
    url = Bonito.online_url(server, "/")
    println("LineCableModels playground listening at $url")
    open_browser && open_default_browser(url)

    Base.exit_on_sigint(false)
    try
        wait(server)
    catch error
        error isa InterruptException || rethrow()
    finally
        close(server)
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
