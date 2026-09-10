module LineCableModelsPlayground

using AWS
using AWSS3
using Bonito
using BonitoWidgets
using Dates
using EzXML
using JSON3
using LineCableModelsPlaygroundProtocol
using NATS
using Observables
using Quarto
using SHA
using URIs
using UUIDs
using ZipFile

const PLAYGROUND_ROOT = normpath(joinpath(@__DIR__, ".."))
const SITE_DIR = joinpath(PLAYGROUND_ROOT, "_site")
const DEFAULT_HOST = "127.0.0.1"
const DEFAULT_PORT = 8080
include(joinpath(PLAYGROUND_ROOT, "common", "PublishedAssets.jl"))
const BRAND_THEME_PATH = joinpath(PLAYGROUND_ROOT, "assets", "brand.css")
include_dependency(BRAND_THEME_PATH)
const BRAND_THEME = read(BRAND_THEME_PATH, String)
const CONTROL_CONTRACT_PATH = joinpath(
    PLAYGROUND_ROOT,
    "assets",
    "control-contract.css"
)
include_dependency(CONTROL_CONTRACT_PATH)
const CONTROL_CONTRACT = read(CONTROL_CONTRACT_PATH, String)
const THEME_INIT_PATH = joinpath(PLAYGROUND_ROOT, "assets", "theme-init.html")
include_dependency(THEME_INIT_PATH)
const THEME_INITIALIZER = String(only(match(r"(?s)^<script>\s*(.*?)\s*</script>\s*$", read(THEME_INIT_PATH, String)).captures))
# Pure browser initialization must not wait for a Bonito session/WebSocket.
# This is trusted repository source, with no interpolation of application data.
theme_script() = DOM.script(THEME_INITIALIZER)

include("diagnostics/ComponentXRay.jl")
include("toolkit/Toolkit.jl")
using .Toolkit
include("workbench/WorkbenchUI.jl")
include("toolbar.jl")
include("ribbon.jl")
include("repeater.jl")
include("uploads.jl")
include("geographic_map.jl")
include("power_system_canvas.jl")
include("console.jl")
include("artifacts.jl")
include("container_runtime.jl")
include("presentation.jl")
include("broker/Subjects.jl")
include("broker/ConnectionState.jl")
include("broker/JobHandle.jl")
include("broker/BrokerClient.jl")
include("broker/RuntimeAdmin.jl")
include("widgets/WorkerStatus.jl")
include("widgets/JobControls.jl")
include("widgets/RuntimeControls.jl")
include("widgets/JuliaTerminal.jl")
include("scientific/ScientificViews.jl")
include("widgets.jl")
include("workbenches/TemplateWorkbench.jl")
include("workbenches/CableStudy.jl")
include("applications/Showcase.jl")
include("applications/UIHost.jl")
include("applications/Catalogue.jl")

export ConsoleEntry,
    serve_owned_ui,
    ConsoleView,
    ComponentXRay,
    ComboBox,
    ConfirmDialog,
    DataTable,
    Dialog,
    DialogAction,
    DialogEvent,
    Disclosure,
    Field,
    Form,
    FormDialog,
    InlineNotice,
    MessageDialog,
    MultiSelect,
    PropertyGrid,
    PropertyItem,
    RadioGroup,
    RangeInput,
    SecretInput,
    SegmentedControl,
    TableColumn,
    TextAreaInput,
    TextInput,
    ToastCenter,
    ToastEntry,
    Toolkit,
    UnitNumberInput,
    ViewportFrame,
    WorkspacePage,
    ArtifactGateway,
    AbstractUploadStore,
    BrokerClient,
    DirectoryUploadStore,
    GeographicMap,
    JobHandle,
    JobPanel,
    KMLUploadSource,
    MapReference,
    PowerSystemCanvas,
    RepeatEntry,
    Repeater,
    Ribbon,
    RibbonGroup,
    RibbonTab,
    UploadedFile,
    UploadField,
    UploadPolicy,
    UploadRegistry,
    UploadTarget,
    TemplateWorkbench,
    Toolbar,
    ToolbarBinding,
    ToolbarButton,
    ToolbarDropdown,
    ToolbarEvent,
    ToolbarNumber,
    ToolbarSeparator,
    ToolbarToggle,
    WorkerStatus,
    RuntimeClient,
    JuliaTerminal,
    WorkerSelector,
    PreparationStatus,
    ScientificJob,
    ScientificResult,
    WorkerDiagnostics,
    WorkerControlPanel,
    WorkbenchUI,
    append_console!,
    add!,
    cancel!,
    clear_console!,
    clear_toasts!,
    close!,
    close_dialog!,
    dismiss_notice!,
    dismiss_toast!,
    field_error!,
    mark_dirty!,
    move!,
    open_dialog!,
    push_toast!,
    reattach!,
    remove!,
    reset_form!,
    set_viewport_state!,
    set_console_status!,
    start!,
    submit!,
    submit_async!,
    duplicate!,
    delete_upload!,
    register_upload_route!,
    snapshot,
    store_upload!,
    submit_payload,
    sweep_upload_staging!,
    upload_path,
    validate_kml_upload,
    validate_power_system_diagram,
    validate_power_system_topology,
    toolbar,
    toolbar_event_name,
    toolbar_icon

const WORKBENCH_ROUTES = (
    "/workbenches/template" => TemplateWorkbench.app,
)

function register_workbench_routes!(server; xray::Bool=false)
    for (route, factory) in WORKBENCH_ROUTES
        Bonito.route!(server, route => factory(; xray))
    end
    return server
end

struct StaticFileHandler
    site::PublishedAssets.PublishedSite
end

# A publication fallback must remain behind routes installed lazily by Bonito
# (notably its content-addressed asset server). A catch-all Regex steals those
# requests when the first interactive page installs its assets after startup.
struct PublishedRoute
    site::PublishedAssets.PublishedSite
end
Bonito.HTTPServer.pattern_priority(::PublishedRoute) = 4
function Bonito.HTTPServer.match_request(route::PublishedRoute, request)
    path = String(URIs.URI(request.target).path)
    return PublishedAssets.published_asset(route.site, path) === nothing ? nothing : path
end

function Bonito.HTTPServer.apply_handler(file::StaticFileHandler, context)
    path = String(URIs.URI(context.request.target).path)
    body = PublishedAssets.published_asset(file.site, path)
    body === nothing && return Bonito.HTTP.Response(404, "Published page not found")
    headers = [
        "Content-Type" => Bonito.file_mimetype(endswith(path, "/") || !occursin('.', basename(path)) ? "index.html" : path),
        "Cache-Control" => "no-store",
        "X-Content-Type-Options" => "nosniff",
    ]
    return Bonito.HTTP.Response(200, headers; body)
end

function register_static_site_routes!(server, directory=SITE_DIR)
    isdir(directory) || throw(ArgumentError("static site directory not found: $directory"))

    handler = StaticFileHandler(PublishedAssets.PublishedSite(directory))
    for route in keys(handler.site.files)
        Bonito.route!(server, route => handler)
    end
    # Exact live routes retain precedence; newly committed static routes and
    # content-hashed assets are resolved by the same validated snapshot.
    Bonito.route!(server, PublishedRoute(handler.site) => handler)
    return server
end

function usage(io::IO=stdout; feature::Union{Nothing,String}=nothing)
    if feature == "playground"
        println(io, "Usage:")
        println(io, "  lcm playground build [--quiet]")
        println(io, "  lcm playground start [options]")
        println(io)
        println(io, "Start options:")
        println(io, "  --host HOST         Bind address (default: HOST or 127.0.0.1)")
        println(io, "  --port PORT         Listening port (default: PORT or 8080)")
        println(io, "  --proxy-url URL     Public Bonito URL (default: PROXY_URL or auto-detected)")
        println(io, "  --no-render         Serve the existing Quarto build")
        println(io, "  --no-open           Do not open the default browser")
        println(io, "  --xray             Enable owned-component diagnostics")
        println(io, "  -h, --help          Show this help")
    elseif feature == "nats"
        println(io, "Usage:")
        println(io, "  lcm nats init")
        println(io, "  lcm nats status")
        println(io)
        println(io, "NATS_CONNECT_URL selects the client endpoint.")
        println(io, "NATS_ADMIN_URL may select a separate administrator endpoint.")
    elseif feature == "container"
        println(io, "Usage:")
        println(io, "  lcm container resolve [--runtime auto|docker|podman]")
        println(io, "  lcm container start [options]")
        println(io, "  lcm container status [options]")
        println(io, "  lcm container logs [options] [SERVICE ...]")
        println(io, "  lcm container stop [options]")
        println(io)
        println(io, "Common options:")
        println(io, "  --runtime RUNTIME   auto (default), docker, or podman")
        println(io, "  --remote            Use the mTLS/S3 remote deployment profile")
        println(io)
        println(io, "Start options:")
        println(io, "  --no-build          Start from existing images")
        println(io, "  --cpu-limits        Add the optional CPU-quota override")
        println(io)
        println(io, "Log options:")
        println(io, "  --follow            Follow log output")
        println(io, "  --tail N            Show the last N lines")
        println(io)
        println(io, "Stop options:")
        println(io, "  --volumes           Also remove persistent deployment volumes")
        println(io)
        println(io, "LCM_CONTAINER_RUNTIME sets the default runtime selection.")
    elseif feature == "presentation"
        println(io, "Usage:")
        println(io, "  lcm presentation build [SOURCE] [--quiet]")
        println(io, "  lcm presentation start [SOURCE] [options]")
        println(io, "  lcm presentation check [SOURCE] [--quiet]")
        println(io, "  lcm presentation export [SOURCE] --pdf [--output FILE] [--quiet]")
        println(io)
        println(io, "SOURCE defaults to presentations/specimen.qmd.")
        println(io)
        println(io, "Start options:")
        println(io, "  --host HOST         Bind address (default: HOST or 127.0.0.1)")
        println(io, "  --port PORT         Listening port (default: PORT or 8080)")
        println(io, "  --proxy-url URL     Public Bonito URL (default: PROXY_URL or auto-detected)")
        println(io, "  --no-render         Serve an existing presentation build")
        println(io, "  --no-open           Do not open the default browser")
        println(io, "  -h, --help          Show this help")
        println(io)
        println(io, "LCM_BROWSER selects Chrome/Chromium for PDF export.")
    else
        println(io, "LineCableModels playground")
        println(io)
        println(io, "Usage: lcm <feature> <action> [options]")
        println(io)
        println(io, "Features:")
        println(io, "  playground  Build or serve the Quarto/Bonito publisher")
        println(io, "  presentation Build, check, serve, or export reusable decks")
        println(io, "  worker      Start an isolated scientific worker")
        println(io, "  nats        Initialize or inspect the JetStream runtime")
        println(io, "  container   Run the isolated stack with Docker or Podman")
        println(io, "  runtime     Check or start the owned application gateway")
        println(io, "  demo        Start, stop, or inspect the private Kubuntu demo")
        println(io)
        println(io, "Run `lcm <feature> --help` for feature-specific usage.")
    end
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
    xray = false

    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument in ("-h", "--help") && return nothing
        if argument == "--no-open"
            open_browser = false
        elseif argument == "--no-render"
            render_before_start = false
        elseif argument == "--xray"
            xray = true
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
    return (; host, port, proxy_url, open_browser, render_before_start, xray)
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
    write(joinpath(SITE_DIR, "assets", "application-catalogue.json"),
        JSON3.write(ApplicationCatalogue.public_entries()))
    PublishedAssets.publish!(SITE_DIR)
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
        render_before_start,
        xray,
        render_target=nothing,
        initial_path="/",
        broker_enabled=true
    )
    if render_before_start
        isnothing(render_target) ? render_site() : render_presentation(render_target)
    elseif !isfile(joinpath(SITE_DIR, "index.html"))
        throw(ArgumentError(
            "No rendered site exists. Run `lcm playground build` first."
        ))
    end

    broker = broker_enabled ? BrokerClient() : nothing
    server = try
        # Bonito.Server opens its listener immediately. Build the complete route
        # table first: an early successful home request must not expose a site
        # whose deeper pages or parser-blocking script assets are still missing.
        routes = Bonito.HTTPServer.Routes()
        register_workbench_routes!(routes; xray)
        register_widget_routes!(routes, broker)
        register_presentation_routes!(routes)
        register_artifact_route!(routes, default_artifact_gateway())
        register_upload_route!(routes)
        register_static_site_routes!(routes)
        Bonito.Server(host, port; proxy_url, routes)
    catch
        isnothing(broker) || close!(broker)
        rethrow()
    end
    url = Bonito.online_url(server, initial_path)
    println("LineCableModels playground listening at $url")
    open_browser && open_default_browser(url)

    server_only_shutdown = server_shutdown(server)
    shutdown = let closed=Ref(false)
        () -> begin
            closed[] && return nothing
            closed[] = true
            isnothing(broker) || close!(broker)
            server_only_shutdown()
            return nothing
        end
    end
    atexit(shutdown)

    supervised = get(ENV, "LCM_SUPERVISED", "") == "1"
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

function parse_presentation_export_options(arguments)
    pdf = false
    quiet = false
    output = nothing
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument in ("-h", "--help") && return nothing
        if argument == "--pdf"
            pdf = true
        elseif argument == "--quiet"
            quiet = true
        elseif argument == "--output" || startswith(argument, "--output=")
            output, index = option_value(arguments, index, "--output")
        else
            throw(ArgumentError("unknown presentation export option: $argument"))
        end
        index += 1
    end
    pdf || throw(ArgumentError("presentation export currently requires --pdf"))
    return (; output, quiet)
end

function run_cli(arguments)
    isempty(arguments) && return usage()
    arguments in (["-h"], ["--help"]) && return usage()

    feature = arguments[1]
    if feature == "nats"
        length(arguments) == 1 && return usage(; feature="nats")
        any(argument -> argument in ("-h", "--help"), arguments[2:end]) &&
            return usage(; feature="nats")
        arguments == ["nats", "init"] && return initialize_runtime()
        arguments == ["nats", "status"] && return runtime_status()
        throw(ArgumentError("unknown nats action: $(arguments[2])"))
    elseif feature == "container"
        length(arguments) == 1 && return usage(; feature="container")
        arguments[2] in ("-h", "--help") && return usage(; feature="container")
        action = arguments[2]
        action in ("resolve", "start", "status", "logs", "stop") ||
            throw(ArgumentError("unknown container action: $action"))
        any(argument -> argument in ("-h", "--help"), arguments[3:end]) &&
            return usage(; feature="container")
        options = parse_container_options(action, arguments[3:end])
        return run_container_action(action, options)
    elseif feature == "presentation"
        length(arguments) == 1 && return usage(; feature="presentation")
        arguments[2] in ("-h", "--help") && return usage(; feature="presentation")
        action = arguments[2]
        action in ("build", "start", "check", "export") ||
            throw(ArgumentError("unknown presentation action: $action"))
        source, options = split_presentation_arguments(arguments[3:end])
        any(option -> option in ("-h", "--help"), options) &&
            return usage(; feature="presentation")

        if action == "build"
            all(==("--quiet"), options) || throw(ArgumentError(
                "unknown presentation build option; only --quiet is supported"
            ))
            render_presentation(source; quiet="--quiet" in options)
            return nothing
        elseif action == "check"
            all(==("--quiet"), options) || throw(ArgumentError(
                "unknown presentation check option; only --quiet is supported"
            ))
            check_presentation(source; quiet="--quiet" in options)
            return nothing
        elseif action == "export"
            export_options = parse_presentation_export_options(options)
            isnothing(export_options) && return usage(; feature="presentation")
            export_presentation(source; export_options...)
            return nothing
        end

        start_options = parse_start_options(options)
        isnothing(start_options) && return usage(; feature="presentation")
        html = presentation_output(source)
        route = "/" * replace(relpath(html, SITE_DIR), '\\' => '/')
        return start_server(;
            start_options...,
            render_target=source,
            initial_path=route,
            broker_enabled=false
        )
    elseif feature != "playground"
        throw(ArgumentError("unknown feature: $feature"))
    end

    length(arguments) == 1 && return usage(; feature="playground")
    arguments[2] in ("-h", "--help") && return usage(; feature="playground")

    command = arguments[2]
    if command == "build"
        options = arguments[3:end]
        any(option -> option in ("-h", "--help"), options) &&
            return usage(; feature="playground")
        all(==("--quiet"), options) || throw(ArgumentError(
            "unknown build option; only --quiet is supported"
        ))
        render_site(; quiet = "--quiet" in options)
        return nothing
    elseif command == "start"
        options = parse_start_options(arguments[3:end])
        isnothing(options) && return usage(; feature="playground")
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
        println(stderr, "Run `lcm --help` for usage.")
        return 2
    end
end

include("Precompile.jl")

end
