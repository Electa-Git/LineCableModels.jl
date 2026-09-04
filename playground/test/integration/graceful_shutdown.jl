using Downloads
using Sockets

function available_port()
    socket = listen(ip"127.0.0.1", 0)
    port = Int(getsockname(socket)[2])
    close(socket)
    return port
end

function responds(url)
    try
        response = Downloads.request(url; timeout=3)
        return response.status == 200
    catch
        return false
    end
end

function await_response(url; timeout_seconds=15)
    return timedwait(
        () -> responds(url),
        timeout_seconds;
        pollint=0.1
    ) == :ok
end

function authored_page_routes(site_directory)
    routes = String[]
    for (root, _, files) in walkdir(site_directory)
        "index.html" in files || continue
        relative = replace(relpath(root, site_directory), '\\' => '/')
        push!(routes, relative == "." ? "/" : "/$relative/")
    end
    return sort!(unique!(routes))
end

const LIVE_WIDGET_ROUTES = [
    "/widgets/slider",
    "/widgets/toggle",
    "/widgets/dropdown",
    "/widgets/text-input",
    "/widgets/number-spinner",
    "/widgets/actions",
    "/widgets/progress",
    "/widgets/console",
    "/widgets/tabs",
    "/widgets/toolbar",
    "/widgets/control-panel",
    "/widgets/form-toolkit",
    "/widgets/overlay-toolkit",
    "/widgets/data-view-toolkit",
    "/widgets/repeater",
    "/widgets/job-panel",
]

script_directory = @__DIR__
playground_directory = normpath(joinpath(script_directory, "..", ".."))
launcher = joinpath(playground_directory, "lcm")
isfile(launcher) || error("graceful shutdown test requires launcher: $launcher")

mktempdir(prefix="lcm-shutdown-") do directory
    publisher_port = available_port()
    broker_port = available_port()
    log_path = joinpath(directory, "publisher.log")
    log = open(log_path, "w")
    command = addenv(
        `$(launcher) playground start --no-open --no-render --host 127.0.0.1 --port $(publisher_port)`,
        "NATS_CONNECT_URL" => "nats://127.0.0.1:$broker_port"
    )
    process = run(pipeline(command; stdout=log, stderr=log); wait=false)
    try
        readiness = timedwait(
            () -> !process_running(process) ||
                  responds("http://127.0.0.1:$publisher_port/"),
            30;
            pollint=0.1
        )
        readiness == :ok || error("publisher did not become ready")
        process_running(process) || error("publisher exited before becoming ready")
        responds("http://127.0.0.1:$publisher_port/") || error(
            "publisher did not serve its static home while NATS was offline"
        )
        routes = vcat(
            authored_page_routes(joinpath(playground_directory, "_site")),
            LIVE_WIDGET_ROUTES
        )
        for route in routes
            await_response("http://127.0.0.1:$publisher_port$route") || error(
                "publisher did not serve $route while NATS was offline"
            )
        end

        kill(process, Base.SIGINT)
        stopped = timedwait(() -> !process_running(process), 20; pollint=0.1)
        stopped == :ok || error("publisher did not stop after SIGINT")
        try
            wait(process)
        catch
        end
        success(process) || error("publisher exited with status $(process.exitcode)")
    finally
        if process_running(process)
            kill(process, Base.SIGTERM)
            timedwait(() -> !process_running(process), 5; pollint=0.1)
        end
        process_running(process) && kill(process, Base.SIGKILL)
        close(log)
    end

    output = read(log_path, String)
    occursin("Stopping LineCableModels playground", output) || error(
        "launcher did not report orderly shutdown:\n$output"
    )
    for forbidden in (
        "error while cleaning up server",
        "InterruptException",
        "Cannot flush send buffer",
        "Connection disconnected after",
    )
        occursin(forbidden, output) && error(
            "shutdown emitted `$forbidden`:\n$output"
        )
    end
end

println("Offline publisher startup and graceful SIGINT shutdown passed")
