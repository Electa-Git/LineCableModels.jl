using LineCableModels
import PowerImpedance
using Bonito
using GraphMakie
using Graphs
using NetworkLayout
using WGLMakie

WGLMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

function browser_target(arguments)
    open_index = findfirst(==("--open"), arguments)
    isnothing(open_index) && return nothing
    target_index = open_index + 1
    target = target_index <= length(arguments) &&
             !startswith(arguments[target_index], "--") ?
             arguments[target_index] : "/#overview"
    startswith(target, "/") || (target = "/$target")
    return target
end

function open_default_browser(url::AbstractString)
    opener = Sys.which("xdg-open")
    if isnothing(opener)
        @warn "Cannot open the default browser because xdg-open is unavailable" url
        return nothing
    end
    try
        return run(`$opener $url`; wait = false)
    catch error
        @warn "Could not open the default browser" url exception = (
            error,
            catch_backtrace()
        )
        return nothing
    end
end

host = get(ENV, "HOST", "127.0.0.1")
port = parse(Int, get(ENV, "PORT", "8080"))
proxy_url = get(ENV, "PROXY_URL", ".")
server = Bonito.Server(host, port; proxy_url)
Bonito.route!(server, cable_routes())
println("Showcase listening at $(Bonito.online_url(server, "/"))")
target = browser_target(ARGS)
isnothing(target) || open_default_browser(Bonito.online_url(server, target))
start_deck_preparations!()

Base.exit_on_sigint(false)
server_task = @async wait(server)
try
    while !istaskdone(server_task)
        sleep(0.25)
    end
    wait(server_task)
catch exception
    exception isa InterruptException || rethrow()
finally
    close(server)
end
