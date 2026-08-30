using LineCableModels
import PowerImpedance
using Bonito
using GraphMakie
using Graphs
using NetworkLayout
using WGLMakie

WGLMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

host = get(ENV, "HOST", "127.0.0.1")
port = parse(Int, get(ENV, "PORT", "8080"))
proxy_url = get(ENV, "PROXY_URL", ".")
server = Bonito.Server(host, port; proxy_url)
Bonito.route!(server, cable_routes())
println("Showcase listening at $(Bonito.online_url(server, "/"))")
start_page_preparations!()

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
