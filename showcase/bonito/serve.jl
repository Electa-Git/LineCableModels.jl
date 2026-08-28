using LineCableModels
using Bonito
using WGLMakie

WGLMakie.activate!()
include(joinpath(@__DIR__, "app.jl"))

host = get(ENV, "HOST", "127.0.0.1")
port = parse(Int, get(ENV, "PORT", "8080"))
proxy_url = get(ENV, "PROXY_URL", ".")
server = Bonito.Server(cable_app(), host, port; proxy_url)
println("Showcase listening at $(Bonito.online_url(server, "/"))")

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
