# Approved mock child: real HTTP process, no Bonito or numerical packages.
ccall(:getppid, Cint, ()) == parse(Int, ENV["LCM_UI_PARENT_PID"]) || exit(70)
using HTTP, JSON3
run_id = ENV["LCM_RUN_ID"]
key = ENV["LCM_UI_HOST_KEY"]
server = HTTP.listen!("127.0.0.1", 0; verbose=-1) do stream
    request = stream.message
    authorized = HTTP.header(request, "X-LCM-Host-Key", "") == key
    if authorized && request.target == "/echo" && HTTP.WebSockets.isupgrade(request)
        HTTP.WebSockets.upgrade(stream; check_origin=(_...) -> true) do socket
            for message in socket
                HTTP.WebSockets.send(socket, message)
            end
        end
    else
        handler = HTTP.streamhandler() do request
            authorized || return HTTP.Response(403)
            return request.target == "/health" ? HTTP.Response(200, run_id) : HTTP.Response(404)
        end
        handler(stream)
    end
end
file = ENV["LCM_UI_READY_FILE"]
write(file * ".pending", JSON3.write((schema_version=1, run_id, pid=getpid(), port=HTTP.port(server))))
mv(file * ".pending", file)
try
    wait(server)
finally
    close(server)
end
