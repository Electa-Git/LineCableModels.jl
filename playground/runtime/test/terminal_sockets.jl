import HTTP, JSON3

@testset "terminal socket frames use strict bounded JSON" begin
    @test RT.read_terminal_browser_frame("{\"action\":\"attach\"}")["action"]=="attach"
    for message in (UInt8[65],"[]","{\"action\":1,\"action\":2}",repeat(" ",65537),"not json")
        @test_throws AccessDenied RT.read_terminal_browser_frame(message)
    end
end

@testset "HTTP terminal queue adapter is local, bounded and closes under flood" begin
    owned=Ref{Any}(nothing);entered=Channel{Nothing}(1);finished=Channel{Nothing}(1)
    server=HTTP.listen!("127.0.0.1",0) do stream
        HTTP.WebSockets.upgrade(stream;maxframesize=65536,maxfragmentation=1) do socket
            owned[]=socket
            try
                RT.bound_terminal_socket!(socket)
                put!(entered,nothing)
                # Deliberately do not consume incoming frames. The codec reader
                # must block at four messages, not grow its default Inf queue.
                take!(finished)
            finally
                RT.abort_terminal_socket!(socket)
            end
        end
    end
    client=nothing
    try
        client=HTTP.WebSockets.open("ws://127.0.0.1:$(HTTP.port(server))/";proxy=nothing,cookies=false,request_timeout=5)
        @test timedwait(()->isready(entered),5)==:ok
        take!(entered)
        for _ in 1:32;HTTP.WebSockets.send(client,"{}");end
        # Julia counts the reader's one blocked put! as available too.
        @test timedwait(()->Base.n_avail(owned[].readchannel)==RT.TERMINAL_SOCKET_QUEUE+1,5)==:ok
        @test lock(()->length(owned[].readchannel.data),owned[].readchannel)==4
        @test client.readchannel.sz_max>4 # unrelated HTTP sockets were not patched
        @test_throws ArgumentError RT.bound_terminal_socket!(owned[])
        close_task=@async RT.abort_terminal_socket!(owned[])
        @test timedwait(()->istaskdone(close_task),7)==:ok
        @test !istaskfailed(close_task)
        @test timedwait(()->istaskdone(owned[].readtask),3)==:ok
    finally
        isready(finished) || put!(finished,nothing)
        client === nothing || close(client)
        close(server)
    end
end
