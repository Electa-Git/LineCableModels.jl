import NATS, Sockets

@testset "disconnected replay policy also covers reconnecting sockets" begin
    connection = NATS.Connection(; url="nats://127.0.0.1:1", info=nothing,
        pong_count=0, pong_received_at=0.0, reconnect_count=0, connect_init_count=0,
        send_buffer_flushed=true, send_enqueue_when_disconnected=false,
        drain_timeout=0.1, drain_poll=0.01)
    for state in (NATS.CONNECTING, NATS.DISCONNECTED, NATS.DRAINED)
        NATS.status(connection, state)
        @test_throws ErrorException NATS.publish(connection, "test.input", "not-replayable")
        @test connection.send_buffer.size == 0
    end
    NATS.status(connection, NATS.CONNECTING)
    @test startswith(NATS.new_inbox(connection), "inbox.")
    connection.inbox_prefix = "lcm.inbox.v2.worker.worker-a."
    @test startswith(NATS.new_inbox(connection), connection.inbox_prefix)
    @test startswith(NATS.new_inbox(connection, "explicit."), "explicit.")
    connection.send_enqueue_when_disconnected = true
    NATS.publish(connection, "test.input", "explicitly-queued")
    @test connection.send_buffer.size > 0
    NATS.status(connection, NATS.CONNECTED)
    sub = NATS.subscribe(connection, "test.cleanup"; channel_size=1)
    @test haskey(connection.sub_data, sub.sid)
    NATS.status(connection, NATS.DRAINED)
    @test_throws ErrorException NATS.unsubscribe(connection, sub)
    @test isempty(connection.sub_data)
    @test_throws ErrorException NATS.subscribe(connection, "test.failed"; channel_size=1)
    @test isempty(connection.sub_data)
    connection.send_enqueue_when_disconnected = false
    NATS.status(connection, NATS.CONNECTING)
    @test_throws ErrorException NATS.request(connection, "test.failed-request"; timeout=0.1)
    @test isempty(connection.sub_data)
    close(connection.send_buffer)
end

@testset "reconnect does not replay previously buffered publications" begin
    for replay in (false, true)
        connection = NATS.Connection(; url="nats://127.0.0.1:1", info=nothing,
            pong_count=0, pong_received_at=0.0, reconnect_count=0, connect_init_count=0,
            send_buffer_flushed=true, send_enqueue_when_disconnected=replay,
            drain_timeout=0.1, drain_poll=0.01)
        NATS.status(connection, NATS.CONNECTED)
        NATS.publish(connection, "test.input", "uncertain-input")
        NATS.status(connection, NATS.CONNECTING)
        NATS.reopen_send_buffer(connection)
        @test occursin("uncertain-input", String(take!(connection.send_buffer))) == replay
        close(connection.send_buffer)
    end
end

@testset "broker handshake has a finite owned-socket lifetime" begin
    for timeout in (0, -1, Inf, NaN, true, 61)
        @test_throws ArgumentError NATS.connect("nats://127.0.0.1:1"; connect_timeout=timeout)
    end
    for prefix in ("", ".", "shared.>", "worker.*.", "worker..", "worker", "a\nb.")
        @test_throws ArgumentError NATS.connect("nats://127.0.0.1:1"; inbox_prefix=prefix)
    end
    for (send_info, require_tls) in ((false, false), (true, false), (true, true))
        listener = Sockets.listen(Sockets.ip"127.0.0.1", 0)
        port = Sockets.getsockname(listener)[2]
        peer = Ref{Union{Nothing,Sockets.TCPSocket}}(nothing)
        received = Ref("")
        task = @async begin
            socket = Sockets.accept(listener)
            peer[] = socket
            try
                if send_info
                    write(socket, "INFO {\"server_id\":\"fixture\",\"server_name\":\"fixture\",\"go\":\"go1.24\",\"version\":\"2.11.0\",\"host\":\"127.0.0.1\",\"port\":4222,\"proto\":1,\"headers\":true,\"max_payload\":262144,\"tls_available\":$require_tls}\r\n")
                    flush(socket)
                end
                received[] = String(read(socket)) # no response; end only when client closes
            finally
                close(socket)
            end
        end
        try
            duration = send_info ? 3.0 : 0.2
            started = time_ns()
            @test_throws Exception NATS.connect("nats://127.0.0.1:$port";
                connect_timeout=duration, retry_on_init_fail=false,
                user=nothing, pass=nothing, auth_token=nothing, jwt=nothing,
                nkey=nothing, nkey_seed=nothing, tls_required=require_tls,
                tls_ca_path=nothing, tls_cert_path=nothing, tls_key_path=nothing)
            @test duration <= (time_ns() - started) / 1e9 < 15 # includes first-call Julia compilation
            @test timedwait(() -> istaskdone(task), 2) == :ok
            @test peer[] !== nothing && !isopen(peer[])
            fetch(task)
            if require_tls
                @test !occursin("CONNECT ", received[])
                @test !isempty(received[]) && first(codeunits(received[])) == 0x16 # TLS handshake, not plaintext auth
            elseif send_info
                @test startswith(received[], "CONNECT ")
                @test occursin("PING\r\n", received[])
            end
        finally
            close(listener)
            peer[] === nothing || close(peer[])
        end
    end
end
