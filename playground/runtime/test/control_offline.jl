import HTTP, JSON3, Sockets

@testset "a stalled broker cannot hold the private gateway startup or health route" begin
    mktempdir() do dir
        password = joinpath(dir, "password")
        write(password, "fixture-not-a-real-secret"); chmod(password, 0o600)
        listener = Sockets.listen(Sockets.ip"127.0.0.1", 0)
        port = Sockets.getsockname(listener)[2]
        sockets = Sockets.TCPSocket[]
        accept_task = @async try
            while isopen(listener)
                push!(sockets, Sockets.accept(listener)) # deliberately never send INFO
            end
        catch error
            isopen(listener) && rethrow()
        end
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        profiles = ProfileRegistry()
        register!(profiles, ProfileDefinition("line-parameters", "/not-loaded", repeat("a",64);
            operations=("system.echo",)))
        config = ControlConfig(BrokerEndpoint("nats://127.0.0.1:$port", password;
            allow_loopback_plaintext=true), profiles,
            [WorkerTrust("worker-a", "credential-a", ("line-parameters",))])
        applications = ApplicationRegistry()
        supervisor = UIHostSupervisor(store, applications, joinpath(dir, "hosts"))
        service = ControlService(config, store, applications)
        start_control!(service)
        policy = LocalIdentity("http://127.0.0.1:8080", Principal("developer"; administrator=true))
        server = start_gateway(supervisor, policy; control=service)
        base = "http://127.0.0.1:$(HTTP.port(server))"
        request(path; timeout=5) = HTTP.get(base * path; headers=["Host"=>"127.0.0.1:8080"],
            proxy=nothing, retry=false, request_timeout=timeout)
        try
            @test timedwait(() -> !isempty(sockets), 10) == :ok
            # Cold HTTP specialization is measured separately from the repeated
            # responsiveness bound; no engine or UI environment is loaded.
            cold = @elapsed @test request("/health"; timeout=15).status == 200
            @info "Cold gateway health specialization (seconds)" elapsed=cold
            # Repeat after compilation, while a new handshake is stalled.
            for _ in 1:3
                elapsed = @elapsed @test request("/health").status == 200
                @test elapsed < 1
            end
            snapshot = JSON3.read(request("/runtime/api/control").body)
            @test snapshot.enabled && snapshot.broker == "unavailable"
            @test isempty(snapshot.workers) && isempty(supervisor.handles)
            @test !occursin(password, JSON3.write(snapshot))
            @test timedwait(() -> length(sockets) >= 2, 8) == :ok
            @test request("/health").status == 200
        finally
            close(server)
            close(service)
            close(listener)
            foreach(close, sockets)
            wait(accept_task)
            close(supervisor)
            close(store)
        end
        @test istaskdone(service.task)
        @test service.connector === nothing || istaskdone(service.connector)
        @test service.link.control === nothing
    end
end
