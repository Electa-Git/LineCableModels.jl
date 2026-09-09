import HTTP, JSON3

@testset "gateway identity, origin and owned HTTP routes" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        registry = ApplicationRegistry()
        register!(registry, LocalApplication(
            ApplicationDefinition("mock", "Mock", :workbench, "/mock"),
            dirname(@__DIR__), joinpath(@__DIR__, "ui_child.jl")))
        supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts");
            limits=RunLimits(max_runs=2, max_runs_per_owner=1, startup_seconds=30, shutdown_seconds=0.5))
        key = repeat("proxy-private-", 4)
        policy = ProxyIdentity("https://lcm.test", ["127.0.0.1"], key)
        server = start_gateway(supervisor, policy)
        base = "http://127.0.0.1:$(HTTP.port(server))"
        alice = ["X-LCM-Proxy-Key"=>key, "X-LCM-Principal"=>"alice"]
        bob = ["X-LCM-Proxy-Key"=>key, "X-LCM-Principal"=>"bob"]
        mutation(headers) = [headers; "Origin"=>"https://lcm.test"; "X-LCM-Request"=>"1"; "Content-Type"=>"application/json"]
        request(method, path, headers=Pair{String,String}[], body=nothing) =
            HTTP.request(method, base * path, headers, body; proxy=nothing,
                status_exception=false, retry=false, request_timeout=10)
        try
            @test request("GET", "/health").status == 200
            catalogue = request("GET", "/runtime/api/applications")
            @test catalogue.status == 200
            @test only(JSON3.read(catalogue.body)).id == "mock"
            @test isempty(supervisor.handles)
            @test request("GET", "/runtime/assets/brand.css").status == 200
            @test request("GET", "/runtime/assets/../Project.toml").status != 200
            @test request("GET", "/runtime/api/runs").status == 401
            @test request("GET", "/runtime/api/runs", ["X-LCM-Principal"=>"alice"]).status == 401
            @test request("GET", "/runtime/api/runs", alice).status == 200
            payload = JSON3.write((application="mock", request_id=string(uuid4())))
            @test request("POST", "/runtime/api/runs", alice, payload).status == 403
            wrong_origin = [alice; "Origin"=>"https://evil.test"; "X-LCM-Request"=>"1"]
            @test request("POST", "/runtime/api/runs", wrong_origin, payload).status == 403
            @test request("POST", "/runtime/api/runs", mutation(alice), "{}").status == 400
            launched = request("POST", "/runtime/api/runs", mutation(alice), payload)
            @test launched.status == 202
            first = UUID(JSON3.read(launched.body).id)
            repeated = request("POST", "/runtime/api/runs", mutation(alice), payload)
            @test UUID(JSON3.read(repeated.body).id) == first
            next_payload = JSON3.write((application="mock", request_id=string(uuid4())))
            @test request("POST", "/runtime/api/runs", mutation(alice), next_payload).status == 429
            launched_bob = request("POST", "/runtime/api/runs", mutation(bob), next_payload)
            @test launched_bob.status == 202
            second = UUID(JSON3.read(launched_bob.body).id)
            @test timedwait(() -> get_run(store, Principal("alice"), first).state == :running, 40) == :ok
            @test timedwait(() -> get_run(store, Principal("bob"), second).state == :running, 40) == :ok
            @test request("GET", "/runtime/api/runs/$first", bob).status == 404
            @test request("GET", "/runtime/runs/$first", bob).status == 404
            surface = request("GET", "/runtime/runs/$first", alice)
            @test surface.status == 200
            @test occursin("/runtime/assets/published-text.css", String(surface.body))
            @test request("GET", "/runtime/api/runs/$(uuid4())", bob).status == 404
            @test request("GET", "/applications/runs/$first/health", bob).status == 404
            @test request("DELETE", "/runtime/api/runs/$first", mutation(bob)).status == 404
            forwarded = request("GET", "/applications/runs/$first/health", alice)
            @test forwarded.status == 200
            @test String(forwarded.body) == string(first)
            @test !HTTP.hasheader(forwarded, "X-LCM-Host-Key")
            @test !occursin(key, String(request("GET", "/runtime/api/runs", alice).body))
            handle = supervisor.handles[first]
            direct = HTTP.get("http://127.0.0.1:$(handle.port)/health";
                status_exception=false, proxy=nothing, request_timeout=5)
            @test direct.status == 403
            socket_url = replace(base, "http:" => "ws:") * "/applications/runs/$first/echo"
            socket_headers = [alice; "Origin"=>"https://lcm.test"]
            HTTP.WebSockets.open(socket_url; headers=socket_headers,
                    proxy=nothing, request_timeout=3, read_idle_timeout=5) do socket
                HTTP.WebSockets.send(socket, "state remains in this run")
                @test HTTP.WebSockets.receive(socket) == "state remains in this run"
                payload_bytes = UInt8[0x00, 0x01, 0x02, 0xff]
                HTTP.WebSockets.send(socket, payload_bytes)
                @test HTTP.WebSockets.receive(socket) == payload_bytes
            end
            @test_throws Exception HTTP.WebSockets.open(socket_url;
                headers=alice, proxy=nothing, request_timeout=3)
            @test_throws Exception HTTP.WebSockets.open(socket_url;
                headers=[bob; "Origin"=>"https://lcm.test"], proxy=nothing, request_timeout=3)
            @test request("DELETE", "/runtime/api/runs/$first", mutation(alice)).status == 200
            @test request("GET", "/applications/runs/$first/health", alice).status == 503
            unavailable = request("GET", "/applications/runs/$first/mock", [alice; "Accept"=>"text/html"])
            @test unavailable.status == 503
            @test occursin("Start a clean run", String(unavailable.body))
            @test request("GET", "/applications/runs/$second/health", bob).status == 200
            @test request("DELETE", "/runtime/api/runs/$second", mutation(bob)).status == 200
        finally
            close(server)
            close(supervisor)
            close(store)
        end
    end
end

@testset "local development rejects rebinding authorities" begin
    mktempdir() do dir
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        supervisor = UIHostSupervisor(store, ApplicationRegistry(), joinpath(dir, "hosts"))
        server = start_gateway(supervisor, LocalIdentity("http://127.0.0.1:8080", Principal("developer")))
        url = "http://127.0.0.1:$(HTTP.port(server))/runtime/api/runs"
        try
            @test HTTP.get(url; headers=["Host"=>"127.0.0.1:8080"], proxy=nothing).status == 200
            @test HTTP.get(url; headers=["Host"=>"rebound.attacker.test"], proxy=nothing,
                status_exception=false).status == 403
        finally
            close(server); close(supervisor); close(store)
        end
    end
end
