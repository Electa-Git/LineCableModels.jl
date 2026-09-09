@testset "broker roles have private reply and stream permissions" begin
    coordinator = CoordinatorIdentity()
    a, b = WorkerIdentity("worker-a"), WorkerIdentity("worker-b")
    @test RT.broker_user(a) != RT.broker_user(b)
    @test RT.broker_inbox(a) != RT.broker_inbox(b)
    @test RT.broker_inbox(a) != RT.broker_inbox(coordinator)
    @test_throws ArgumentError WorkerIdentity("worker.*")
    ap = broker_permissions(a)
    @test "lcm.science.v2.worker-a.report" in ap.publish
    @test "lcm.science.v2.worker-a.command" in ap.subscribe
    @test !("lcm.science.v2.worker-a.command" in ap.publish)
    @test "lcm.terminal.v2.worker-a.*.*.*.report" in ap.publish
    @test "lcm.terminal.v2.worker-a.*.*.*.command" in ap.subscribe
    @test !("lcm.terminal.v2.worker-a.*.*.*.command" in ap.publish)
    @test "\$JS.API.CONSUMER.MSG.NEXT.LCM_V2_JOBS_worker-a.agent" in ap.publish
    @test !any(occursin("worker-b", p) for p in vcat(ap.publish, ap.subscribe))
    @test !any(occursin("CONSUMER.CREATE", p) for p in ap.publish)
    @test !any(occursin("inbox.*", p) for p in ap.subscribe)
    cp = broker_permissions(coordinator, ("worker-a", "worker-b"))
    @test "lcm.science.v2.*.command" in cp.publish
    @test "lcm.science.v2.*.report" in cp.subscribe
    @test "lcm.terminal.v2.worker-a.*.*.*.command" in cp.publish
    @test "lcm.terminal.v2.worker-b.*.*.*.report" in cp.subscribe
    @test !any(occursin("lcm.terminal.v2.*.",p) for p in [cp.publish;cp.subscribe])
    @test "\$JS.API.CONSUMER.CREATE.LCM_V2_JOBS_worker-b.agent" in cp.publish
    @test !any(occursin("LCM_JOBS", p) for p in cp.publish)
    @test_throws ArgumentError broker_permissions(coordinator)
    @test_throws ArgumentError broker_user_config(a; password_environment="SECRET\nANYTHING")
    stanza = broker_user_config(a; password_environment="LCM_TEST_A_PASSWORD")
    @test occursin("password: \$LCM_TEST_A_PASSWORD", stanza)
    @test !occursin("test-password-value", stanza)
    @test_throws AccessDenied RT.control_subject(a,
        RT.Protocol.WorkerProbe("2.0", "worker-a", string(uuid4()), string(uuid4())))
end

@testset "broker endpoint is inert, strict and redacted" begin
    endpoint = BrokerEndpoint("tls://broker.example:4222", "/private/password";
        ca_file="/private/ca", certificate_file="/private/cert", key_file="/private/key")
    @test !occursin("/private", repr(endpoint))
    @test endpoint.url == "tls://broker.example:4222"
    @test BrokerEndpoint("nats://127.0.0.1:4222", "/private/password";
        allow_loopback_plaintext=true).url == "nats://127.0.0.1:4222"
    for url in ("nats://broker.example:4222", "nats://localhost:4222", "http://broker.example",
                "tls://user:password@broker.example", "tls://broker.example/?query=secret",
                "tls://broker.example#fragment", "tls://broker.example:70000",
                "tls://a,tls://b", "tls://broker.example/path")
        @test_throws ArgumentError BrokerEndpoint(url, "/private/password"; allow_loopback_plaintext=true)
    end
    @test_throws ArgumentError BrokerEndpoint("nats://127.0.0.1", "/private/password")
    @test_throws ArgumentError BrokerEndpoint("tls://broker.example", "/private/password";
        certificate_file="/private/cert")
    @test_throws ArgumentError BrokerEndpoint("tls://broker.example", "/private/password";
        server_name="bad\nname")
    @test_throws ArgumentError BrokerEndpoint("tls://broker.example", "")
    @test_throws ArgumentError RT.connect_broker(endpoint, CoordinatorIdentity())
end
