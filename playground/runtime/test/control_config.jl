function control_config_data()
    Dict{String,Any}("schema_version"=>1,
        "broker"=>Dict("url"=>"nats://127.0.0.1:4222", "password_file"=>"coordinator.password",
            "allow_loopback_plaintext"=>true),
        "profiles"=>[Dict("id"=>"line-parameters", "environment"=>"native-environment",
            "fingerprint"=>repeat("a",64), "operations"=>["system.echo"])],
        "workers"=>[Dict("id"=>"worker-a", "credential_ref"=>"credential-a",
            "profiles"=>["line-parameters"], "capacity"=>1)])
end

@testset "server-owned control configuration is inert and strict" begin
    mktempdir() do dir
        secret = joinpath(dir, "coordinator.password")
        write(secret, "configuration-test-secret")
        chmod(secret, 0o600)
        path = joinpath(dir, "control.toml")
        function parse_control(data)
            open(path, "w") do io
                TOML.print(io, data)
            end
            read_control_config(path)
        end
        config = parse_control(control_config_data())
        @test only(values(config.workers)).worker_id == "worker-a"
        @test config.profiles.definitions["line-parameters"].environment == joinpath(dir, "native-environment")
        @test !ispath(joinpath(dir, "native-environment"))
        @test config.limits.total == 128
        @test !occursin("configuration-test-secret", repr(config))
        @test !occursin(dir, repr(config))
        @test !ispath(joinpath(dir, "state"))
        @test config.artifacts===nothing
        data=control_config_data()
        data["artifacts"]=Dict("backend"=>"filesystem","root"=>"private-results")
        @test parse_control(data).artifacts.root==joinpath(dir,"private-results")
        @test !ispath(joinpath(dir,"private-results"))
        for key in ("commands", "allow_eval", "trust_browser")
            data = control_config_data(); data[key] = true
            @test_throws ArgumentError parse_control(data)
        end
        for key in ("password", "user", "auth_token")
            data = control_config_data(); data["broker"][key] = "in-browser-is-forbidden"
            @test_throws ArgumentError parse_control(data)
        end
        data = control_config_data(); data["schema_version"] = true
        @test_throws ArgumentError parse_control(data)
        data = control_config_data(); data["workers"][1]["profiles"] = ["missing"]
        @test_throws ArgumentError parse_control(data)
        data = control_config_data(); push!(data["workers"], copy(first(data["workers"])))
        @test_throws ArgumentError parse_control(data)
        data = control_config_data(); data["profiles"][1]["environment"] = ""
        @test_throws ArgumentError parse_control(data)
        data = control_config_data(); data["profiles"][1]["operations"] = ["system.echo", "system.echo"]
        @test_throws ArgumentError parse_control(data)
        data = control_config_data(); data["workers"][1]["capacity"] = true
        @test_throws ArgumentError parse_control(data)
        data = control_config_data(); data["profiles"][1]["budget"] = Dict("cpus"=>0)
        @test_throws ArgumentError parse_control(data)
        chmod(secret, 0o644)
        @test_throws ArgumentError parse_control(control_config_data())
        chmod(secret, 0o600)
        parse_control(control_config_data())
        runtimepath = joinpath(dir, "runtime.toml")
        open(runtimepath, "w") do io
            TOML.print(io, Dict("schema_version"=>1,
                "gateway"=>Dict("public_origin"=>"http://127.0.0.1:8080"),
                "identity"=>Dict("mode"=>"local-development"),
                "control"=>Dict("config_file"=>"control.toml")))
        end
        @test read_config(runtimepath).control isa ControlConfig
        output = joinpath(dir, "permissions.txt")
        open(output, "w") do io
            redirect_stdout(io) do
                runtime_cli(["runtime", "permissions", "--config", runtimepath])
            end
        end
        permissions = read(output, String)
        @test occursin("LCM_V2_COORDINATOR_PASSWORD", permissions)
        @test occursin("lcm-worker-worker-a", permissions)
        @test !occursin("configuration-test-secret", permissions)
        @test !occursin(dir, permissions)
        @test !ispath(joinpath(dir, "state"))
    end
end

@testset "control events are bounded, owner-filtered and not arbitrary log text" begin
    events = ControlEvents(; capacity=2)
    alice, bob = Principal("alice"), Principal("bob")
    RT.record_event!(events, :control_connected)
    RT.record_event!(events, :control_unavailable; owner="alice")
    result = control_events(events, bob)
    @test length(result.records) == 1 && !result.gap
    @test !haskey(first(result.records), :owner)
    RT.record_event!(events, :control_rejected; worker_id="worker-a")
    @test control_events(events, alice).gap
    @test control_events(events, alice).evicted == 1
    @test length(control_events(events, alice).records) == 2
    @test length(control_events(events, bob).records) == 1
    @test isempty(control_events(events, alice; after=3).records)
    @test typeof(control_events(events, alice).records) ===
        typeof(control_events(events, alice; after=3).records) === Vector{NamedTuple}
    @test control_events(events, alice; after=20).gap
    previous = control_events(events, alice)
    replacement = ControlEvents()
    for _ in 1:5
        RT.record_event!(replacement, :control_connected)
    end
    replaced = control_events(replacement, alice; after=previous.cursor, epoch=previous.epoch)
    @test replaced.gap && length(replaced.records) == 5
    @test replaced.epoch != previous.epoch
    @test !control_events(replacement, alice; after=5, epoch=replaced.epoch).gap
    @test_throws ArgumentError control_events(events, alice; after=-1)
    @test_throws ArgumentError RT.record_event!(events, Symbol("raw private error"))
    @test_throws ArgumentError RT.record_event!(events, :control_rejected; worker_id="bad.>")
end
