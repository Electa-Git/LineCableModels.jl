using Test, TOML
include("configure.jl")

@testset "Private demo configuration" begin
    mktempdir() do root
        private = joinpath(root,"private","file.toml")
        @test DemoSetup.save(private,Dict("value"=>1)) == private
        @test filemode(private) & 0o777 == 0o600
        @test filemode(dirname(private)) & 0o777 == 0o700
        @test_throws ErrorException DemoSetup.save(private,"replace")
        dangling = joinpath(root,"dangling")
        symlink(joinpath(root,"absent"),dangling)
        @test_throws ErrorException DemoSetup.save(dangling,"replace")
        @test !ispath(joinpath(root,"absent"))
        @test_throws ErrorException DemoSetup.operator_path("relative/path")
        @test_throws ErrorException DemoSetup.operator_path("/tmp/unsafe\npath")
        images = Dict(id=>"localhost/lcm-$id@sha256:" * repeat("a",64)
            for id in ("line-parameters","power-flow","julia-terminal"))
        profiles = DemoSetup.profiles(images)
        @test length(profiles) == 3
        @test all(p["isolation"]=="container" for p in profiles)
        @test profiles[3]["budget"]["pids"] == 64
        @test profiles[1]["budget"]["prepare_seconds"] == 600
        @test_throws ErrorException DemoSetup.profiles(merge(images,Dict("line-parameters"=>"localhost/lcm-line:latest")))
        credentials = joinpath(root,"client","coordinator")
        for file in ("coordinator.password","ca.pem","client-cert.pem","client-key.pem")
            DemoSetup.save(joinpath(credentials,file),"test-file-content")
        end
        DemoSetup.save(joinpath(credentials,"artifact-coordinator.toml"),Dict("access_key_id"=>"test","secret_access_key"=>"test-secret"))
        DemoSetup.save(joinpath(credentials,"images.toml"),images)
        config = Dict("root"=>dirname(credentials),"source"=>root,"depot"=>root,
            "julia"=>joinpath(Sys.BINDIR,Base.julia_exename()),"ssh_adapter"=>"/usr/bin/ts","host"=>"kubuntu")
        @test DemoSetup.client(config) == dirname(credentials)
        gateway = TOML.parsefile(joinpath(dirname(credentials),"gateway.toml"))
        @test gateway["gateway"]["port"] == 8081
        @test gateway["gateway"]["listen_host"] == "127.0.0.1"
        @test gateway["limits"]["disconnect_grace_seconds"] == 60
        control = TOML.parsefile(joinpath(dirname(credentials),"control.toml"))
        @test control["broker"]["server_name"] == "localhost"
        @test startswith(control["artifacts"]["endpoint"],"https://")
        @test Set(w["id"] for w in control["workers"]) == Set(DemoSetup.WORKERS)
        unit = read(joinpath(dirname(credentials),"units/lcm-demo-gateway.service"),String)
        @test occursin("KillMode=mixed",unit)
        @test !occursin("--xray",unit)
        @test !occursin("Restart=always",unit)
        tunnel = read(joinpath(dirname(credentials),"units/lcm-demo-tunnel.service"),String)
        @test occursin("127.0.0.1:24222:127.0.0.1:14222",tunnel)
        @test occursin("BatchMode=yes",tunnel)
        @test occursin("ExitOnForwardFailure=yes",tunnel)
        @test !occursin("PartOf=lcm-demo-local.target",tunnel)
        @test_throws ErrorException DemoSetup.client(config)
    end
end

@testset "Stopped infrastructure ownership checks" begin
    mktempdir() do root
        DemoSetup.save(joinpath(root,"instance"),"lcm-demo-v1\n")
        script = DemoSetup.save(joinpath(root,"engine"),raw"""#!/usr/bin/env bash
set -eu
if [[ "$1 $2" == 'container ls' ]]; then
    [[ "$LCM_DEMO_TEST_CASE" != offline ]] || exit 1
    [[ "$LCM_DEMO_TEST_CASE" == absent ]] || printf '0123456789abcdef lcm-demo-nats\n'
elif [[ "$1" == inspect ]]; then
    case "$3" in
        *Labels*) [[ "$LCM_DEMO_TEST_CASE" == foreign ]] && printf 'not-ours\n' || printf 'lcm-demo-v1\n' ;;
        *Running*) [[ "$LCM_DEMO_TEST_CASE" == running ]] && printf 'true\n' || printf 'false\n' ;;
        *Mounts*) [[ "$LCM_DEMO_TEST_CASE" == wrong-mount ]] && printf '/another/demo\n' || printf '%s/broker-data\n' "$LCM_DEMO_TEST_ROOT" ;;
        *) exit 64 ;;
    esac
elif [[ "$1 $2" == 'container rm' && "$#" == 3 && "$3" == 0123456789abcdef ]]; then
    printf 'removed\n' > "$LCM_DEMO_TEST_ROOT/removed"
else
    exit 64
fi
""")
        chmod(script,0o700)
        cleanup = `bash $(joinpath(@__DIR__,"cleanup-infrastructure")) $script $root lcm-demo-nats`
        for scenario in ("absent","offline","running","foreign","wrong-mount","owned")
            result = success(addenv(cleanup,"LCM_DEMO_TEST_CASE"=>scenario,"LCM_DEMO_TEST_ROOT"=>root))
            @test result == (scenario in ("absent","owned"))
            @test isfile(joinpath(root,"removed")) == (scenario=="owned")
        end
    end
end
