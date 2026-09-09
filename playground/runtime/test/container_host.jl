using Test
using LineCableModelsRuntime
const HostRT = LineCableModelsRuntime

host_json(value) = HostRT.JSON3.read(HostRT.JSON3.write(value))
function podman_info_fixture()
    Dict("host"=>Dict("os"=>"linux", "serviceIsRemote"=>false, "cgroupVersion"=>"v2",
        "cgroupControllers"=>["cpu", "memory", "pids"],
        "security"=>Dict("seccompEnabled"=>true, "rootless"=>true)))
end

function docker_info_fixture()
    Dict("OSType"=>"linux", "CgroupVersion"=>"2", "MemoryLimit"=>true,
        "SwapLimit"=>true, "CpuCfsPeriod"=>true, "CpuCfsQuota"=>true,
        "PidsLimit"=>true, "SecurityOptions"=>["name=seccomp,profile=builtin", "name=cgroupns"])
end

@testset "live host rechecks use one fresh pinned response, never rediscovery" begin
    runner = CommandRunner()
    try
        for kind in (:docker,:podman)
            info = kind==:docker ? merge(docker_info_fixture(),Dict("ID"=>"daemon-one")) :
                merge(podman_info_fixture(),Dict("store"=>Dict("graphRoot"=>"/owned/store")))
            engine = ContainerEngine(kind,string(kind),true)
            original = kind==:docker ? HostRT.docker_host_check(engine,host_json(info),"unix:///owned/docker.sock") :
                HostRT.podman_host_check(engine,host_json(info))
            history = Vector{String}[]
            probe = args->begin
                push!(history,copy(args)); (true,HostRT.JSON3.write(info))
            end
            options = (;probe,machine_id=repeat("1",32),storage_identity=path->(path,"1","42"))
            current,scope = HostRT.recheck_container_host(runner,original;options...)
            @test isempty(current.failures) && current.command==original.command
            @test length(history)==1 && history[1][1:length(original.command)]==original.command
            @test scope==container_scope(runner,original;options...)
            kind==:docker ? (info["ID"]="daemon-two") : (info["store"]["graphRoot"]="/other/store")
            @test last(HostRT.recheck_container_host(runner,original;options...))!=scope
            kind==:docker ? (info["PidsLimit"]=false) : (info["host"]["cgroupControllers"]=["cpu","memory"])
            @test :pids_controller_missing in first(HostRT.recheck_container_host(runner,original;options...)).failures
            kind==:docker ? (info["OSType"]="windows") : (info["host"]["serviceIsRemote"]=true)
            @test_throws ArgumentError HostRT.recheck_container_host(runner,original;options...)
            @test_throws ArgumentError HostRT.recheck_container_host(runner,original;
                probe=_->(true,"bad JSON"))
            @test_throws ArgumentError HostRT.recheck_container_host(runner,original;
                probe=_->(false,"private failure"))
        end
        for (kind,command) in ((:podman,["podman","--remote=true"]),
                (:docker,["docker","--host","tcp://foreign:2375"]))
            host = ContainerHostCheck(ContainerEngine(kind,string(kind),true),command,true,())
            @test_throws ArgumentError HostRT.recheck_container_host(runner,host;probe=_->error("must not invoke"))
        end
        @test isempty(runner.active)
    finally
        close(runner)
    end
end
@testset "container prerequisites are explicit and never imply physical readiness" begin
    podman = ContainerEngine(:podman, "/usr/bin/podman", true)
    info = podman_info_fixture()
    check = HostRT.podman_host_check(podman, host_json(info))
    @test isempty(check.failures)
    @test check.rootless
    @test check.command == ["/usr/bin/podman", "--remote=false"]
    info["host"]["cgroupControllers"] = ["memory", "pids"]
    @test HostRT.podman_host_check(podman, host_json(info)).failures == (:cpu_controller_missing,)
    for (field, value, reason) in (
            ("serviceIsRemote", true, :local_engine_required),
            ("serviceIsRemote", "false", :local_engine_required),
            ("cgroupVersion", "v1", :cgroup_v2_required),
            ("os", "windows", :linux_required),
            ("cgroupControllers", nothing, :controllers_unverified),
            ("cgroupControllers", ["cpu"], :memory_controller_missing),
            ("cgroupControllers", ["cpu", "memory"], :pids_controller_missing))
        info = podman_info_fixture()
        info["host"][field] = value
        @test reason in HostRT.podman_host_check(podman, host_json(info)).failures
    end
    info = podman_info_fixture()
    info["host"]["security"]["seccompEnabled"] = false
    @test :seccomp_unavailable in HostRT.podman_host_check(podman, host_json(info)).failures
    @test_throws ArgumentError HostRT.podman_host_check(podman, host_json(Dict()))
    docker = ContainerEngine(:docker, "/usr/bin/docker", false)
    info = docker_info_fixture()
    check = HostRT.docker_host_check(docker, host_json(info), "unix:///run/docker.sock")
    @test isempty(check.failures)
    @test check.command == ["/usr/bin/docker", "--host", "unix:///run/docker.sock"]
    for endpoint in ("tcp://remote:2375", "ssh://remote", "unix:///tmp/../socket",
            "unix://relative", "unix:///tmp/socket?secret", "unix:///tmp/socket\n--bad", nothing)
        check = HostRT.docker_host_check(docker, host_json(info), endpoint)
        @test :local_engine_required in check.failures
        @test isempty(check.command)
    end
    for (field, value, reason) in (
            ("MemoryLimit", false, :memory_controller_missing),
            ("SwapLimit", false, :swap_limit_unavailable),
            ("CpuCfsQuota", false, :cpu_controller_missing),
            ("PidsLimit", 1, :pids_controller_missing),
            ("CgroupVersion", "1", :cgroup_v2_required),
            ("OSType", "windows", :linux_required),
            ("SecurityOptions", ["name=seccomp-disabled"], :seccomp_unavailable))
        info = docker_info_fixture()
        info[field] = value
        @test reason in HostRT.docker_host_check(docker, host_json(info), "unix:///run/docker.sock").failures
    end
end

@testset "shared engine selection needs no Compose or resource creation" begin
    history = Vector{String}[]
    responses = Dict(
        "docker version" => (true, "Docker Engine"),
        "docker info" => (true, "Docker Engine"),
        "docker context inspect --format {{json .Endpoints.docker.Host}}" => (true, "\"unix:///run/docker.sock\""),
        "docker --host unix:///run/docker.sock info --format {{json .}}" => (true, HostRT.JSON3.write(docker_info_fixture())),
        "podman info" => (true, "Podman"),
        "podman info --format {{.Host.Security.Rootless}}" => (true, "true\n"),
        "podman --remote=false info --format json" => (true, HostRT.JSON3.write(podman_info_fixture())))
    which = name -> name in ("docker", "podman") ? name : nothing
    probe = args -> begin
        push!(history, copy(args))
        get(responses, join(args, " "), (false, ""))
    end
    runner = CommandRunner()
    try
        check = check_container_host(runner; which, probe)
        @test check.engine.name == :docker
        @test isempty(check.failures)
        responses["docker version"] = (true, "Client: Podman Engine")
        check = check_container_host(runner; which, probe)
        @test check.engine.name == :podman
        @test isempty(check.failures)
        @test_throws ArgumentError check_container_host(runner; requested="docker", which, probe)
        check = check_container_host(runner; requested="podman", which, probe)
        @test check.engine.name == :podman
        responses["podman --remote=false info --format json"] = (true, "private invalid json")
        @test_throws ArgumentError check_container_host(runner; requested="podman", which, probe)
        @test isempty(runner.active)
        @test !any(args -> any(in(("create", "run", "start", "pull", "compose")), args), history)
        @test_throws ArgumentError check_container_host(runner; requested="containerd", which, probe)
    finally
        close(runner)
    end
end
