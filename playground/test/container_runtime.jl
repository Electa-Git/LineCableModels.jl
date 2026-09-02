function fake_container_host(executables, responses)
    which = name -> name in executables ? name : nothing
    probe = arguments -> get(responses, join(arguments, " "), (false, ""))
    return (; which, probe)
end

@testset "container runtime resolution" begin
    docker_host = fake_container_host(
        Set(["docker"]),
        Dict(
            "docker version" => (true, "Docker Engine 28"),
            "docker info" => (true, "Server: Docker Engine"),
            "docker compose version" => (true, "Docker Compose v2"),
        ),
    )
    runtime = LineCableModelsPlayground.resolve_container_runtime(
        "auto";
        docker_host...,
    )
    @test runtime.name == :docker
    @test runtime.compose_command == ["docker", "compose"]
    @test !runtime.rootless

    podman_host = fake_container_host(
        Set(["docker", "podman", "podman-compose"]),
        Dict(
            "docker version" => (true, "Client: Podman Engine"),
            "podman info" => (true, "host is available"),
            "podman compose version" => (true, "podman-compose 1.5"),
            "podman info --format {{.Host.Security.Rootless}}" => (true, "true\n"),
        ),
    )
    runtime = LineCableModelsPlayground.resolve_container_runtime(
        "auto";
        podman_host...,
    )
    @test runtime.name == :podman
    @test runtime.compose_command == ["podman", "compose"]
    @test runtime.rootless

    error = try
        LineCableModelsPlayground.resolve_container_runtime(
            "docker";
            podman_host...,
        )
        nothing
    catch caught
        caught
    end
    @test error isa ArgumentError
    @test occursin("Podman compatibility shim", sprint(showerror, error))
    @test occursin("--runtime podman", sprint(showerror, error))

    standalone_host = fake_container_host(
        Set(["podman", "podman-compose"]),
        Dict(
            "podman info" => (true, "host is available"),
            "podman compose version" => (false, "provider unavailable"),
            "podman-compose version" => (true, "podman-compose 1.5"),
            "podman info --format {{.Host.Security.Rootless}}" => (true, "false"),
        ),
    )
    runtime = LineCableModelsPlayground.resolve_container_runtime(
        "podman";
        standalone_host...,
    )
    @test runtime.compose_command == ["podman-compose"]
    @test !runtime.rootless

    @test_throws ArgumentError LineCableModelsPlayground.resolve_container_runtime(
        "containerd";
        docker_host...,
    )
end

@testset "container command construction" begin
    options = LineCableModelsPlayground.parse_container_options(
        "start",
        ["--runtime=podman", "--no-build", "--cpu-limits"],
    )
    @test options.runtime == "podman"
    @test !options.build
    @test options.cpu_limits
    @test !options.remote

    mktempdir() do directory
        compose = joinpath(directory, "compose.yml")
        cpu_limits = joinpath(directory, "compose.cpu-limits.yml")
        environment = joinpath(directory, ".env")
        write(compose, "services: {}\n")
        write(cpu_limits, "services: {}\n")
        write(environment, "LCM_PUBLISHER_PORT=8080\n")
        deployment = (;
            directory,
            compose,
            cpu_limits,
            env=environment,
            env_example=joinpath(directory, ".env.example"),
            project="lcm-test",
        )
        runtime = LineCableModelsPlayground.ContainerRuntime(
            :podman,
            ["podman", "compose"],
            true,
        )
        command, returned = LineCableModelsPlayground.container_compose_command(
            "start",
            options,
            runtime;
            deployment,
        )
        @test returned === deployment
        @test command.dir == directory
        @test command.exec == [
            "podman",
            "compose",
            "--env-file",
            environment,
            "--project-name",
            "lcm-test",
            "--file",
            compose,
            "--file",
            cpu_limits,
            "up",
            "--detach",
        ]
    end

    logs = LineCableModelsPlayground.parse_container_options(
        "logs",
        ["--remote", "--follow", "--tail", "50", "worker", "publisher"],
    )
    @test logs.remote
    @test logs.follow
    @test logs.tail == 50
    @test logs.services == ["worker", "publisher"]
    @test_throws ArgumentError LineCableModelsPlayground.parse_container_options(
        "status",
        ["--volumes"],
    )
    @test_throws ArgumentError LineCableModelsPlayground.parse_container_options(
        "stop",
        ["worker"],
    )
end
