struct ContainerRuntime
    name::Symbol
    compose_command::Vector{String}
    rootless::Bool
end

function _container_probe(arguments::Vector{String})
    output = IOBuffer()
    process = run(pipeline(
        ignorestatus(Cmd(arguments));
        stdout=output,
        stderr=output,
    ))
    return success(process), String(take!(output))
end

function _compose_command(executable, runtime::Symbol; which, probe)
    integrated = [executable, "compose"]
    ok, _ = probe([integrated..., "version"])
    ok && return integrated

    standalone_name = runtime == :docker ? "docker-compose" : "podman-compose"
    standalone = which(standalone_name)
    isnothing(standalone) && return nothing
    ok, _ = probe([standalone, "version"])
    return ok ? [standalone] : nothing
end

function _docker_runtime(; which, probe)
    executable = which("docker")
    isnothing(executable) && return nothing, "the docker command was not found", false

    version_ok, version_output = probe([executable, "version"])
    version_ok || return nothing, "docker version failed", false
    podman_shim = occursin("podman", lowercase(version_output))
    podman_shim && return nothing,
        "the docker command is a Podman compatibility shim, not Docker Engine",
        true

    info_ok, _ = probe([executable, "info"])
    info_ok || return nothing, "Docker Engine is not reachable", false
    compose = _compose_command(executable, :docker; which, probe)
    isnothing(compose) && return nothing,
        "neither docker compose nor docker-compose is usable",
        false
    return ContainerRuntime(:docker, compose, false), "", false
end

function _podman_runtime(; which, probe)
    executable = which("podman")
    isnothing(executable) && return nothing, "the podman command was not found"

    info_ok, _ = probe([executable, "info"])
    info_ok || return nothing, "the Podman service is not usable"
    compose = _compose_command(executable, :podman; which, probe)
    isnothing(compose) && return nothing,
        "neither podman compose nor podman-compose is usable"

    rootless_ok, rootless_output = probe([
        executable,
        "info",
        "--format",
        "{{.Host.Security.Rootless}}",
    ])
    rootless = rootless_ok && lowercase(strip(rootless_output)) == "true"
    return ContainerRuntime(:podman, compose, rootless), ""
end

"""
    resolve_container_runtime(requested="auto"; which=Sys.which, probe=_container_probe)

Resolve a working Compose frontend. `auto` prefers a real Docker Engine on a
Docker-first host, but recognizes Docker commands implemented by a Podman shim
and selects native Podman instead.
"""
function resolve_container_runtime(
        requested::AbstractString=get(ENV, "LCM_CONTAINER_RUNTIME", "auto");
        which=Sys.which,
        probe=_container_probe,
    )
    choice = lowercase(strip(requested))
    choice in ("auto", "docker", "podman") || throw(ArgumentError(
        "container runtime must be auto, docker, or podman (received: $requested)"
    ))

    if choice in ("auto", "docker")
        runtime, reason, shim = _docker_runtime(; which, probe)
        !isnothing(runtime) && return runtime
        if choice == "docker"
            hint = shim ? "; use --runtime podman on this host" : ""
            throw(ArgumentError("Docker is unavailable: $reason$hint"))
        end
    end

    if choice in ("auto", "podman")
        runtime, reason = _podman_runtime(; which, probe)
        !isnothing(runtime) && return runtime
        choice == "podman" && throw(ArgumentError(
            "Podman is unavailable: $reason"
        ))
    end

    throw(ArgumentError(
        "no usable container runtime was found; install Docker Compose or Podman Compose, " *
        "or select one explicitly with --runtime"
    ))
end

function parse_container_options(action::String, arguments)
    runtime = get(ENV, "LCM_CONTAINER_RUNTIME", "auto")
    remote = false
    cpu_limits = false
    build = true
    volumes = false
    follow = false
    tail = nothing
    services = String[]

    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--runtime" || startswith(argument, "--runtime=")
            runtime, index = option_value(arguments, index, "--runtime")
        elseif argument == "--remote"
            remote = true
        elseif argument == "--cpu-limits"
            cpu_limits = true
        elseif argument == "--no-build"
            build = false
        elseif argument == "--volumes"
            volumes = true
        elseif argument == "--follow"
            follow = true
        elseif argument == "--tail" || startswith(argument, "--tail=")
            value, index = option_value(arguments, index, "--tail")
            tail = tryparse(Int, value)
            (isnothing(tail) || tail < 0) && throw(ArgumentError(
                "--tail must be a non-negative integer"
            ))
        elseif startswith(argument, "-")
            throw(ArgumentError("unknown container option: $argument"))
        else
            push!(services, argument)
        end
        index += 1
    end

    action == "start" || build == true || throw(ArgumentError(
        "--no-build is only valid with container start"
    ))
    action == "start" || cpu_limits == false || throw(ArgumentError(
        "--cpu-limits is only valid with container start"
    ))
    action == "stop" || volumes == false || throw(ArgumentError(
        "--volumes is only valid with container stop"
    ))
    action == "logs" || follow == false || throw(ArgumentError(
        "--follow is only valid with container logs"
    ))
    action == "logs" || isnothing(tail) || throw(ArgumentError(
        "--tail is only valid with container logs"
    ))
    action == "logs" || isempty(services) || throw(ArgumentError(
        "service names are only valid with container logs"
    ))

    return (;
        runtime,
        remote,
        cpu_limits,
        build,
        volumes,
        follow,
        tail,
        services,
    )
end

function container_deployment(remote::Bool)
    directory = remote ? joinpath(PLAYGROUND_ROOT, "deploy", "remote") :
        joinpath(PLAYGROUND_ROOT, "deploy")
    return (;
        directory,
        compose=joinpath(directory, "compose.yml"),
        cpu_limits=joinpath(directory, "compose.cpu-limits.yml"),
        env=joinpath(directory, ".env"),
        env_example=joinpath(directory, ".env.example"),
        project=remote ? "lcm-playground-remote" : "lcm-playground",
    )
end

function require_container_environment(deployment)
    isfile(deployment.env) && return deployment.env
    throw(ArgumentError(
        "container environment not found: $(deployment.env). " *
        "Create it with `cp $(deployment.env_example) $(deployment.env)`, then replace every placeholder."
    ))
end

function container_compose_command(
        action::String,
        options,
        runtime::ContainerRuntime;
        deployment=container_deployment(options.remote),
    )
    require_container_environment(deployment)

    command = [
        runtime.compose_command...,
        "--env-file",
        deployment.env,
        "--project-name",
        deployment.project,
        "--file",
        deployment.compose,
    ]
    if options.cpu_limits
        isfile(deployment.cpu_limits) || throw(ArgumentError(
            "CPU-limit override not found: $(deployment.cpu_limits)"
        ))
        append!(command, ["--file", deployment.cpu_limits])
    end

    if action == "start"
        append!(command, ["up", "--detach"])
        options.build && push!(command, "--build")
    elseif action == "stop"
        push!(command, "down")
        options.volumes && push!(command, "--volumes")
    elseif action == "status"
        push!(command, "ps")
    elseif action == "logs"
        push!(command, "logs")
        !isnothing(options.tail) && append!(command, ["--tail", string(options.tail)])
        options.follow && push!(command, "--follow")
        append!(command, options.services)
    else
        throw(ArgumentError("unknown container action: $action"))
    end
    return Cmd(Cmd(command); dir=deployment.directory), deployment
end

function print_container_runtime(io::IO, runtime::ContainerRuntime)
    println(io, "Container runtime: ", runtime.name)
    println(io, "Compose command: ", join(runtime.compose_command, " "))
    runtime.name == :podman && println(
        io,
        "Podman mode: ",
        runtime.rootless ? "rootless" : "rootful or undetermined",
    )
    return nothing
end

function run_container_action(action::String, options; io::IO=stdout)
    runtime = resolve_container_runtime(options.runtime)
    print_container_runtime(io, runtime)
    action == "resolve" && return nothing

    if options.cpu_limits && runtime.name == :podman && runtime.rootless
        @warn "CPU quotas were explicitly requested under rootless Podman; this requires a delegated CPU cgroup controller"
    end
    command, deployment = container_compose_command(action, options, runtime)
    println(io, "Deployment: ", options.remote ? "remote" : "local")
    println(io, "Compose project: ", deployment.project)
    run(command)
    return nothing
end
