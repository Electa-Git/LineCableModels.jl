"""
    ContainerHostCheck

Describe engine-level prerequisites without allocating a container. Passing this
check does not certify the limits of a future process: the physical driver must
also inspect the created container and its effective kernel limits.
"""
struct ContainerHostCheck
    "Resolved engine identity, independent of Compose."
    engine::ContainerEngine
    "Explicit local command prefix; Docker is pinned to its inspected Unix socket."
    command::Vector{String}
    "Whether the runtime reported rootless operation."
    rootless::Bool
    "Finite prerequisite failure codes, not raw engine diagnostics."
    failures::Tuple{Vararg{Symbol}}
end
Base.show(io::IO, check::ContainerHostCheck) = print(io, "ContainerHostCheck(",
    check.engine.name, ", prerequisites=", isempty(check.failures) ? "present" : "unavailable", ")")

function inspected_container_json(runner, arguments; probe=arguments -> container_probe(runner, arguments))
    ok, output = probe(arguments)
    ok || throw(ArgumentError("container host inspection failed"))
    ncodeunits(output) <= 4 * 1024^2 || throw(ArgumentError("container host inspection exceeds its byte limit"))
    parsed = try
        JSON3.read(output)
    catch
        nothing
    end
    parsed === nothing && throw(ArgumentError("container host inspection is not valid JSON"))
    return parsed
end

function podman_host_check(engine::ContainerEngine, info)
    host = get(info, :host, nothing)
    host isa JSON3.Object || throw(ArgumentError("Podman host inspection has an unsupported schema"))
    failures = Symbol[]
    get(host, :os, nothing) == "linux" || push!(failures, :linux_required)
    get(host, :serviceIsRemote, nothing) === false || push!(failures, :local_engine_required)
    get(host, :cgroupVersion, nothing) == "v2" || push!(failures, :cgroup_v2_required)
    controllers = get(host, :cgroupControllers, nothing)
    if controllers isa AbstractVector && all(c -> c isa AbstractString, controllers)
        "cpu" in controllers || push!(failures, :cpu_controller_missing)
        "memory" in controllers || push!(failures, :memory_controller_missing)
        "pids" in controllers || push!(failures, :pids_controller_missing)
    else
        push!(failures, :controllers_unverified)
    end
    security = get(host, :security, nothing)
    security isa JSON3.Object || throw(ArgumentError("Podman security inspection has an unsupported schema"))
    get(security, :seccompEnabled, nothing) === true || push!(failures, :seccomp_unavailable)
    rootless = get(security, :rootless, nothing)
    rootless isa Bool || push!(failures, :rootless_mode_unverified)
    return ContainerHostCheck(engine, [engine.executable, "--remote=false"],
        rootless === true, Tuple(failures))
end

function local_docker_endpoint(endpoint)
    return endpoint isa AbstractString && ncodeunits(endpoint) <= 1024 &&
        startswith(endpoint, "unix:///") && !occursin(r"[\x00\r\n]", endpoint) &&
        !occursin(r"/\.\.?(/|$)", endpoint) && !occursin('?', endpoint) && !occursin('#', endpoint)
end

function docker_host_check(engine::ContainerEngine, info, endpoint)
    info isa JSON3.Object || throw(ArgumentError("Docker host inspection has an unsupported schema"))
    failures = Symbol[]
    get(info, :OSType, nothing) == "linux" || push!(failures, :linux_required)
    get(info, :CgroupVersion, nothing) == "2" || push!(failures, :cgroup_v2_required)
    get(info, :MemoryLimit, nothing) === true || push!(failures, :memory_controller_missing)
    get(info, :SwapLimit, nothing) === true || push!(failures, :swap_limit_unavailable)
    (get(info, :CpuCfsPeriod, nothing) === true && get(info, :CpuCfsQuota, nothing) === true) ||
        push!(failures, :cpu_controller_missing)
    get(info, :PidsLimit, nothing) === true || push!(failures, :pids_controller_missing)
    security = get(info, :SecurityOptions, nothing)
    valid_security = security isa AbstractVector && all(s -> s isa AbstractString, security)
    valid_security && any(s -> occursin(r"^name=seccomp(?:,|$)", s), security) ||
        push!(failures, :seccomp_unavailable)
    rootless = valid_security && any(s -> s == "name=rootless", security)
    local_endpoint = local_docker_endpoint(endpoint)
    local_endpoint || push!(failures, :local_engine_required)
    command = local_endpoint ? [engine.executable, "--host", String(endpoint)] : String[]
    return ContainerHostCheck(engine, command, rootless, Tuple(failures))
end

"""
    check_container_host(runner; requested="auto", which=Sys.which) -> ContainerHostCheck

Read engine prerequisites through bounded, credential-filtered commands. Auto
uses the same Docker/Podman detection as the deployment CLI but does not require
Compose. Remote engine contexts, non-Linux engines, missing cgroup-v2 controllers,
unverified schema and missing seccomp are not accepted. A remote *agent* is still
supported: its engine must be local to that agent.

This diagnostic never pulls an image, creates a container or declares an executor
ready. Inspecting effective per-container limits remains a mandatory launch step.
"""
function check_container_host(runner::CommandRunner; requested::AbstractString="auto",
        which=Sys.which, probe=arguments -> container_probe(runner, arguments))
    engine = resolve_container_engine(requested; which, probe)
    if engine.name == :podman
        info = inspected_container_json(runner,
            [engine.executable, "--remote=false", "info", "--format", "json"]; probe)
        return podman_host_check(engine, info)
    end
    endpoint = inspected_container_json(runner,
        [engine.executable, "context", "inspect", "--format", "{{json .Endpoints.docker.Host}}"]; probe)
    # Inspect the same explicit socket that future lifecycle commands will use.
    local_endpoint = endpoint isa AbstractString && startswith(endpoint, "unix:///") &&
        ncodeunits(endpoint) <= 1024 && !occursin(r"[\x00\r\n]", endpoint)
    info_command = local_endpoint ? [engine.executable, "--host", String(endpoint), "info", "--format", "{{json .}}"] :
        [engine.executable, "info", "--format", "{{json .}}"]
    info = inspected_container_json(runner, info_command; probe)
    return docker_host_check(engine, info, endpoint)
end

"""
    recheck_container_host(runner, host) -> (current, scope)

Reinspect the already selected local engine with one bounded info request. Use
that same response to validate host prerequisites and derive its storage/daemon
identity. Do not repeat discovery or follow a subsequently changed Docker context.
No prior prerequisite response is cached; container policy and effective kernel
limits remain separate mandatory checks.
"""
function recheck_container_host(runner::CommandRunner,host::ContainerHostCheck;
        probe=arguments->container_probe(runner,arguments),scope_options...)
    engine = host.engine
    if engine.name == :podman
        host.command == [engine.executable,"--remote=false"] ||
            throw(ArgumentError("Podman command is not pinned locally"))
    elseif engine.name == :docker
        length(host.command)==3 && host.command[1:2]==[engine.executable,"--host"] &&
            local_docker_endpoint(last(host.command)) ||
            throw(ArgumentError("Docker command is not pinned locally"))
    else
        throw(ArgumentError("unsupported container engine"))
    end
    format = engine.name == :podman ? "json" : "{{json .}}"
    info = inspected_container_json(runner,[host.command;"info";"--format";format];probe)
    current = engine.name == :podman ? podman_host_check(engine,info) :
        docker_host_check(engine,info,last(host.command))
    current.command == host.command || throw(ArgumentError("local container command changed"))
    return current,container_scope_identity(current,info;scope_options...)
end
