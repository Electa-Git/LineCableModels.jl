"""
    ContainerEngine

Identify a reachable container engine independently of a Compose frontend.
Discovery is not evidence that the engine enforces an executor resource budget.
"""
struct ContainerEngine
    "Docker Engine or native Podman, never a Docker-compatible Podman alias."
    name::Symbol
    "Resolved operator-owned command."
    executable::String
    "Whether discovery identified rootless Podman; Docker policy is checked separately."
    rootless::Bool
end

function _docker_engine(; which, probe)
    executable = which("docker")
    isnothing(executable) && return nothing, "the docker command was not found", false
    version_ok, version_output = probe([executable, "version"])
    version_ok || return nothing, "docker version failed", false
    shim = occursin("podman", lowercase(version_output))
    shim && return nothing, "the docker command is a Podman compatibility shim, not Docker Engine", true
    info_ok, _ = probe([executable, "info"])
    info_ok || return nothing, "Docker Engine is not reachable", false
    return ContainerEngine(:docker, executable, false), "", false
end

function _podman_engine(; which, probe)
    executable = which("podman")
    isnothing(executable) && return nothing, "the podman command was not found"
    info_ok, _ = probe([executable, "info"])
    info_ok || return nothing, "the Podman service is not usable"
    rootless_ok, rootless_output = probe([executable, "info", "--format", "{{.Host.Security.Rootless}}"])
    rootless = rootless_ok && lowercase(strip(rootless_output)) == "true"
    return ContainerEngine(:podman, executable, rootless), ""
end

"""
    resolve_container_engine(requested="auto"; which=Sys.which, probe)

Find a reachable engine through the supplied bounded command probe. Auto prefers
real Docker Engine, recognizes its Podman compatibility shim, then tries native
Podman. No Compose dependency, image pull or resource allocation is introduced.
An explicit selection never silently changes engines. Probe receives argument
vectors and returns a success Boolean and text; it must not invoke a shell.
"""
function resolve_container_engine(requested::AbstractString="auto"; which=Sys.which, probe)
    choice = lowercase(strip(requested))
    choice in ("auto", "docker", "podman") ||
        throw(ArgumentError("container runtime must be auto, docker, or podman"))
    if choice in ("auto", "docker")
        engine, reason, shim = _docker_engine(; which, probe)
        engine === nothing || return engine
        if choice == "docker"
            hint = shim ? "; use --runtime podman on this host" : ""
            throw(ArgumentError("Docker is unavailable: $reason$hint"))
        end
    end
    if choice in ("auto", "podman")
        engine, reason = _podman_engine(; which, probe)
        engine === nothing || return engine
        choice == "podman" && throw(ArgumentError("Podman is unavailable: $reason"))
    end
    throw(ArgumentError("no usable Docker Engine or Podman command was found"))
end
