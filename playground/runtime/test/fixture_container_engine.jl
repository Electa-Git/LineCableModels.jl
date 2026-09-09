# Test-only bridge for shell harnesses. Production owns engine detection and
# command-environment filtering; this file adds no deployment launch policy.
module FixtureContainerEngine

using LineCableModelsRuntime
const RT = LineCableModelsRuntime

shell_word(value::AbstractString) = "'" * replace(value, "'" => "'\"'\"'") * "'"

function launcher(command::Vector{String}, environment::AbstractDict)
    isempty(command) && throw(ArgumentError("An inspected local command is required"))
    isabspath(first(command)) || throw(ArgumentError("The engine executable must be absolute"))
    env = Sys.which("env")
    env === nothing && throw(ArgumentError("The env executable is required"))
    words = [env, "-i"]
    append!(words, ["$key=$value" for (key, value) in sort!(collect(environment); by=first)])
    append!(words, command)
    return "#!/bin/sh\nexec " * join(shell_word.(words), " ") * " \"\$@\"\n"
end

function main(arguments)
    length(arguments) == 2 || throw(ArgumentError("Expected runtime name and new launcher path"))
    requested, destination = arguments
    requested in ("podman", "docker") || throw(ArgumentError("Expected podman or docker"))
    parent = dirname(abspath(destination))
    isdir(parent) && !islink(parent) && realpath(parent) == parent &&
        stat(parent).uid == Libc.getuid() && stat(parent).mode & 0o077 == 0 ||
        throw(ArgumentError("The launcher requires an owned private directory"))
    !ispath(destination) && !islink(destination) || throw(ArgumentError("The launcher must be new"))
    runner = CommandRunner()
    try
        host = check_container_host(runner; requested)
        :local_engine_required in host.failures && throw(ArgumentError("Test engine must be local"))
        # Transport fixtures do not claim executor quota certification. Preserve
        # the exact inspected prefix and the same filtered environment for every
        # command, including cleanup and Julia-triggered pause/unpause actions.
        write(destination, launcher(host.command, RT.container_command_environment()))
        chmod(destination, 0o700)
        println("Actual test engine: ", host.engine.name, " (pinned local target)")
    finally
        close(runner)
    end
end

end

abspath(PROGRAM_FILE) == (@__FILE__) && FixtureContainerEngine.main(ARGS)
