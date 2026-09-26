const PSCAD_TIMING_SCOPE = "PSCAD compile call; excludes output-readiness wait and transfer"

"""
$(TYPEDEF)

Station connection and filesystem mapping for native PSCAD execution.
`local_root` and `shared_root` name the same directory on the caller and station.
`remote_root` is scratch space on the station. Construction from field values
performs no I/O. Construction from a TOML filename only reads that file.

$(TYPEDFIELDS)
"""
struct RemoteConfig
    local_root::String
    host::String
    shared_root::String
    remote_root::String
    julia_executable::String
    python_executable::String
    pscad_version::String
    transport::Symbol
    timeout_seconds::Int
    "Argument array for `:command` transport; exact `{host}` arguments are replaced."
    command::Vector{String}
end

function RemoteConfig(
        host::AbstractString,
        shared_root::AbstractString,
        remote_root::AbstractString,
        julia_executable::AbstractString,
        python_executable::AbstractString;
        local_root::AbstractString,
        pscad_version::AbstractString = "5.1.0",
        transport::Symbol = :ssh,
        verbosity = nothing,
        timeout_seconds::Integer = 1800,
        command::AbstractVector{<:AbstractString} = String[]
)
    verbosity === nothing || throw(ArgumentError(
        "RemoteConfig no longer owns verbosity; set it with " *
        "options=(verbosity=(default=0, PSCAD=2),) in compute",
    ))
    isempty(strip(host)) && throw(ArgumentError("PSCAD host cannot be empty"))
    isempty(strip(local_root)) && throw(ArgumentError("PSCAD local root cannot be empty"))
    isempty(strip(shared_root)) && throw(ArgumentError(
        "PSCAD shared root cannot be empty",
    ))
    isempty(strip(remote_root)) && throw(ArgumentError("PSCAD remote root cannot be empty"))
    isempty(strip(julia_executable)) && throw(ArgumentError(
        "PSCAD-host Julia executable cannot be empty",
    ))
    isempty(strip(python_executable)) && throw(ArgumentError(
        "PSCAD-host Python executable cannot be empty",
    ))
    pscad_version == "5.1.0" || throw(ArgumentError(
        "this PSCAD adapter supports version 5.1.0 only",
    ))
    timeout_seconds > 0 || throw(ArgumentError(
        "PSCAD timeout_seconds must be positive",
    ))
    if transport === :command
        isempty(command) && throw(ArgumentError("PSCAD command transport requires an argument array"))
        isempty(strip(first(command))) && throw(ArgumentError("PSCAD command executable cannot be empty"))
    else
        isempty(command) || throw(ArgumentError("PSCAD command arguments require transport=:command"))
    end
    return RemoteConfig(
        abspath(local_root),
        String(host),
        String(shared_root),
        String(remote_root),
        String(julia_executable),
        String(python_executable),
        String(pscad_version),
        transport,
        Int(timeout_seconds),
        String.(command)
    )
end

"""
$(TYPEDSIGNATURES)

Read a user-selected TOML file into a station configuration. No directory is
created and no command is executed. The backend does not discover configuration
files or read environment variables.

# Arguments

- `path`: TOML filename. Required fields are `host`, `local_root`, `shared_root`,
  `remote_root`, `julia_executable`, and `python_executable`. Optional fields are
  `pscad_version`, `transport`, `timeout_seconds`, and `command`.

# Returns

- A `RemoteConfig`. A relative `local_root` is resolved against the file's
  directory; station-side paths are retained verbatim.

# Notes

`transport` is a string naming a transport method. For `"command"`, `command`
is an argument array, with exact `{host}` arguments replaced by `host`.
Encoded PowerShell arguments are appended without shell evaluation.
"""
function RemoteConfig(path::AbstractString)
    values = TOML.parsefile(path)
    required = ("host", "local_root", "shared_root", "remote_root",
        "julia_executable", "python_executable")
    optional = ("pscad_version", "transport", "timeout_seconds", "command")
    unknown = setdiff(keys(values), (required..., optional...))
    isempty(unknown) || throw(ArgumentError("unknown PSCAD configuration fields: $(sort!(collect(unknown)))"))
    for name in required
        get(values, name, nothing) isa AbstractString || throw(ArgumentError(
            "PSCAD configuration requires string field $name"))
    end
    for name in ("pscad_version", "transport")
        !haskey(values, name) || values[name] isa AbstractString || throw(ArgumentError(
            "PSCAD configuration field $name must be a string"))
    end
    timeout = get(values, "timeout_seconds", 1800)
    timeout isa Integer && !(timeout isa Bool) || throw(ArgumentError(
        "PSCAD timeout_seconds must be an integer"))
    command = get(values, "command", String[])
    command isa AbstractVector && all(value -> value isa AbstractString, command) ||
        throw(ArgumentError("PSCAD command must be an array of strings"))
    local_root = values["local_root"]
    isempty(strip(local_root)) && throw(ArgumentError("PSCAD local root cannot be empty"))
    return RemoteConfig(values["host"], values["shared_root"], values["remote_root"],
        values["julia_executable"], values["python_executable"];
        local_root = isabspath(local_root) ? local_root : joinpath(dirname(abspath(path)), local_root),
        pscad_version = get(values, "pscad_version", "5.1.0"),
        transport = Symbol(get(values, "transport", "ssh")),
        timeout_seconds = timeout, command = String[value for value in command])
end
