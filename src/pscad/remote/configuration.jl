"""
The PSCAD version that native runs require, compared with the version that the automation
server reports.
"""
const PSCAD_VERSION = "5.1.0"

"""
$(TYPEDEF)

The station that runs PSCAD for native runs and the directory it shares with the caller. A
run reads the station's listing of its TCP sockets with the `ssh` command, and borrows a
PSCAD that listens on `port`. Otherwise it starts `pscad_executable` over SSH with
`/startup:au /port:<port>` and its working directory at `shared_root`, and quits it at the
end. It forwards PSCAD's automation port to the caller with the same command. `local_root` and `shared_root` name one directory on
the caller and on the station. Construction from field values does not perform I/O.
Construction from a TOML filename only reads that file.

$(TYPEDFIELDS)
"""
struct RemoteConfig
    "SSH destination of the station."
    host::String
    "The command that connects to the station, such as `[\"ssh\"]` or `[\"ts\", \"ssh\"]`."
    ssh::Vector{String}
    "PSCAD's automation port on the station, from 1024 to 65535."
    port::Int
    "The station's path of `Pscad.exe`."
    pscad_executable::String
    "The shared directory on the caller."
    local_root::String
    "The shared directory on the station, and PSCAD's working directory."
    shared_root::String
    "Bound of every wait of a run \\[s\\]."
    timeout::Int
end

function RemoteConfig(host::AbstractString, port::Integer, shared_root::AbstractString,
        pscad_executable::AbstractString; local_root::AbstractString,
        ssh::AbstractVector{<:AbstractString} = ["ssh"], timeout::Integer = 1800)
    isempty(strip(host)) && throw(ArgumentError("PSCAD host cannot be empty"))
    isempty(ssh) && throw(ArgumentError("PSCAD ssh command cannot be empty"))
    port in 1024:65535 || throw(ArgumentError("PSCAD port must be from 1024 to 65535"))
    isempty(strip(pscad_executable)) && throw(ArgumentError(
        "PSCAD executable cannot be empty"))
    isempty(strip(local_root)) && throw(ArgumentError("PSCAD local root cannot be empty"))
    isempty(strip(shared_root)) && throw(ArgumentError("PSCAD shared root cannot be empty"))
    timeout > 0 || throw(ArgumentError("PSCAD timeout must be positive"))
    return RemoteConfig(String(host), String.(ssh), Int(port), String(pscad_executable),
        abspath(local_root), String(shared_root), Int(timeout))
end

"""
$(TYPEDSIGNATURES)

Read a user-selected TOML file into a PSCAD configuration.

# Arguments

- `path`: TOML filename. Required fields are `host`, `port`, `local_root`, `shared_root`
  and `pscad_executable`. Optional fields are `ssh`, a non-empty array of strings, and
  `timeout`.

# Returns

- A `RemoteConfig`. A relative `local_root` is resolved against the file's
  directory. `shared_root` and `pscad_executable` are retained verbatim.
"""
function RemoteConfig(path::AbstractString)
    values = TOML.parsefile(path)
    required = ("host", "port", "local_root", "shared_root", "pscad_executable")
    unknown = setdiff(keys(values), (required..., "ssh", "timeout"))
    isempty(unknown) || throw(ArgumentError("unknown PSCAD configuration fields: $(sort!(collect(unknown)))"))
    for name in ("host", "local_root", "shared_root", "pscad_executable")
        get(values, name, nothing) isa AbstractString || throw(ArgumentError(
            "PSCAD configuration requires string field $name"))
    end
    get(values, "port", nothing) isa Int || throw(ArgumentError(
        "PSCAD configuration requires integer field port"))
    ssh = get(values, "ssh", ["ssh"])
    ssh isa Vector{String} && !isempty(ssh) || throw(ArgumentError(
        "PSCAD ssh must be a non-empty array of strings"))
    timeout = get(values, "timeout", 1800)
    timeout isa Int || throw(ArgumentError("PSCAD timeout must be an integer"))
    local_root = values["local_root"]
    isempty(strip(local_root)) && throw(ArgumentError("PSCAD local root cannot be empty"))
    return RemoteConfig(values["host"], values["port"], values["shared_root"],
        values["pscad_executable"];
        local_root = isabspath(local_root) ? local_root : joinpath(dirname(abspath(path)), local_root),
        ssh, timeout)
end
