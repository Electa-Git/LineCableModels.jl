const PSCAD_TIMING_SCOPE = "PSCAD compile call; excludes output-readiness wait and transfer"

"""
$(TYPEDEF)

Station connection and filesystem mapping for native PSCAD execution.
`local_root` and `shared_root` name the same directory on the caller and station.
`remote_root` is scratch space on the station. Construction performs no I/O.

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
        timeout_seconds::Integer = 1800
)
    verbosity === nothing || throw(ArgumentError(
        "RemoteConfig no longer owns verbosity; set it with " *
        "options=(verbosity=(default=0, PSCAD=2),) in compute",
    ))
    isempty(strip(host)) && throw(ArgumentError("PSCAD host cannot be empty"))
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
    return RemoteConfig(
        abspath(local_root),
        String(host),
        String(shared_root),
        String(remote_root),
        String(julia_executable),
        String(python_executable),
        String(pscad_version),
        transport,
        Int(timeout_seconds)
    )
end
