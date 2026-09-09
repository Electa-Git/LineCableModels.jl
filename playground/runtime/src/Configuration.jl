"""
    RunLimits(; max_runs=8, max_runs_per_owner=2, startup_seconds=60,
              shutdown_seconds=5, disconnect_grace_seconds=60)

Bound application admission and process lifetimes. Durations are in seconds.
"""
struct RunLimits
    "Maximum simultaneously reserved or running UI hosts."
    max_runs::Int
    "Maximum simultaneously reserved or running UI hosts per owner."
    max_runs_per_owner::Int
    "Maximum UI-host startup duration in seconds."
    startup_seconds::Float64
    "Graceful shutdown duration before forced termination, in seconds."
    shutdown_seconds::Float64
    "Maximum disconnected-run retention duration in seconds."
    disconnect_grace_seconds::Float64

    function RunLimits(; max_runs=8, max_runs_per_owner=2, startup_seconds=60,
            shutdown_seconds=5, disconnect_grace_seconds=60)
        max_runs isa Integer && !(max_runs isa Bool) && 1 <= max_runs <= 256 ||
            throw(ArgumentError("max_runs must be an integer in 1:256"))
        max_runs_per_owner isa Integer && !(max_runs_per_owner isa Bool) &&
            1 <= max_runs_per_owner <= max_runs ||
            throw(ArgumentError("max_runs_per_owner must be in 1:max_runs"))
        durations = (startup_seconds, shutdown_seconds, disconnect_grace_seconds)
        all(v -> v isa Real && !(v isa Bool) && isfinite(v) && 0 < v <= 3600, durations) ||
            throw(ArgumentError("run durations must be finite and in (0, 3600] seconds"))
        return new(max_runs, max_runs_per_owner, Float64.(durations)...)
    end
end

"""
    RuntimeConfig

Hold validated, server-owned runtime configuration. Creation does not start
workers, connect to a broker, open a database or create directories.
"""
struct RuntimeConfig{P<:AbstractIdentityPolicy,C}
    "Explicit opt-in for runtime routes."
    enabled::Bool
    "Loopback bind address for the private gateway."
    listen_host::String
    "Configured listener port; zero requests a test-only ephemeral port."
    port::Int
    "Caller identity and origin policy."
    identity::P
    "Absolute SQLite database path."
    database::String
    "Absolute directory reserved for owned runtime resources."
    scratch_root::String
    "Optional already-built static publication directory."
    site_directory::Union{Nothing,String}
    "Admission and lifecycle limits."
    limits::RunLimits
    "Optional validated worker-control configuration; nothing leaves it disabled."
    control::C
end

function strict_keys(table::AbstractDict, allowed, label)
    isempty(setdiff(Set(keys(table)), Set(allowed))) ||
        throw(ArgumentError("unknown configuration key in $label"))
    return table
end

function config_table(data, name)
    value = get(data, name, Dict{String,Any}())
    value isa AbstractDict || throw(ArgumentError("$name must be a TOML table"))
    return value
end

function config_path(base, value)
    value isa AbstractString && !isempty(value) ||
        throw(ArgumentError("storage and credential paths must be nonempty strings"))
    resolved = normpath(joinpath(base, value))
    return resolved == "/" ? resolved : String(rstrip(resolved, '/'))
end

"""
    read_config(path) -> RuntimeConfig

Read strict schema-version 1 TOML without starting any services. Relative paths
resolve against the file directory. Only implemented tables are accepted.

# Errors

Raise ArgumentError for unknown keys, unsupported versions, unsafe listeners,
invalid resource limits or an unreadable/insecure proxy-key file.
"""
function read_config(path::AbstractString)
    data = TOML.parsefile(path)
    strict_keys(data, ("schema_version", "enabled", "gateway", "identity", "storage", "limits", "publisher", "control"), "root")
    get(data, "schema_version", nothing) === 1 ||
        throw(ArgumentError("unsupported runtime configuration version"))
    enabled = get(data, "enabled", false)
    enabled isa Bool || throw(ArgumentError("enabled must be boolean"))
    base = dirname(abspath(path))
    gateway = strict_keys(config_table(data, "gateway"),
        ("listen_host", "port", "public_origin"), "gateway")
    host = get(gateway, "listen_host", "127.0.0.1")
    host in ("127.0.0.1", "::1") ||
        throw(ArgumentError("v1 private gateway must bind a literal loopback address"))
    port = get(gateway, "port", 8080)
    port isa Integer && !(port isa Bool) && 0 <= port <= 65535 ||
        throw(ArgumentError("invalid gateway port"))
    origin = get(gateway, "public_origin", "")
    origin isa AbstractString || throw(ArgumentError("public_origin must be a string"))
    identity = config_table(data, "identity")
    mode = get(identity, "mode", "proxy")
    policy = if mode == "local-development"
        strict_keys(identity, ("mode", "principal", "administrator"), "identity")
        admin = get(identity, "administrator", false)
        admin isa Bool || throw(ArgumentError("administrator must be boolean"))
        LocalIdentity(origin, Principal(get(identity, "principal", "developer"); administrator=admin))
    elseif mode == "proxy"
        strict_keys(identity, ("mode", "proxy_peers", "proxy_key_file", "administrators"), "identity")
        keyfile = config_path(base, get(identity, "proxy_key_file", ""))
        isfile(keyfile) && !islink(keyfile) && filesize(keyfile) <= 4096 &&
            iszero(filemode(keyfile) & 0o077) ||
            throw(ArgumentError("proxy key must be a private regular file (mode 0600 or stricter)"))
        ProxyIdentity(origin, get(identity, "proxy_peers", ["127.0.0.1"]),
            strip(read(keyfile, String)); administrators=get(identity, "administrators", String[]))
    else
        throw(ArgumentError("identity mode must be proxy or local-development"))
    end
    storage = strict_keys(config_table(data, "storage"), ("database", "scratch_root"), "storage")
    database = config_path(base, get(storage, "database", "state/runtime.sqlite"))
    scratch = config_path(base, get(storage, "scratch_root", "state/runs"))
    scratch in ("/", base, dirname(base), homedir()) &&
        throw(ArgumentError("scratch_root must be a dedicated runtime directory"))
    limits = strict_keys(config_table(data, "limits"),
        ("max_runs", "max_runs_per_owner", "startup_seconds", "shutdown_seconds",
         "disconnect_grace_seconds"), "limits")
    parsed_limits = RunLimits(; (Symbol(k) => v for (k,v) in limits)...)
    publisher = strict_keys(config_table(data, "publisher"), ("site_directory",), "publisher")
    site = haskey(publisher, "site_directory") ? config_path(base, publisher["site_directory"]) : nothing
    control = strict_keys(config_table(data, "control"), ("config_file",), "control")
    parsed_control = isempty(control) ? nothing : read_control_config(config_path(base, control["config_file"]))
    return RuntimeConfig(enabled, String(host), port, policy, database, scratch, site, parsed_limits, parsed_control)
end
