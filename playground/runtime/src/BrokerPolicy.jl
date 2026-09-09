"""
    AbstractBrokerIdentity

Select a fixed broker permission role. These identities describe server-provisioned
credentials; constructing one does not approve or authenticate a worker.
"""
abstract type AbstractBrokerIdentity end

"""Identify the single trusted runtime coordinator on the control connection."""
struct CoordinatorIdentity <: AbstractBrokerIdentity end

"""Identify one provisioned worker whose subjects are restricted by the broker."""
struct WorkerIdentity <: AbstractBrokerIdentity
    "Literal worker subject token."
    worker_id::String
    WorkerIdentity(id::AbstractString) = new(Protocol.runtime_token(id))
end

broker_user(::CoordinatorIdentity) = "lcm-coordinator"
broker_user(identity::WorkerIdentity) = "lcm-worker-" * identity.worker_id
broker_inbox(::CoordinatorIdentity) = "lcm.inbox.v2.coordinator."
broker_inbox(identity::WorkerIdentity) = "lcm.inbox.v2.worker.$(identity.worker_id)."
job_stream(id::AbstractString) = "LCM_V2_JOBS_" * Protocol.runtime_token(id)
result_stream(id::AbstractString) = "LCM_V2_RESULTS_" * Protocol.runtime_token(id)

"""
    broker_permissions(identity, worker_ids=())

Return literal least-privilege NATS grants for the version-2 transport. Each
worker has its own job/result stream and private reply prefix. Only the
coordinator creates consumers; workers cannot widen consumer filters or pull
another worker's stream. Version-1 subjects and administrative stream deletion
are not granted.
"""
function broker_permissions(identity::WorkerIdentity, worker_ids=())
    id = identity.worker_id
    jobs = job_stream(id)
    results = result_stream(id)
    return (
        publish=["lcm.report.v2.$id", "lcm.ack.v2.$id", "lcm.science.v2.$id.report", "lcm.events.v2.$id.>",
            "lcm.terminal.v2.$id.*.*.*.report",
            "lcm.results.v2.$id.>", "\$JS.API.CONSUMER.MSG.NEXT.$jobs.agent",
            "\$JS.API.CONSUMER.INFO.$jobs.agent", "\$JS.ACK.$jobs.agent.>",
            "\$JS.API.STREAM.INFO.$results", "\$JS.API.DIRECT.GET.$results", "\$JS.API.STREAM.MSG.GET.$results"],
        subscribe=["lcm.control.v2.$id.probe", "lcm.control.v2.$id.lease", "lcm.science.v2.$id.command",
            "lcm.terminal.v2.$id.*.*.*.command",broker_inbox(identity) * "*"],
    )
end

function broker_permissions(identity::CoordinatorIdentity, worker_ids=())
    ids = sort!(unique(Protocol.runtime_token.(collect(worker_ids))))
    isempty(ids) && throw(ArgumentError("coordinator permissions require provisioned worker identities"))
    length(ids) <= 4096 || throw(ArgumentError("too many provisioned broker identities"))
    publish = ["lcm.control.v2.*.probe", "lcm.control.v2.*.lease", "lcm.science.v2.*.command", "lcm.jobs.v2.>"]
    terminal_reports=String[]
    for id in ids
        push!(publish,"lcm.terminal.v2.$id.*.*.*.command")
        push!(terminal_reports,"lcm.terminal.v2.$id.*.*.*.report")
        jobs, results = job_stream(id), result_stream(id)
        for stream in (jobs, results), action in ("CREATE", "UPDATE", "INFO")
            push!(publish, "\$JS.API.STREAM.$action.$stream")
        end
        append!(publish, ["\$JS.API.CONSUMER.CREATE.$jobs.agent",
            "\$JS.API.CONSUMER.INFO.$jobs.agent", "\$JS.API.CONSUMER.DELETE.$jobs.agent",
            "\$JS.API.STREAM.MSG.GET.$results", "\$JS.API.DIRECT.GET.$results"])
    end
    return (publish=publish,
        subscribe=vcat(["lcm.report.v2.*", "lcm.ack.v2.*", "lcm.science.v2.*.report", "lcm.events.v2.>", "lcm.results.v2.>",
            broker_inbox(identity) * "*"],terminal_reports))
end

"""
    broker_user_config(identity; password_environment, worker_ids=()) -> String

Render a NATS user stanza with an environment-variable reference, never a secret
value. Operators install this stanza in their broker configuration and provision
the matching private client password file. Rendering does not modify a broker.
"""
function broker_user_config(identity::AbstractBrokerIdentity;
        password_environment::AbstractString, worker_ids=())
    occursin(r"^[A-Z][A-Z0-9_]{0,95}$", password_environment) ||
        throw(ArgumentError("invalid broker password environment name"))
    permissions = broker_permissions(identity, worker_ids)
    grants = JSON3.write((publish=(allow=permissions.publish,), subscribe=(allow=permissions.subscribe,)))
    return "{ user: \"$(broker_user(identity))\", password: \$$password_environment, permissions: $grants }"
end

"""
    BrokerEndpoint(url, password_file; ca_file=nothing, certificate_file=nothing,
        key_file=nothing, server_name=nothing, allow_loopback_plaintext=false)

Hold server-owned connection configuration without opening a connection.
Credentials are read only when connecting. Remote endpoints require verified
TLS; plaintext is an explicit literal-loopback development exception. Passwords,
query strings and fragments are forbidden in URLs.
"""
struct BrokerEndpoint
    "Credential-free NATS or TLS endpoint."
    url::String
    "Private password file, never forwarded to a browser."
    password_file::String
    "Optional operator-provided certificate authority file."
    ca_file::Union{Nothing,String}
    "Optional client certificate, paired with key_file."
    certificate_file::Union{Nothing,String}
    "Optional private client key file."
    key_file::Union{Nothing,String}
    "Optional verified TLS certificate hostname."
    server_name::Union{Nothing,String}
end

function BrokerEndpoint(url::AbstractString, password_file::AbstractString;
        ca_file=nothing, certificate_file=nothing, key_file=nothing, server_name=nothing,
        allow_loopback_plaintext::Bool=false)
    endpoint = try
        URIs.URI(url)
    catch
        throw(ArgumentError("invalid broker endpoint"))
    end
    endpoint.scheme in ("tls", "nats") && !isempty(endpoint.host) &&
        isempty(endpoint.userinfo) && isempty(endpoint.path) && isempty(endpoint.query) &&
        isempty(endpoint.fragment) && !any(isspace, url) && !occursin(',', url) ||
        throw(ArgumentError("broker endpoint must contain only scheme, host and optional port"))
    port = isempty(endpoint.port) ? 4222 : tryparse(Int, endpoint.port)
    port !== nothing && 1 <= port <= 65535 || throw(ArgumentError("invalid broker port"))
    endpoint.scheme == "tls" || (allow_loopback_plaintext && endpoint.host in ("127.0.0.1", "[::1]", "::1")) ||
        throw(ArgumentError("broker requires TLS except for explicit literal-loopback development"))
    (certificate_file === nothing) == (key_file === nothing) ||
        throw(ArgumentError("broker client certificate and key must be provided together"))
    if server_name !== nothing
        server_name isa AbstractString && occursin(r"^[A-Za-z0-9][A-Za-z0-9.:-]{0,252}$", server_name) ||
            throw(ArgumentError("invalid broker TLS server name"))
    end
    endpoint.scheme == "tls" || all(isnothing, (ca_file, certificate_file, key_file, server_name)) ||
        throw(ArgumentError("TLS files require a TLS endpoint"))
    paths = (password_file, ca_file, certificate_file, key_file)
    all(value -> value === nothing || (value isa AbstractString && !isempty(value)), paths) ||
        throw(ArgumentError("broker credential paths must be nonempty strings"))
    normalized = map(value -> value === nothing ? nothing : abspath(value), paths)
    return BrokerEndpoint(String(url), normalized..., server_name === nothing ? nothing : String(server_name))
end

Base.show(io::IO, ::BrokerEndpoint) = print(io, "BrokerEndpoint(server-owned credentials)")

function broker_file(path::AbstractString; private=false, max_bytes=1024^2)
    isfile(path) && !islink(path) && 0 < filesize(path) <= max_bytes ||
        throw(ArgumentError("broker credential must be a bounded regular file"))
    !private || iszero(filemode(path) & 0o077) ||
        throw(ArgumentError("broker secret file must have mode 0600 or stricter"))
    return path
end
