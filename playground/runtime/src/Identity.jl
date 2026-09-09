"""
    AccessDenied(status, reason)

Report a public access failure without disclosing credentials or resource owners.
"""
struct AccessDenied <: Exception
    "HTTP response status."
    status::Int
    "Fixed, non-sensitive public reason."
    reason::String
end
Base.showerror(io::IO, error::AccessDenied) = print(io, error.reason)

function checked_id(value::AbstractString)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.@-]{0,95}$", value) ||
        throw(ArgumentError("identity must contain 1–96 safe ASCII characters"))
    return String(value)
end

"""
    Principal(id; administrator=false)

Represent an authenticated caller. Browser payloads cannot assign administrator
authority; the configured identity policy owns that decision.
"""
struct Principal
    "Stable authenticated identity."
    id::String
    "Whether the configured operator policy grants administration."
    administrator::Bool
    Principal(id::AbstractString; administrator::Bool=false) =
        new(checked_id(id), administrator)
end

"""
    AbstractIdentityPolicy

Determine caller identity from an owned connection boundary.
"""
abstract type AbstractIdentityPolicy end

function checked_origin(value::AbstractString; local_only::Bool=false)
    uri = URIs.URI(value)
    isempty(uri.userinfo) && isempty(uri.query) && isempty(uri.fragment) &&
        isempty(uri.path) && !isempty(uri.host) &&
        lowercase(uri.scheme) in ("http", "https") ||
        throw(ArgumentError("public_origin must be an absolute origin without path or credentials"))
    local_only && !(uri.host in ("127.0.0.1", "[::1]", "::1")) &&
        throw(ArgumentError("local development requires a literal loopback origin"))
    uri.scheme == "http" && !local_only &&
        throw(ArgumentError("proxy identity requires HTTPS"))
    isempty(uri.port) || (tryparse(Int, uri.port) !== nothing &&
        1 <= parse(Int, uri.port) <= 65535) ||
        throw(ArgumentError("invalid origin port"))
    return String(value)
end

"""
    ProxyIdentity(origin, peers, key; administrators=())

Trust an asserted principal only when the connection peer and private proxy key
both match. The key contains at least 32 bytes; only its SHA-256 digest is kept.
Display never includes the digest or credential.
"""
struct ProxyIdentity <: AbstractIdentityPolicy
    "Exact externally configured HTTPS origin."
    origin::String
    "Literal addresses allowed to assert proxy identities."
    peers::Set{String}
    "Digest of the private proxy credential; never diagnostic metadata."
    key_digest::Vector{UInt8}
    "Identities granted explicit administration rights."
    administrators::Set{String}

    function ProxyIdentity(origin, peers, key::AbstractString; administrators=())
        ncodeunits(key) >= 32 || throw(ArgumentError("proxy key needs at least 32 bytes"))
        addresses = Set(String.(peers))
        isempty(addresses) && throw(ArgumentError("at least one proxy peer is required"))
        all(literal_ip, addresses) ||
            throw(ArgumentError("proxy peers must be literal IP addresses"))
        return new(checked_origin(origin), addresses, sha256(key),
            Set(checked_id.(administrators)))
    end
end
Base.show(io::IO, ::ProxyIdentity) = print(io, "ProxyIdentity(<redacted>)")
Base.show(io::IO, ::MIME"text/plain", policy::ProxyIdentity) = show(io, policy)

function literal_ip(value)
    try
        parse(IPAddr, value)
        return true
    catch error
        error isa ArgumentError || rethrow()
        return false
    end
end

"""
    LocalIdentity(origin, principal)

Enable explicitly configured loopback-only development identity. Client identity
headers are never used. This is not a remote deployment authentication method.
"""
struct LocalIdentity <: AbstractIdentityPolicy
    "Exact loopback browser origin."
    origin::String
    "Fixed development identity."
    principal::Principal
    LocalIdentity(origin, principal::Principal) =
        new(checked_origin(origin; local_only=true), principal)
end

function only_header(headers, name::AbstractString)
    values = [String(value) for (key, value) in headers if lowercase(key) == lowercase(name)]
    length(values) <= 1 || throw(AccessDenied(400, "Duplicate security header"))
    return isempty(values) ? "" : only(values)
end

function same_digest(a, b)
    length(a) == length(b) || return false
    difference = UInt8(0)
    for i in eachindex(a, b)
        difference |= a[i] ⊻ b[i]
    end
    return iszero(difference)
end

"""
    authenticate(policy, headers, peer) -> Principal

Authenticate request headers received from the actual connection peer. Unknown
peers, duplicate security headers, absent credentials and malformed asserted
identities raise AccessDenied.
"""
function authenticate(policy::ProxyIdentity, headers, peer::AbstractString)
    peer in policy.peers || throw(AccessDenied(401, "Untrusted proxy connection"))
    key = only_header(headers, "X-LCM-Proxy-Key")
    same_digest(sha256(key), policy.key_digest) ||
        throw(AccessDenied(401, "Proxy authentication required"))
    id = only_header(headers, "X-LCM-Principal")
    try
        return Principal(id; administrator=id in policy.administrators)
    catch error
        error isa ArgumentError || rethrow()
        throw(AccessDenied(401, "Authenticated identity required"))
    end
end

function authenticate(policy::LocalIdentity, headers, peer::AbstractString)
    peer in ("127.0.0.1", "::1") ||
        throw(AccessDenied(401, "Local development identity is loopback-only"))
    any(name -> !isempty(only_header(headers, name)),
        ("X-LCM-Principal", "X-LCM-Proxy-Key")) &&
        throw(AccessDenied(400, "Identity headers are not accepted in local development"))
    return policy.principal
end

"""
    authorize_request(policy, headers, peer; method="GET", websocket=false) -> Principal

Authenticate and enforce the exact configured Origin for mutations/WebSockets.
Mutations additionally require the non-simple X-LCM-Request header. The gateway
must not grant cross-origin CORS access. Resource ownership remains a separate
mandatory check after this request boundary.
"""
function authorize_request(policy::AbstractIdentityPolicy, headers, peer;
        method::AbstractString="GET", websocket::Bool=false)
    principal = authenticate(policy, headers, peer)
    mutation = !(method in ("GET", "HEAD", "OPTIONS"))
    if mutation || websocket
        only_header(headers, "Origin") == policy.origin ||
            throw(AccessDenied(403, "Request origin is not permitted"))
    end
    if mutation
        only_header(headers, "X-LCM-Request") == "1" ||
            throw(AccessDenied(403, "Explicit same-origin request required"))
    end
    return principal
end
