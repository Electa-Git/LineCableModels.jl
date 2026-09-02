const ARTIFACT_ROUTE_PATTERN = r"^/artifacts/sha256/[0-9a-f]{64}$"

abstract type AbstractArtifactBackend end

struct FilesystemArtifactBackend <: AbstractArtifactBackend
    directory::String
end

struct S3ArtifactBackend <: AbstractArtifactBackend
    config::AWS.AbstractAWSConfig
    bucket::String
    prefix::String
end

"Minimal AWS.jl configuration adapter for a user-supplied S3-compatible endpoint."
struct S3EndpointConfig <: AWS.AbstractAWSConfig
    endpoint::URIs.URI
    region::String
    credentials::AWS.AWSCredentials
end

AWS.region(config::S3EndpointConfig) = config.region
AWS.credentials(config::S3EndpointConfig) = config.credentials

function AWS.generate_service_url(
        config::S3EndpointConfig,
        service::String,
        resource::String
    )
    service == "s3" || throw(ArgumentError(
        "S3 endpoint configuration cannot serve $service"
    ))
    return string(config.endpoint, resource)
end

"Expose local or S3-backed artifacts through an opaque same-origin route."
struct ArtifactGateway{Backend<:AbstractArtifactBackend}
    backend::Backend
end

ArtifactGateway(directory::AbstractString) = ArtifactGateway(
    FilesystemArtifactBackend(abspath(directory))
)

function normalize_s3_prefix(prefix::AbstractString)
    normalized = strip(string(prefix), '/')
    any(==(".."), split(normalized, '/')) && throw(ArgumentError(
        "S3 artifact prefix cannot contain `..` path segments"
    ))
    return normalized
end

function S3EndpointConfig(
        endpoint::AbstractString,
        access_key::AbstractString,
        secret_key::AbstractString;
        region::AbstractString="us-east-1",
        allow_insecure::Bool=false
    )
    isempty(access_key) && throw(ArgumentError("S3 access key cannot be empty"))
    isempty(secret_key) && throw(ArgumentError("S3 secret key cannot be empty"))
    uri = URIs.URI(rstrip(string(endpoint), '/'))
    uri.scheme in ("http", "https") || throw(ArgumentError(
        "S3 endpoint must use http or https"
    ))
    uri.scheme == "https" || allow_insecure || throw(ArgumentError(
        "Plain HTTP S3 endpoints require LCM_S3_ALLOW_INSECURE=1"
    ))
    isempty(uri.host) && throw(ArgumentError("S3 endpoint must include a host"))
    return S3EndpointConfig(
        uri,
        string(region),
        AWS.AWSCredentials(string(access_key), string(secret_key))
    )
end

function S3ArtifactGateway(
        endpoint::AbstractString,
        bucket::AbstractString,
        access_key::AbstractString,
        secret_key::AbstractString;
        prefix::AbstractString="linecablemodels",
        region::AbstractString="us-east-1",
        allow_insecure::Bool=false
    )
    isempty(bucket) && throw(ArgumentError("S3 bucket cannot be empty"))
    config = S3EndpointConfig(
        endpoint,
        access_key,
        secret_key;
        region,
        allow_insecure
    )
    return ArtifactGateway(S3ArtifactBackend(
        config,
        string(bucket),
        normalize_s3_prefix(prefix)
    ))
end

function default_artifact_directory()
    return abspath(get(
        ENV,
        "LCM_ARTIFACT_ROOT",
        joinpath(
            get(ENV, "XDG_DATA_HOME", joinpath(homedir(), ".local", "share")),
            "linecablemodels-worker",
            "artifacts"
        )
    ))
end

function required_artifact_environment(name::AbstractString)
    value = get(ENV, string(name), "")
    isempty(value) && throw(ArgumentError("$name is required for the S3 artifact backend"))
    return value
end

function default_artifact_gateway()
    backend = lowercase(get(ENV, "LCM_ARTIFACT_BACKEND", "filesystem"))
    if backend == "filesystem"
        return ArtifactGateway(default_artifact_directory())
    elseif backend == "s3"
        return S3ArtifactGateway(
            required_artifact_environment("LCM_S3_ENDPOINT"),
            required_artifact_environment("LCM_S3_BUCKET"),
            required_artifact_environment("AWS_ACCESS_KEY_ID"),
            required_artifact_environment("AWS_SECRET_ACCESS_KEY");
            prefix=get(ENV, "LCM_S3_PREFIX", "linecablemodels"),
            region=get(ENV, "AWS_REGION", "us-east-1"),
            allow_insecure=get(ENV, "LCM_S3_ALLOW_INSECURE", "0") == "1"
        )
    end
    throw(ArgumentError(
        "LCM_ARTIFACT_BACKEND must be `filesystem` or `s3`, got `$backend`"
    ))
end

function artifact_digest_from_target(target::AbstractString)
    path = first(split(string(target), '?'; limit=2))
    matched = match(ARTIFACT_ROUTE_PATTERN, path)
    isnothing(matched) && return nothing
    return last(split(path, '/'))
end

function artifact_object_key(
        backend::S3ArtifactBackend,
        kind::AbstractString,
        digest::AbstractString
    )
    suffix = "$kind/$digest"
    return isempty(backend.prefix) ? suffix : "$(backend.prefix)/$suffix"
end

function validate_artifact_metadata(document, digest::AbstractString)
    string(document.sha256) == digest || return nothing
    size = try
        Int(document.size)
    catch
        return nothing
    end
    0 <= size <= (Int64(1) << 40) || return nothing
    media_type = string(document.media_type)
    occursin(
        r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$",
        media_type
    ) || return nothing
    return (media_type, size)
end

function artifact_metadata(
        backend::FilesystemArtifactBackend,
        digest::AbstractString
    )
    metadata_path = joinpath(backend.directory, "$digest.metadata.json")
    isfile(metadata_path) || return nothing
    return validate_artifact_metadata(JSON3.read(read(metadata_path, String)), digest)
end

function is_missing_s3_artifact(error)
    return error isa AWS.AWSException && error.code in ("404", "NoSuchKey")
end

function with_s3_transport_retry(f; attempts::Integer=4)
    attempts > 0 || throw(ArgumentError("S3 transport attempts must be positive"))
    for attempt in 1:attempts
        try
            return f()
        catch error
            # The AWS layer owns service-level retry decisions. Retry only
            # transport failures that escaped it; artifact GETs are idempotent.
            (error isa AWS.AWSException || attempt == attempts) && rethrow()
            sleep(0.1 * 2.0^(attempt - 1))
        end
    end
    error("unreachable S3 transport retry state")
end

function artifact_metadata(backend::S3ArtifactBackend, digest::AbstractString)
    bytes = try
        with_s3_transport_retry() do
            AWSS3.s3_get(
                backend.config,
                backend.bucket,
                artifact_object_key(backend, "metadata", "$digest.json");
                raw=true,
                retry=false
            )
        end
    catch error
        is_missing_s3_artifact(error) && return nothing
        rethrow()
    end
    return validate_artifact_metadata(JSON3.read(bytes), digest)
end

function parse_artifact_range(value::AbstractString, size::Int)
    isempty(value) && return nothing
    matched = match(r"^bytes=(\d+)-(\d*)$", strip(value))
    isnothing(matched) && return :invalid
    first_byte = parse(Int, matched.captures[1])
    last_byte = isempty(matched.captures[2]) ? size - 1 : parse(Int, matched.captures[2])
    (0 <= first_byte <= last_byte < size) || return :invalid
    return (first_byte, last_byte)
end

function artifact_bytes(
        backend::FilesystemArtifactBackend,
        digest::AbstractString,
        byte_range
    )
    path = joinpath(backend.directory, digest)
    isfile(path) || return nothing
    isnothing(byte_range) && return read(path)
    first_byte, last_byte = byte_range
    return open(path, "r") do io
        seek(io, first_byte)
        read(io, last_byte - first_byte + 1)
    end
end

function artifact_bytes(
        backend::S3ArtifactBackend,
        digest::AbstractString,
        byte_range
    )
    range = isnothing(byte_range) ? nothing :
        (byte_range[1] + 1):(byte_range[2] + 1)
    return try
        with_s3_transport_retry() do
            AWSS3.s3_get(
                backend.config,
                backend.bucket,
                artifact_object_key(backend, "sha256", digest);
                raw=true,
                retry=false,
                byte_range=range
            )
        end
    catch error
        is_missing_s3_artifact(error) && return nothing
        rethrow()
    end
end

function artifact_response(gateway::ArtifactGateway, request, digest::AbstractString)
    metadata = artifact_metadata(gateway.backend, digest)
    isnothing(metadata) && return Bonito.HTTP.Response(404, "Artifact not found")
    media_type, size = metadata
    range_value = request.method == "GET" ?
        Bonito.HTTP.header(request, "Range", "") : ""
    byte_range = parse_artifact_range(range_value, size)
    byte_range === :invalid && return Bonito.HTTP.Response(
        416,
        ["Content-Range" => "bytes */$size"];
        body="Invalid artifact byte range"
    )

    response_size = isnothing(byte_range) ? size : byte_range[2] - byte_range[1] + 1
    headers = [
        "Content-Type" => media_type,
        "Content-Length" => string(response_size),
        "Accept-Ranges" => "bytes",
        "Cache-Control" => "public, max-age=31536000, immutable",
        "ETag" => "\"sha256:$digest\"",
        "X-Content-Type-Options" => "nosniff",
    ]
    if !isnothing(byte_range)
        push!(headers, "Content-Range" => "bytes $(byte_range[1])-$(byte_range[2])/$size")
    end
    request.method == "HEAD" && return Bonito.HTTP.Response(200, headers; body=UInt8[])

    body = artifact_bytes(gateway.backend, digest, byte_range)
    isnothing(body) && return Bonito.HTTP.Response(404, "Artifact not found")
    length(body) == response_size || return Bonito.HTTP.Response(
        409,
        "Artifact metadata does not match stored content"
    )
    return Bonito.HTTP.Response(
        isnothing(byte_range) ? 200 : 206,
        headers;
        body
    )
end

function Bonito.HTTPServer.apply_handler(gateway::ArtifactGateway, context)
    request = context.request
    request.method in ("GET", "HEAD") || return Bonito.HTTP.Response(
        405,
        ["Allow" => "GET, HEAD"];
        body="Method not allowed"
    )
    digest = artifact_digest_from_target(request.target)
    isnothing(digest) && return Bonito.HTTP.Response(404, "Artifact not found")
    try
        return artifact_response(gateway, request, digest)
    catch error
        diagnostic_id = string(uuid4())
        @error "Artifact backend request failed" diagnostic_id exception=(
            error,
            catch_backtrace()
        )
        return Bonito.HTTP.Response(
            503,
            "Artifact storage is temporarily unavailable (diagnostic $diagnostic_id)"
        )
    end
end

function register_artifact_route!(server, gateway::ArtifactGateway)
    if gateway.backend isa FilesystemArtifactBackend
        mkpath(gateway.backend.directory)
    end
    Bonito.route!(server, ARTIFACT_ROUTE_PATTERN => gateway)
    return gateway
end
