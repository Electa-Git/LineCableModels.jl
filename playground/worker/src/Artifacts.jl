abstract type AbstractArtifactStore end

"Local content-addressed storage used by the native single-machine profile."
mutable struct ArtifactStore <: AbstractArtifactStore
    directory::String
    max_inline_bytes::Int
end

"S3-compatible content-addressed storage used by container and remote workers."
struct S3ArtifactStore <: AbstractArtifactStore
    config::AWS.AbstractAWSConfig
    bucket::String
    prefix::String
    max_inline_bytes::Int
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

const ARTIFACT_ROUTE_PREFIX = "/artifacts/sha256"

artifact_retrieval_reference(digest::AbstractString) =
    "$ARTIFACT_ROUTE_PREFIX/$digest"

artifact_metadata_path(store::ArtifactStore, digest::AbstractString) =
    joinpath(store.directory, "$digest.metadata.json")

function validate_artifact_media_type(media_type::AbstractString)
    occursin(
        r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$",
        media_type
    ) || throw(ArgumentError("invalid artifact media type"))
    return string(media_type)
end

function ArtifactStore(directory::AbstractString; max_inline_bytes::Integer=64 * 1024)
    max_inline_bytes > 0 || throw(ArgumentError("inline result limit must be positive"))
    mkpath(directory)
    return ArtifactStore(abspath(directory), Int(max_inline_bytes))
end

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

function S3ArtifactStore(
        endpoint::AbstractString,
        bucket::AbstractString,
        access_key::AbstractString,
        secret_key::AbstractString;
        prefix::AbstractString="linecablemodels",
        region::AbstractString="us-east-1",
        allow_insecure::Bool=false,
        max_inline_bytes::Integer=64 * 1024
    )
    max_inline_bytes > 0 || throw(ArgumentError("inline result limit must be positive"))
    isempty(bucket) && throw(ArgumentError("S3 bucket cannot be empty"))
    config = S3EndpointConfig(
        endpoint,
        access_key,
        secret_key;
        region,
        allow_insecure
    )
    return S3ArtifactStore(
        config,
        string(bucket),
        normalize_s3_prefix(prefix),
        Int(max_inline_bytes)
    )
end

function required_artifact_environment(name::AbstractString)
    value = get(ENV, string(name), "")
    isempty(value) && throw(ArgumentError("$name is required for the S3 artifact backend"))
    return value
end

function configured_artifact_store(data_directory::AbstractString)
    limit = parse(Int, get(ENV, "LCM_ARTIFACT_INLINE_BYTES", string(64 * 1024)))
    backend = lowercase(get(ENV, "LCM_ARTIFACT_BACKEND", "filesystem"))
    if backend == "filesystem"
        directory = get(
            ENV,
            "LCM_ARTIFACT_ROOT",
            joinpath(data_directory, "artifacts")
        )
        return ArtifactStore(directory; max_inline_bytes=limit)
    elseif backend == "s3"
        return S3ArtifactStore(
            required_artifact_environment("LCM_S3_ENDPOINT"),
            required_artifact_environment("LCM_S3_BUCKET"),
            required_artifact_environment("AWS_ACCESS_KEY_ID"),
            required_artifact_environment("AWS_SECRET_ACCESS_KEY");
            prefix=get(ENV, "LCM_S3_PREFIX", "linecablemodels"),
            region=get(ENV, "AWS_REGION", "us-east-1"),
            allow_insecure=get(ENV, "LCM_S3_ALLOW_INSECURE", "0") == "1",
            max_inline_bytes=limit
        )
    end
    throw(ArgumentError(
        "LCM_ARTIFACT_BACKEND must be `filesystem` or `s3`, got `$backend`"
    ))
end

function artifact_object_key(
        store::S3ArtifactStore,
        kind::AbstractString,
        digest::AbstractString
    )
    suffix = "$kind/$digest"
    return isempty(store.prefix) ? suffix : "$(store.prefix)/$suffix"
end

function store_artifact!(
        store::ArtifactStore,
        bytes::Vector{UInt8};
        media_type::AbstractString="application/json"
    )
    media_type = validate_artifact_media_type(media_type)
    digest = bytes2hex(SHA.sha256(bytes))
    artifact_id = "sha256:$digest"
    path = joinpath(store.directory, digest)
    if !isfile(path)
        temporary, io = mktemp(store.directory)
        try
            write(io, bytes)
            flush(io)
            close(io)
            mv(temporary, path)
        catch
            isopen(io) && close(io)
            isfile(temporary) && rm(temporary; force=true)
            rethrow()
        end
    end
    metadata_path = artifact_metadata_path(store, digest)
    if !isfile(metadata_path)
        temporary, io = mktemp(store.directory)
        try
            write(io, JSON3.write(Dict(
                "media_type" => string(media_type),
                "size" => length(bytes),
                "sha256" => digest,
            )))
            flush(io)
            close(io)
            mv(temporary, metadata_path)
        catch
            isopen(io) && close(io)
            isfile(temporary) && rm(temporary; force=true)
            rethrow()
        end
    end
    return ArtifactReference(
        artifact_id,
        string(media_type),
        length(bytes),
        digest,
        "local_filesystem",
        artifact_retrieval_reference(digest)
    )
end

function with_s3_transport_retry(f; attempts::Integer=4)
    attempts > 0 || throw(ArgumentError("S3 transport attempts must be positive"))
    for attempt in 1:attempts
        try
            return f()
        catch error
            # AWS service errors include permanent authorization and request
            # failures. The AWS stack already classifies and retries its
            # retryable service responses, so only retry transport failures
            # that escaped that layer. Content-addressed PUTs are idempotent.
            (error isa AWS.AWSException || attempt == attempts) && rethrow()
            sleep(0.1 * 2.0^(attempt - 1))
        end
    end
    error("unreachable S3 transport retry state")
end

function store_artifact!(
        store::S3ArtifactStore,
        bytes::Vector{UInt8};
        media_type::AbstractString="application/json"
    )
    media_type = validate_artifact_media_type(media_type)
    digest = bytes2hex(SHA.sha256(bytes))
    metadata = collect(codeunits(JSON3.write(Dict(
        "media_type" => string(media_type),
        "size" => length(bytes),
        "sha256" => digest,
    ))))
    with_s3_transport_retry() do
        AWSS3.s3_put(
            store.config,
            store.bucket,
            artifact_object_key(store, "sha256", digest),
            bytes,
            string(media_type);
            parse_response=false
        )
    end
    # Metadata is the commit marker. A partial upload is therefore never exposed
    # through the publisher gateway.
    with_s3_transport_retry() do
        AWSS3.s3_put(
            store.config,
            store.bucket,
            artifact_object_key(store, "metadata", "$digest.json"),
            metadata,
            "application/json";
            parse_response=false
        )
    end
    return ArtifactReference(
        "sha256:$digest",
        string(media_type),
        length(bytes),
        digest,
        "s3",
        artifact_retrieval_reference(digest)
    )
end

artifact_inline_limit(store::AbstractArtifactStore) = store.max_inline_bytes

function split_result(store::AbstractArtifactStore, result::AbstractDict)
    normalized = normalize_wire(result)
    bytes = collect(codeunits(JSON3.write(normalized)))
    if length(bytes) <= artifact_inline_limit(store)
        return (inline=normalized, artifact=nothing)
    end
    artifact = store_artifact!(store, bytes)
    return (inline=nothing, artifact)
end
