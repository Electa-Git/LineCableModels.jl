# Shared content-addressed storage contract for the publisher, legacy worker
# and private assigned-runtime path. Including this file performs no I/O.
"""Adapt an operator-supplied S3 endpoint to the installed AWS signing interface."""
struct S3EndpointConfig <: AWS.AbstractAWSConfig
    "S3-compatible service endpoint."
    endpoint::URIs.URI
    "S3 signing region."
    region::String
    "Server-only signing credentials, excluded from display."
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


Base.show(io::IO, ::S3EndpointConfig) = print(io, "S3EndpointConfig(server-owned credentials)")

"""Return the existing content/metadata object key below an operator-owned prefix."""
function artifact_storage_key(prefix::AbstractString, kind::AbstractString, digest::AbstractString)
    suffix = "$kind/$digest"
    return isempty(prefix) ? suffix : "$prefix/$suffix"
end

"""Describe stored bytes using the shared content-addressed metadata format."""
artifact_metadata_document(digest::AbstractString, media_type::AbstractString, size::Integer) =
    Dict("media_type"=>String(media_type), "size"=>size, "sha256"=>String(digest))
