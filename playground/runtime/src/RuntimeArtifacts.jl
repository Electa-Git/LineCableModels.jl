include(joinpath(@__DIR__, "..", "..", "common", "artifact_contract.jl"))

"""Maximum private result size accepted by storage and verified downloads, in bytes."""
const ASSIGNED_ARTIFACT_BYTES = 4 * 1024 * 1024
"""Result-value size above which configured private storage replaces inline delivery, in bytes."""
const ASSIGNED_INLINE_BYTES = 64 * 1024

"""Describe operator-owned storage; construction performs no storage request."""
abstract type AbstractRuntimeArtifacts end

"""Use a dedicated private filesystem root for trusted single-machine execution."""
struct LocalRuntimeArtifacts <: AbstractRuntimeArtifacts
    "Dedicated private root, separate from public artifacts and user uploads."
    root::String
    function LocalRuntimeArtifacts(root::AbstractString)
        path=abspath(root)
        path in ("/",homedir(),pwd(),dirname(pwd()),tempdir()) &&
            throw(ArgumentError("Private artifacts require a dedicated directory"))
        new(path)
    end
end

"""
    S3RuntimeArtifacts(endpoint, bucket, prefix, credentials_file;
        region="us-east-1", ca_file=nothing, allow_loopback_plaintext=false)

Locate private job-scoped content in the existing S3 content/metadata format.
Credentials remain in a bounded private TOML file and are read only for storage
requests. Remote storage requires verified TLS; plaintext is restricted to an
explicit literal-loopback test/development endpoint. No URL is browser supplied.
"""
struct S3RuntimeArtifacts <: AbstractRuntimeArtifacts
    "Credential-free approved object-storage endpoint."
    endpoint::String
    "Operator-provisioned private bucket."
    bucket::String
    "Private runtime prefix; worker/run/job scopes are appended by the owner."
    prefix::String
    "Private TOML containing access_key_id and secret_access_key only."
    credentials_file::String
    "Approved S3 signing region."
    region::String
    "Optional private-network CA bundle; hostname verification remains enabled."
    ca_file::Union{Nothing,String}
    function S3RuntimeArtifacts(endpoint::AbstractString,bucket::AbstractString,prefix::AbstractString,
            credentials_file::AbstractString;region::AbstractString="us-east-1",ca_file=nothing,
            allow_loopback_plaintext::Bool=false)
        uri=URIs.URI(endpoint)
        uri.scheme in ("https","http") && !isempty(uri.host) && isempty(uri.userinfo) &&
            isempty(uri.path) && isempty(uri.query) && isempty(uri.fragment) && !any(isspace,endpoint) ||
            throw(ArgumentError("Artifact endpoint requires only scheme, host and optional port"))
        port=isempty(uri.port) ? 443 : tryparse(Int,uri.port)
        port!==nothing && 1<=port<=65535 || throw(ArgumentError("Invalid artifact endpoint port"))
        uri.scheme=="https" || (allow_loopback_plaintext && uri.host in ("127.0.0.1","::1","[::1]")) ||
            throw(ArgumentError("Remote artifact storage requires verified HTTPS"))
        occursin(r"^[a-z0-9][a-z0-9.-]{1,61}[a-z0-9]$",bucket) || throw(ArgumentError("Invalid artifact bucket"))
        occursin(r"^[a-zA-Z0-9_-]+(?:/[a-zA-Z0-9_-]+)*$",prefix) && ncodeunits(prefix)<=128 ||
            throw(ArgumentError("Private artifact prefix must be explicit and bounded"))
        occursin(r"^[a-z0-9-]{1,64}$",region) || throw(ArgumentError("Invalid artifact region"))
        isempty(credentials_file) && throw(ArgumentError("Artifact credentials require a private file"))
        ca_file===nothing || (uri.scheme=="https" && ca_file isa AbstractString && !isempty(ca_file)) ||
            throw(ArgumentError("Artifact CA requires HTTPS and a nonempty path"))
        new(String(endpoint),String(bucket),String(prefix),abspath(credentials_file),String(region),
            ca_file===nothing ? nothing : abspath(ca_file))
    end
end
Base.show(io::IO,::AbstractRuntimeArtifacts)=print(io,"RuntimeArtifacts(server-owned private storage)")

"""Report finite storage failure without retaining URLs, credentials or server bodies."""
struct ArtifactUnavailable <: Exception end
Base.showerror(io::IO,::ArtifactUnavailable)=print(io,"Private result storage is unavailable")

function configured_runtime_artifacts(data,base)
    table=config_table(data,"artifacts")
    isempty(table) && return nothing
    backend=get(table,"backend",nothing)
    if backend=="filesystem"
        strict_keys(table,("backend","root"),"artifacts")
        root=config_path(base,get(table,"root",""))
        root in (base,dirname(base)) && throw(ArgumentError("Private artifacts require a dedicated directory"))
        return LocalRuntimeArtifacts(root)
    elseif backend=="s3"
        strict_keys(table,("backend","endpoint","bucket","prefix","credentials_file","region","ca_file",
            "allow_loopback_plaintext"),"artifacts")
        ca=haskey(table,"ca_file") ? config_path(base,table["ca_file"]) : nothing
        result=S3RuntimeArtifacts(get(table,"endpoint",""),get(table,"bucket",""),get(table,"prefix",""),
            config_path(base,get(table,"credentials_file",""));region=get(table,"region","us-east-1"),ca_file=ca,
            allow_loopback_plaintext=get(table,"allow_loopback_plaintext",false))
        broker_file(result.credentials_file;private=true,max_bytes=8192)
        ca===nothing || broker_file(ca)
        return result
    end
    throw(ArgumentError("Artifact backend must be filesystem or s3"))
end

function artifact_scope(fence::Protocol.AssignmentFence,job_id::AbstractString)
    Protocol.validate(fence); Protocol.runtime_uuid(job_id)
    join(("workers",fence.worker_id,"runs",fence.run_id,fence.worker_boot,fence.lease_id,string(fence.generation),job_id),'/')
end

function checked_artifact_reference(reference::Protocol.ArtifactReference)
    Protocol.validate(reference)
    reference.media_type=="application/json" && 0<reference.size<=ASSIGNED_ARTIFACT_BYTES ||
        throw(AccessDenied(409,"Assigned artifact violates its result contract"))
    reference
end

function private_artifact_directory(location::LocalRuntimeArtifacts,scope::String;create=false)
    path=joinpath(location.root,split(scope,'/')...)
    if !ispath(location.root)
        create || return nothing
        mkpath(location.root;mode=0o700)
    end
    realpath(location.root)==location.root || throw(ArtifactUnavailable())
    receipt_private_stat(lstat(location.root);directory=true)
    for segment in split(scope,'/')
        # Validate existing ancestors before creating any nested directory.
        location=LocalRuntimeArtifacts(joinpath(location.root,segment))
        if !ispath(location.root)
            islink(location.root) && throw(ArtifactUnavailable())
            create || return nothing
            mkdir(location.root;mode=0o700)
        end
        receipt_private_stat(lstat(location.root);directory=true)
    end
    path
end

function local_artifact_read(path::String,limit::Int)
    ispath(path) || (islink(path) && throw(ArtifactUnavailable()); return nothing)
    receipt_private_stat(lstat(path))
    file=receipt_open(path)
    try
        before=stat(file); receipt_private_stat(before)
        before.size<=limit || throw(ArtifactUnavailable())
        bytes=read(file,limit+1)
        length(bytes)==before.size && eof(file) || throw(ArtifactUnavailable())
        after=stat(file)
        (before.size,before.mtime,before.ctime)==(after.size,after.mtime,after.ctime) || throw(ArtifactUnavailable())
        bytes
    finally
        close(file)
    end
end

function local_artifact_write(path::String,bytes::Vector{UInt8})
    temporary=joinpath(dirname(path),".pending-"*string(uuid4()))
    file=receipt_open(temporary;create=true,exclusive=true,writable=true)
    try
        write(file,bytes)==length(bytes) || throw(ArtifactUnavailable())
        receipt_sync(file);close(file)
        (ispath(path)||islink(path)) && receipt_private_stat(lstat(path))
        Base.Filesystem.rename(temporary,path)
        directory=receipt_open(dirname(path);directory=true)
        try receipt_sync(directory) finally close(directory) end
    finally
        isopen(file) && close(file)
        if isfile(temporary) && !islink(temporary)
            receipt_private_stat(lstat(temporary));rm(temporary)
        end
    end
    nothing
end

function artifact_object_path(location::LocalRuntimeArtifacts,scope,kind,digest;create=false)
    parent=private_artifact_directory(location,scope;create)
    parent===nothing ? nothing : joinpath(parent,kind=="metadata" ? "$digest.metadata.json" : digest)
end

function artifact_object_request(location::LocalRuntimeArtifacts,scope,kind,digest,method,bytes=UInt8[];limit=ASSIGNED_ARTIFACT_BYTES)
    path=artifact_object_path(location,scope,kind,digest;create=method=="PUT")
    path===nothing && return nothing
    if method=="GET"
        return local_artifact_read(path,limit)
    elseif method=="PUT"
        local_artifact_write(path,bytes);return UInt8[]
    elseif method=="DELETE"
        if ispath(path)||islink(path)
            receipt_private_stat(lstat(path));rm(path)
        end
        return UInt8[]
    end
    throw(ArgumentError("Unsupported artifact operation"))
end

function artifact_object_request(location::S3RuntimeArtifacts,args...;kwargs...)
    try
        return s3_artifact_object_request(location,args...;kwargs...)
    catch
    end
    # Discard SDK/TOML/transport exception contexts before throwing. Otherwise a
    # TaskFailedException can print credentials or a private URL as a root cause.
    throw(ArtifactUnavailable())
end

function s3_artifact_object_request(location::S3RuntimeArtifacts,scope,kind,digest,method,bytes=UInt8[];limit=ASSIGNED_ARTIFACT_BYTES)
    broker_file(location.credentials_file;private=true,max_bytes=8192)
    credentials=TOML.parsefile(location.credentials_file)
    Set(keys(credentials))==Set(("access_key_id","secret_access_key")) || throw(ArtifactUnavailable())
    all(v->v isa String && 0<ncodeunits(v)<=4096 && !any(iscntrl,v),values(credentials)) || throw(ArtifactUnavailable())
    config=S3EndpointConfig(location.endpoint,credentials["access_key_id"],credentials["secret_access_key"];
        region=location.region,allow_insecure=startswith(location.endpoint,"http:"))
    suffix=kind=="metadata" ? "$digest.json" : digest
    key=artifact_storage_key(location.prefix*"/"*scope,kind,suffix)
    url=location.endpoint*"/"*location.bucket*"/"*key
    uri=URIs.URI(location.endpoint)
    authority=uri.host*(isempty(uri.port) ? "" : ":"*uri.port)
    request=AWS.Request(service="s3",api_version="2006-03-01",request_method=method,url=url,
        headers=Dict("Host"=>authority,"Content-Type"=>"application/json"),content=bytes)
    AWS.sign!(config,request) # Reuse the installed SDK; do not implement signing.
    location.ca_file===nothing || broker_file(location.ca_file)
    tls=HTTP.TLS.Config(;ca_file=location.ca_file,verify_peer=true,verify_hostname=true)
    client=HTTP.Client(transport=HTTP.Transport(;tls_config=tls,proxy=nothing,max_idle_per_host=1,max_idle_total=1))
    output=BoundedResponseBody(IOBuffer(),limit)
    try
        response=Logging.with_logger(Logging.NullLogger()) do
            HTTP.request(method,url,collect(pairs(request.headers)),bytes;client,
                response_stream=output,request_timeout=5,connect_timeout=2,read_idle_timeout=2,
                redirect=false,retry=false,cookies=false,decompress=false,status_exception=false)
        end
        response.status==404 && method in ("GET","DELETE") && return nothing
        response.status in (200,201,204) || throw(ArtifactUnavailable())
        take!(output.buffer)
    finally
        close(client)
    end
end

"""Read and hash-check private bytes using the saved job's scope, never a supplied URL."""
function read_job_artifact(location::AbstractRuntimeArtifacts,job::Protocol.AssignedJob,reference::Protocol.ArtifactReference)
    try
        return checked_job_artifact_bytes(location,job,reference)
    catch
    end
    throw(ArtifactUnavailable())
end

function checked_job_artifact_bytes(location::AbstractRuntimeArtifacts,job::Protocol.AssignedJob,reference::Protocol.ArtifactReference)
    checked_artifact_reference(reference)
    scope=artifact_scope(job.fence,job.request.job_id)
    metadata=artifact_object_request(location,scope,"metadata",reference.sha256,"GET";limit=4096)
    metadata===nothing && return nothing
    document=JSON3.read(metadata)
    Set(keys(document))==Set((:sha256,:media_type,:size)) && document.sha256==reference.sha256 &&
        document.media_type==reference.media_type && document.size isa Integer && !(document.size isa Bool) &&
        document.size==reference.size || throw(ArtifactUnavailable())
    bytes=artifact_object_request(location,scope,"sha256",reference.sha256,"GET";limit=reference.size)
    bytes!==nothing && length(bytes)==reference.size && bytes2hex(SHA.sha256(bytes))==reference.sha256 ||
        throw(ArtifactUnavailable())
    bytes
end

prune_artifact_directories(::S3RuntimeArtifacts,scope)=nothing
function prune_artifact_directories(location::LocalRuntimeArtifacts,scope)
    path=private_artifact_directory(location,scope)
    path===nothing && return nothing
    while path!=location.root && isdir(path) && !islink(path) && isempty(readdir(path))
        receipt_private_stat(lstat(path);directory=true)
        rm(path);path=dirname(path)
    end
    nothing
end

"""
    store_job_artifact!(location, job, bytes) -> ArtifactReference

Store bounded JSON bytes in a job-specific namespace using the shared digest and
metadata format. Metadata is written last. A failed write attempts removal of
only this job's exact two objects; unavailable remote cleanup still requires the
operator's private-prefix retention policy. No upload staging file is retained.
The returned reference does not confer access; gateway retrieval requires an
owned receipt and matching durable result.
"""
function store_job_artifact!(location::AbstractRuntimeArtifacts,job::Protocol.AssignedJob,bytes::Vector{UInt8})
    0<length(bytes)<=ASSIGNED_ARTIFACT_BYTES || throw(ArtifactUnavailable())
    digest=bytes2hex(SHA.sha256(bytes))
    reference=Protocol.ArtifactReference("sha256:"*digest,"application/json",length(bytes),digest,
        location isa LocalRuntimeArtifacts ? "local_filesystem" : "s3","/artifacts/sha256/"*digest)
    # Remote agents retain write/delete-only credentials. Immutable PUT retries
    # use the same job scope and digest; the durable job owner prevents replay.
    if location isa LocalRuntimeArtifacts
        read_job_artifact(location,job,reference)===nothing || return reference
    end
    scope=artifact_scope(job.fence,job.request.job_id)
    stored=false
    try
        artifact_object_request(location,scope,"sha256",digest,"PUT",bytes;limit=4096)
        metadata=collect(codeunits(JSON3.write(artifact_metadata_document(digest,"application/json",length(bytes)))))
        artifact_object_request(location,scope,"metadata",digest,"PUT",metadata;limit=4096)
        stored=true
    catch
        for kind in ("metadata","sha256")
            try artifact_object_request(location,scope,kind,digest,"DELETE";limit=4096) catch end
        end
        try prune_artifact_directories(location,scope) catch end
    end
    stored || throw(ArtifactUnavailable())
    reference
end
