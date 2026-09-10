
import TOML
using JLD2: JLD2, jldopen
import LineCableModels
using LineCableModels: basis, domain, observe
using LineCableModels.Engine: LineParameters, LineParametersBenchmark, RMSError,
                              absolute_error, compare, frequencies,
                              relative_error, Z, Y
using Pkg.Artifacts: Artifacts, archive_artifact, bind_artifact!, create_artifact
using SHA: SHA, sha256
using Dates: Dates, UTC, now
import Downloads

export ARTIFACT_ROOT, ARTIFACTS_TOML, SNAPSHOT_SCHEMA_VERSION,
       artifact_name, benchmark_stage, bind_published_artifact,
       cleanup_work, collection_archive_name, collection_release,
       collection_stage, finalize_staging, gauntlet_instrumented,
       package_collection, prepare_staging, release_tag, read_collection,
       MomentResult, MomentBenchmark, extract_moments,
       moment_comparison_passes, moment_error_summary, read_moments

export BenchmarkCalculation, BenchmarkDefinition,
       benchmark_definition, compare_saved, read_benchmark,
       read_calculation

const GAUNTLET_ROOT = @__DIR__
const ARTIFACT_ROOT = joinpath(GAUNTLET_ROOT, ".artifacts")
const ARTIFACTS_TOML = joinpath(GAUNTLET_ROOT, "Artifacts.toml")
const WORK_ROOT = joinpath(GAUNTLET_ROOT, ".work")
const SNAPSHOT_SCHEMA_VERSION = 2

function _collection_name(collection::Symbol)
    name = string(collection)
    occursin(r"^[a-z][a-z0-9_]*$", name) || throw(ArgumentError(
        "Gauntlet collection names must use lowercase letters, digits, and underscores; " *
        "got $(repr(name))",
    ))
    return name
end

function _release_version(version::VersionNumber)
    isempty(version.prerelease) || throw(ArgumentError(
        "Gauntlet releases cannot use prerelease versions: $version",
    ))
    isempty(version.build) || throw(ArgumentError(
        "Gauntlet releases cannot use build metadata: $version",
    ))
    version >= v"1.0.0" || throw(ArgumentError(
        "Gauntlet release versions start at 1.0.0; got $version",
    ))
    return version
end

artifact_name(collection::Symbol) = "gauntlet_$(_collection_name(collection))"

function release_tag(collection::Symbol, version::VersionNumber)
    return "gauntlet-$(_collection_name(collection))-v$(_release_version(version))"
end

function collection_archive_name(collection::Symbol, version::VersionNumber)
    return "benchmarks-$(_collection_name(collection))-v$(_release_version(version)).tar.gz"
end

function collection_stage(
        collection::Symbol;
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    return joinpath(artifact_root, "staging", _collection_name(collection))
end

function benchmark_stage(
        collection::Symbol,
        benchmark_id::Symbol;
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    return joinpath(
        collection_stage(collection; artifact_root),
        "benchmarks",
        string(benchmark_id)
    )
end

function collection_release(
        collection::Symbol,
        version::VersionNumber;
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    return joinpath(
        artifact_root,
        "releases",
        _collection_name(collection),
        "v$(_release_version(version))"
    )
end

include("comparisons/uq_moments.jl")
include("definitions.jl")
include("fingerprints.jl")
include("comparisons/saved.jl")
include("read.jl")

function gauntlet_instrumented()
    options = Base.JLOptions()
    return !iszero(options.code_coverage) || !iszero(options.malloc_log)
end

function cleanup_work(; work_root::AbstractString = WORK_ROOT)
    validate(Base.write,work_root;recursive=true)
    ispath(work_root) && rm(work_root; recursive = true, force = true)
    return work_root
end

function prepare_staging(
        ; artifact_root::AbstractString = ARTIFACT_ROOT,
        force::Bool = false
)
    staging_root = joinpath(artifact_root, "staging")
    validate(Base.write,staging_root;recursive=true)
    occupied = isdir(staging_root) && !isempty(readdir(staging_root))
    occupied && !force &&
        throw(ArgumentError(
            "Gauntlet staging is not empty: $staging_root\n" *
            "Pass force=true explicitly to replace it.",
        ))
    occupied && rm(staging_root; recursive = true, force = true)
    mkpath(staging_root)
    return staging_root
end

function _staged_collections(; artifact_root::AbstractString = ARTIFACT_ROOT)
    staging_root = joinpath(artifact_root, "staging")
    isdir(staging_root) || return Symbol[]
    collections = Symbol[]
    for entry in readdir(staging_root)
        benchmarks = joinpath(staging_root, entry, "benchmarks")
        isdir(benchmarks) && !isempty(readdir(benchmarks)) &&
            push!(collections, Symbol(entry))
    end
    return sort!(collections; by = string)
end

function finalize_staging(; artifact_root::AbstractString = ARTIFACT_ROOT)
    collections = _staged_collections(; artifact_root)
    isempty(collections) && throw(ArgumentError(
        "record mode produced no staged Gauntlet snapshots",
    ))
    return map(collections) do collection
        stage = collection_stage(collection; artifact_root)
        snapshots = read_collection(stage; collection)
        (
            collection,
            path = stage,
            schema_version = SNAPSHOT_SCHEMA_VERSION,
            benchmarks = length(snapshots)
        )
    end
end

function _write_toml(path::AbstractString, document::AbstractDict)
    validate(Base.write,path)
    mkpath(dirname(path))
    temporary=tempname(dirname(path))
    try
        open(temporary, "w") do io
            TOML.print(io, document; sorted = true)
        end
        mv(temporary, path; force = true)
    finally
        isfile(temporary) && rm(temporary)
    end
    return path
end

"""
    package_collection(collection, version; bundles, reason, output)

Package explicitly accepted bundles as a versioned Julia artifact. Inputs are
verified before archival. Existing identical packages are verified and reused;
a different release definition requires another version. No computation,
illustration, upload or Git operation is performed.
"""
function package_collection(collection::Symbol,version::VersionNumber;
        bundles,reason::AbstractString,output::AbstractString)
    name=_collection_name(collection)
    release_version=_release_version(version)
    isempty(strip(reason)) && throw(ArgumentError("release description cannot be empty"))
    isempty(bundles) && throw(ArgumentError("package needs explicit accepted bundles"))
    inputs=Dict{String,Any}[]
    seen=Set{Tuple{String,String}}()
    for path in bundles
        root=abspath(path)
        isfile(joinpath(root,"bundle.toml")) || throw(ArgumentError("packaging requires locked bundles; drafts cannot be packaged"))
        records=read_campaign(root)
        document=TOML.parsefile(joinpath(root,"bundle.toml"))
        for record in records
            key=(first(record.analyses)["case_id"],string(record.id))
            key in seen && throw(ArgumentError("duplicate or conflicting benchmark in release: $key"))
            push!(seen,key)
        end
        push!(inputs,Dict("identity"=>document["identity"],"path"=>root))
    end
    allunique(entry["identity"] for entry in inputs) || throw(ArgumentError("duplicate accepted bundle"))
    identities=sort!([entry["identity"] for entry in inputs])
    release=Dict("schema"=>3,"collection"=>name,"version"=>string(release_version),
        "description"=>String(reason),"bundles"=>identities,"tag"=>release_tag(collection,release_version))
    signature=semantic_sha256(release)
    destination=abspath(output)
    validate(Base.write,destination)
    if ispath(destination)
        path=joinpath(destination,"package.toml")
        isfile(path) || throw(ArgumentError("release destination already exists"))
        existing=TOML.parsefile(path)
        existing["signature"] == signature || throw(ArgumentError("release contents changed; choose a new version"))
        archive=joinpath(destination,existing["artifact"]["archive"])
        bytes2hex(open(sha256,archive)) == existing["artifact"]["archive_sha256"] || throw(ArgumentError("existing release archive changed"))
        return (collection,version=release_version,path=destination,archive,
            archive_sha256=existing["artifact"]["archive_sha256"],tree_hash=existing["artifact"]["tree_hash"],package_path=path)
    end
    hash=create_artifact() do directory
        for entry in inputs
            target=joinpath(directory,"bundles",entry["identity"])
            mkpath(dirname(target))
            cp(entry["path"],target)
            read_campaign(target)
        end
        _write_toml(joinpath(directory,"release.toml"),release)
        read_collection(directory;collection)
    end
    mkpath(dirname(destination))
    staging=mktempdir(dirname(destination))
    try
        archive_name=collection_archive_name(collection,release_version)
        archive_sha256=archive_artifact(hash,joinpath(staging,archive_name))
        # Read the actual archived files rather than relying only on the source tree.
        mktempdir() do extracted
            unpack(joinpath(staging,archive_name),extracted)
            read_collection(extracted;collection)
        end
        artifact=Dict("archive"=>archive_name,"archive_sha256"=>archive_sha256,
            "tree_hash"=>string(hash),"name"=>artifact_name(collection)*"_v"*replace(string(release_version),'.'=>'_'))
        _write_toml(joinpath(staging,"package.toml"),Dict("signature"=>signature,"release"=>release,"artifact"=>artifact))
        mv(staging,destination)
        return (collection,version=release_version,path=destination,
            archive=joinpath(destination,archive_name),archive_sha256,tree_hash=string(hash),
            package_path=joinpath(destination,"package.toml"))
    finally
        isdir(staging) && rm(staging;recursive=true)
    end
end

"""Read an explicit TOML release definition and package only its selected accepted bundles."""
function package_collection(definition::AbstractString;output::AbstractString)
    source=abspath(definition)
    document=TOML.parsefile(source)
    Set(keys(document)) == Set(("collection","version","description","bundles")) ||
        throw(ArgumentError("release definition requires collection, version, description and bundles"))
    bundles=[normpath(joinpath(dirname(source),path)) for path in document["bundles"]]
    return package_collection(Symbol(document["collection"]),VersionNumber(document["version"]);
        bundles,reason=document["description"],output)
end

"""
    bind_published_artifact(package, url; artifacts_toml=ARTIFACTS_TOML, current=false)

Verify an uploaded archive and bind its immutable version-specific download.
The served archive checksum and unpacked artifact tree must match the package.
Updating the convenience current-version binding is explicit. This does not upload.
"""
function bind_published_artifact(package::AbstractString,url::AbstractString;
        artifacts_toml::AbstractString=ARTIFACTS_TOML,current::Bool=false)
    validate(Base.write,artifacts_toml)
    path=isdir(package) ? joinpath(package,"package.toml") : abspath(package)
    document=TOML.parsefile(path)
    artifact=document["artifact"]
    release=document["release"]
    semantic_sha256(release)==document["signature"] || throw(ArgumentError("release definition changed"))
    expected_name=artifact_name(Symbol(release["collection"]))*"_v"*replace(release["version"],'.'=>'_')
    artifact["name"]==expected_name || throw(ArgumentError("artifact name differs from the release definition"))
    archive=joinpath(dirname(path),artifact["archive"])
    bytes2hex(open(sha256,archive)) == artifact["archive_sha256"] || throw(ArgumentError("local archive changed"))
    isempty(strip(url)) && throw(ArgumentError("published URL cannot be empty"))
    name=artifact["name"]
    mkpath(dirname(abspath(artifacts_toml)))
    lease=open(artifacts_toml*".lock","a+")
    acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
        ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
    acquired || (close(lease);throw(ArgumentError("another process is updating artifact bindings")))
    try
        prior=isfile(artifacts_toml) ? TOML.parsefile(artifacts_toml) : Dict()
        if haskey(prior,name)
            prior[name]["git-tree-sha1"] == artifact["tree_hash"] &&
                all(item -> item["sha256"] == artifact["archive_sha256"],get(prior[name],"download",[])) ||
                throw(ArgumentError("published collection/version cannot be rebound to different contents"))
        end
        mktempdir() do directory
            downloaded=joinpath(directory,"download.tar.gz")
            Downloads.download(String(url),downloaded)
            bytes2hex(open(sha256,downloaded)) == artifact["archive_sha256"] || throw(ArgumentError("published archive checksum mismatch"))
            hash=create_artifact() do extracted
                unpack(downloaded,extracted)
                TOML.parsefile(joinpath(extracted,"release.toml"))==release ||
                    throw(ArgumentError("served release definition differs from the package"))
                read_collection(extracted;collection=Symbol(release["collection"]))
            end
            string(hash) == artifact["tree_hash"] || throw(ArgumentError("published artifact tree mismatch"))
            mkpath(dirname(abspath(artifacts_toml)))
            temporary=tempname(dirname(abspath(artifacts_toml)))
            try
                isfile(artifacts_toml) && cp(artifacts_toml,temporary)
                bind_artifact!(temporary,name,hash;download_info=[(String(url),artifact["archive_sha256"])],lazy=true,force=true)
                if current
                    bind_artifact!(temporary,artifact_name(Symbol(release["collection"])),hash;
                        download_info=[(String(url),artifact["archive_sha256"])],lazy=true,force=true)
                end
                mv(temporary,artifacts_toml;force=true)
            finally
                isfile(temporary) && rm(temporary)
            end
        end
    finally
        close(lease)
    end
    return (artifact=name,version=release["version"],tree_hash=artifact["tree_hash"],
        archive_sha256=artifact["archive_sha256"],url=String(url))
end

"""
    lock_campaign(directory, destination; benchmarks=nothing, expected=nothing, note="", illustrations=())

Accept selected complete benchmarks into a new checksummed bundle. Omitting
`benchmarks` requires every declared campaign member. `expected` can name the
inspected snapshot identity for a single benchmark. No calculation or plot runs.
"""
function lock_campaign(directory::AbstractString,destination::AbstractString;
        benchmarks=nothing,expected=nothing,note::AbstractString="",illustrations=())
    source=abspath(directory)
    isfile(joinpath(source,"bundle.toml")) && throw(ArgumentError("source is already a locked bundle; package that bundle directly"))
    target=abspath(destination)
    validate(Base.write,target)
    ispath(target) && throw(ArgumentError("bundle destination already exists: $target"))
    (target == source || startswith(target,source*Base.Filesystem.path_separator)) &&
        throw(ArgumentError("bundle destination must be outside its campaign"))
    manifest=TOML.parsefile(joinpath(source,"campaign.toml"))
    declared=String.(manifest["benchmarks"])
    selected=benchmarks === nothing ? declared : benchmarks isa Union{Symbol,AbstractString} ? [string(benchmarks)] : string.(benchmarks)
    !isempty(selected) && allunique(selected) && all(id -> id in declared,selected) ||
        throw(ArgumentError("select distinct declared benchmark IDs"))
    expected === nothing || length(selected)==1 || throw(ArgumentError("expected snapshot identity requires one selected benchmark"))
    leases=IO[]
    staging=nothing
    try
        for id in sort(selected)
            isfile(joinpath(source,id,"state.toml")) || throw(ArgumentError("benchmark $id has no completed attempt"))
            path=manifest["schema"] == 3 ? joinpath(source,id,"execution.lock") : joinpath(source,"execution.lock")
            manifest["schema"] == 3 || isempty(leases) || continue
            lease=open(path,"a+")
            acquired=Sys.iswindows() ? ccall(:_locking,Cint,(Cint,Cint,Clong),fd(lease),2,1)==0 :
                ccall(:flock,Cint,(Cint,Cint),fd(lease),6)==0
            acquired || (close(lease);throw(ArgumentError("another process owns benchmark $id")))
            push!(leases,lease)
        end
        roots=Dict{String,String}()
        identities=Dict{String,String}()
        inventories=Dict{String,Vector{String}}()
        for id in selected
            state_path=joinpath(source,id,"state.toml")
            isfile(state_path) || throw(ArgumentError("benchmark $id has no completed attempt"))
            state=TOML.parsefile(state_path)
            state["state"] == "complete" || throw(ArgumentError("only completed benchmarks can be locked: $id"))
            root=manifest["schema"] == 3 ? joinpath(source,id,state["current"]) : joinpath(source,id)
            loaded=read_benchmark(root)
            files=String[]
            declaration=joinpath(root,"declarations.jld2")
            if isfile(declaration)
                bytes2hex(open(sha256,declaration)) == strip(read(declaration*".sha256",String)) ||
                    throw(ArgumentError("campaign declaration integrity check failed"))
                append!(files,[declaration,declaration*".sha256"])
            end
            for role in ("reference","candidate")
                operand=joinpath(root,role,"calculation.jld2")
                append!(files,[operand,operand*".sha256"])
                marker=joinpath(dirname(operand),"complete.toml")
                isfile(marker) && push!(files,marker)
                evidence=jldopen(operand,"r") do file
                    haskey(file,"retained_files") ? file["retained_files"] : ()
                end
                append!(files,[joinpath(dirname(operand),entry.path) for entry in evidence])
            end
            for (folder,_,names) in walkdir(joinpath(root,"analyses"))
                "snapshot.jld2" in names || continue
                snapshot=joinpath(folder,"snapshot.jld2")
                read_benchmark(snapshot)
                append!(files,[snapshot,snapshot*".sha256"])
            end
            identity=semantic_sha256(read_benchmark,root)
            expected === nothing || identity == expected || throw(ArgumentError("draft differs from the inspected snapshot identity"))
            roots[id]=root
            identities[id]=identity
            inventories[id]=sort!(unique(files))
        end
        mkpath(dirname(target))
        staging=mktempdir(dirname(target))
        files=Dict{String,String}()
        for id in selected
            for path in inventories[id]
                relative=relpath(path,roots[id])
                first(splitpath(relative)) == ".." && throw(ArgumentError("bundle dependency escapes benchmark"))
                islink(path) && throw(ArgumentError("bundle evidence must be regular files"))
                target_path=joinpath(staging,id,relative)
                mkpath(dirname(target_path))
                cp(path,target_path)
                digest=bytes2hex(open(sha256,path))
                bytes2hex(open(sha256,target_path)) == digest || throw(ArgumentError("benchmark changed during locking"))
                files[relpath(target_path,staging)]=digest
            end
            state_path=joinpath(staging,id,"state.toml")
            _write_toml(state_path,Dict("state"=>"complete","identity"=>identities[id]))
            files[relpath(state_path,staging)]=bytes2hex(open(sha256,state_path))
        end
        if manifest["schema"] == 2
            for name in ("declarations.jld2","declarations.jld2.sha256")
                path=joinpath(source,name)
                isfile(path) || continue
                cp(path,joinpath(staging,name))
                files[name]=bytes2hex(open(sha256,path))
            end
        end
        figures=Dict{String,Any}[]
        for (index,entry) in enumerate(illustrations)
            string(entry.benchmark) in selected || throw(ArgumentError("illustration names an unselected benchmark"))
            haskey(entry,:selection) && entry.selection isa Union{NamedTuple,AbstractDict} && !isempty(entry.selection) ||
                throw(ArgumentError("illustration requires its explicit problem, quantities and display selections"))
            selection=Dict(string(key)=>value for (key,value) in pairs(entry.selection))
            all(key -> haskey(selection,key),("problem","quantities")) ||
                throw(ArgumentError("illustration selection requires problem and quantities"))
            isfile(entry.path) || throw(ArgumentError("illustration is missing: $(entry.path)"))
            relative=joinpath("illustrations",string(index)*splitext(entry.path)[2])
            mkpath(dirname(joinpath(staging,relative)))
            cp(entry.path,joinpath(staging,relative))
            files[relative]=bytes2hex(open(sha256,joinpath(staging,relative)))
            push!(figures,Dict("path"=>relative,"benchmark"=>string(entry.benchmark),"caption"=>entry.caption,"selection"=>selection))
        end
        campaign=joinpath(staging,"campaign.toml")
        _write_toml(campaign,Dict("schema"=>2,"benchmarks"=>selected))
        files["campaign.toml"]=bytes2hex(open(sha256,campaign))
        bundle=Dict("schema"=>2,"files"=>files,"snapshots"=>identities,
            "accepted"=>string(now(UTC)),"note"=>String(note),"illustrations"=>figures)
        identity=semantic_sha256(bundle)
        bundle["identity"]=identity
        _write_toml(joinpath(staging,"bundle.toml"),bundle)
        read_campaign(staging)
        mv(staging,target)
        return (;path=target,identity)
    finally
        staging === nothing || (isdir(staging) && rm(staging;recursive=true))
        foreach(close,leases)
    end
end

"""Read and verify complete campaign results without loading declarations or solvers."""
function read_campaign(directory::AbstractString)
    root=abspath(directory)
    bundle=joinpath(root,"bundle.toml")
    if isfile(bundle)
        document=TOML.parsefile(bundle)
        document["schema"] in (1,2) || throw(ArgumentError("unsupported bundle schema"))
        recorded_identity=document["identity"]
        inventory=copy(document)
        delete!(inventory,"identity")
        semantic_sha256(document["schema"] == 1 ? document["files"] : inventory) == recorded_identity ||
            throw(ArgumentError("bundle inventory changed"))
        actual=String[]
        for (folder,children,names) in walkdir(root)
            any(child -> islink(joinpath(folder,child)),children) && throw(ArgumentError("bundle directories must not be symlinks"))
            for name in names
                relative=relpath(joinpath(folder,name),root)
                relative == "bundle.toml" || push!(actual,relative)
            end
        end
        Set(actual) == Set(keys(document["files"])) || throw(ArgumentError("bundle inventory contains missing or unexpected files"))
        for (relative,digest) in document["files"]
            (isabspath(relative) || first(splitpath(normpath(relative))) == "..") && throw(ArgumentError("invalid bundle path"))
            path=joinpath(root,relative)
            isfile(path) && !islink(path) && bytes2hex(open(sha256,path)) == digest || throw(ArgumentError("bundle file is missing or changed: $relative"))
        end
    end
    manifest=TOML.parsefile(joinpath(root,"campaign.toml"))
    manifest["schema"] in (2,3) || throw(ArgumentError("unsupported campaign schema"))
    return map(manifest["benchmarks"]) do id
        path=joinpath(root,id)
        state=TOML.parsefile(joinpath(path,"state.toml"))
        state["state"] == "complete" || throw(ArgumentError("benchmark is incomplete: $id"))
        read_benchmark(path)
    end
end

export lock_campaign,read_campaign

"""Reject writes into accepted bundles; recursive deletion also checks contained bundles."""
function validate(::typeof(Base.write),path::AbstractString;recursive::Bool=false)
    if recursive && isdir(path)
        for (directory,_,files) in walkdir(path)
            "bundle.toml" in files && throw(ArgumentError(
                "locked bundles are immutable; recursive deletion would remove an accepted bundle"))
        end
    end
    directory=abspath(path)
    while true
        if isdir(directory)
            directory=realpath(directory)
            isfile(joinpath(directory,"bundle.toml")) && throw(ArgumentError(
                "locked bundles are immutable; write a new draft outside the vault"))
        end
        parent=dirname(directory)
        parent == directory && break
        directory=parent
    end
    return nothing
end
