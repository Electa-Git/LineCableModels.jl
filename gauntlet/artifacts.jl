
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

export ARTIFACT_ROOT, ARTIFACTS_TOML, SNAPSHOT_SCHEMA_VERSION,
       artifact_name, benchmark_stage, bind_published_artifact,
       cleanup_work, collection_archive_name, collection_release,
       collection_stage, finalize_staging, gauntlet_instrumented,
       package_collection, prepare_staging, release_tag, read_collection,
       MomentResult, MomentBenchmark, extract_moments,
       moment_comparison_passes, moment_error_summary, read_moments

export BenchmarkCalculation, BenchmarkDefinition,
       benchmark_definition, compare_saved, read_benchmark,
       read_calculation, benchmark_comparisons

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
    ispath(work_root) && rm(work_root; recursive = true, force = true)
    return work_root
end

function prepare_staging(
        ; artifact_root::AbstractString = ARTIFACT_ROOT,
        force::Bool = false
)
    staging_root = joinpath(artifact_root, "staging")
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

function _package_collision(path::AbstractString)
    return ArgumentError(
        "Gauntlet release package already exists: $path\n" *
        "Pass --force to the external packaging command to replace this local package.",
    )
end

function _write_toml(path::AbstractString, document::AbstractDict)
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
    package_collection(collection::Symbol, version::VersionNumber; reason, git_commit, kwargs...)

Archive completed numerical records, their comparisons and retained evidence. Packaging does not run calculations, publish files, create Git
tags or approve numerical references.

# Arguments

- `collection`: Lowercase collection identifier.
- `version`: Release version, starting at `v"1.0.0"`.

# Keywords

- `reason`: Nonempty release description.
- `git_commit`: Full Git object ID of the packaging checkout.
- `artifact_root`: Local staging and release root.
- `force`: Replace an existing local package; defaults to `false`.

# Returns

- Archive path, SHA-256, artifact tree hash and package metadata path.

# Errors

Invalid release identifiers, missing or corrupted snapshots, and existing
packages without `force=true` raise `ArgumentError`.
"""
function package_collection(
        collection::Symbol,
        version::VersionNumber;
        reason::AbstractString,
        git_commit::AbstractString,
        artifact_root::AbstractString = ARTIFACT_ROOT,
        force::Bool = false
)
    name = _collection_name(collection)
    release_version = _release_version(version)
    release_reason = strip(reason)
    isempty(release_reason) && throw(ArgumentError("release reason cannot be empty"))
    commit = lowercase(strip(git_commit))
    occursin(r"^[0-9a-f]{40}([0-9a-f]{24})?$", commit) || throw(ArgumentError(
        "git_commit must be a full hexadecimal Git object ID",
    ))
    stage = collection_stage(collection; artifact_root)
    read_collection(stage; collection)
    destination = collection_release(collection, release_version; artifact_root)
    ispath(destination) && !force && throw(_package_collision(destination))
    tag = release_tag(collection, release_version)
    release_document = Dict(
        "collection" => name,
        "git_commit" => commit,
        "reason" => release_reason,
        "schema_version" => SNAPSHOT_SCHEMA_VERSION,
        "tag" => tag,
        "version" => string(release_version)
    )
    hash = create_artifact() do directory
        # Retain whole modern bundles; old snapshots remain byte-for-byte intact.
        for entry in
            sort!(filter(isdir, readdir(joinpath(stage, "benchmarks"); join = true)))
            target=joinpath(directory, "benchmarks", basename(entry))
            mkpath(dirname(target))
            cp(entry, target; follow_symlinks = true)
        end
        read_collection(directory; collection)
        _write_toml(joinpath(directory, "release.toml"), release_document)
    end
    mkpath(dirname(destination))
    staging = mktempdir(dirname(destination))
    try
        archive_name = collection_archive_name(collection, release_version)
        archive_sha256 = archive_artifact(hash, joinpath(staging, archive_name))
        package_document = Dict(
            "artifact" => Dict(
                "archive" => archive_name,
                "archive_sha256" => archive_sha256,
                "name" => artifact_name(collection),
                "tree_hash" => string(hash)
            ),
            "release" => release_document
        )
        _write_toml(joinpath(staging, "package.toml"), package_document)
        # Only an explicitly forced local package may be replaced. A failed
        # archive operation never leaves a release that appears complete.
        force && ispath(destination) && rm(destination; recursive = true)
        mv(staging, destination; force = false)
        return (
            collection,
            version = release_version,
            tag,
            reason = release_reason,
            git_commit = commit,
            path = destination,
            archive = joinpath(destination, archive_name),
            archive_sha256,
            tree_hash = string(hash),
            artifact = artifact_name(collection),
            package_path = joinpath(destination, "package.toml")
        )
    finally
        isdir(staging) && rm(staging; recursive = true)
    end
end

function bind_published_artifact(
        collection::Symbol,
        version::VersionNumber,
        url::AbstractString;
        artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    download_url = strip(url)
    isempty(download_url) && throw(ArgumentError("artifact download URL cannot be empty"))
    destination = collection_release(collection, version; artifact_root)
    package_path = joinpath(destination, "package.toml")
    isfile(package_path) || throw(ArgumentError(
        "Gauntlet package metadata is missing: $package_path",
    ))
    package = TOML.parsefile(package_path)
    release = package["release"]
    artifact = package["artifact"]
    release["collection"] == _collection_name(collection) || throw(ArgumentError(
        "Gauntlet package collection does not match $collection",
    ))
    release["version"] == string(_release_version(version)) || throw(ArgumentError(
        "Gauntlet package version does not match $version",
    ))
    archive = joinpath(destination, artifact["archive"])
    isfile(archive) || throw(ArgumentError("Gauntlet archive is missing: $archive"))
    archive_sha256 = bytes2hex(sha256(read(archive)))
    archive_sha256 == artifact["archive_sha256"] || throw(ArgumentError(
        "Gauntlet archive digest does not match $package_path",
    ))
    name = artifact_name(collection)
    artifact["name"] == name || throw(ArgumentError(
        "Gauntlet package artifact name does not match $name",
    ))
    mkpath(dirname(artifacts_toml))
    bind_artifact!(
        artifacts_toml,
        name,
        Base.SHA1(artifact["tree_hash"]);
        download_info = [(String(download_url), archive_sha256)],
        lazy = true,
        force = true
    )
    return (
        collection,
        version = _release_version(version),
        tag = release["tag"],
        artifact = name,
        tree_hash = artifact["tree_hash"],
        archive_sha256,
        url = String(download_url)
    )
end

"""
    lock_campaign(directory, destination)

Copy a completed campaign, its numerical operands, analyses, runtime source bytes
and backend-declared raw files into an immutable checksummed bundle. All internal
file bindings remain relative to the bundle. Publication is a separate action.
"""
function lock_campaign(directory::AbstractString, destination::AbstractString)
    source=abspath(directory);
    target=abspath(destination)
    ispath(target) && throw(ArgumentError("bundle destination already exists: $target"))
    (target == source || startswith(target, source*Base.Filesystem.path_separator)) &&
        throw(ArgumentError("bundle destination must be outside its campaign"))
    lease=open(joinpath(source, "execution.lock"), "a+")
    acquired=Sys.iswindows() ?
             ccall(:_locking, Cint, (Cint, Cint, Clong), fd(lease), 2, 1) == 0 :
             ccall(:flock, Cint, (Cint, Cint), fd(lease), 6) == 0
    acquired ||
        (close(lease); throw(ArgumentError("another process owns campaign $source")))
    staging=nothing
    try
        states=campaign_status(source)
        !isempty(states) && all(row -> row.state === :complete, states) ||
            throw(ArgumentError("only completed campaigns can be locked"))
        read_campaign(source)
        mkpath(dirname(target))
        staging=mktempdir(dirname(target))
        files=Dict{String, String}()
        for (directory, children, names) in walkdir(source)
            any(child -> islink(joinpath(directory, child)), children) &&
                throw(ArgumentError("bundle directories must not be symlinks"))
            for name in sort(names)
                name in ("execution.lock", "bundle.toml") && continue
                path=joinpath(directory, name)
                islink(path) &&
                    throw(ArgumentError("bundle evidence must be files, not symlinks: $path"))
                relative=relpath(path, source)
                output=joinpath(staging, relative)
                mkpath(dirname(output));
                cp(path, output)
                digest=bytes2hex(open(sha256, path))
                bytes2hex(open(sha256, output)) == digest ||
                    throw(ArgumentError("campaign changed during locking"))
                files[relative]=digest
            end
        end
        identity=semantic_sha256(files)
        _write_toml(joinpath(staging, "bundle.toml"),
            Dict("schema"=>1, "identity"=>identity, "files"=>files))
        read_campaign(staging)
        mv(staging, target)
        return (; path = target, identity)
    finally
        staging === nothing || (isdir(staging) && rm(staging; recursive = true))
        close(lease)
    end
end

"""
    read_campaign(directory)

Read completed numerical operands and their comparisons without loading a
declaration constructor or starting a solver. Verify a locked bundle's full
inventory before reading individual checksummed calculation payloads.
"""
function read_campaign(directory::AbstractString)
    root=abspath(directory)
    bundle=joinpath(root, "bundle.toml")
    if isfile(bundle)
        document=TOML.parsefile(bundle)
        document["schema"] == 1 || throw(ArgumentError("unsupported bundle schema"))
        semantic_sha256(document["files"]) == document["identity"] ||
            throw(ArgumentError("bundle inventory changed"))
        for (relative, digest) in document["files"]
            isabspath(relative) || first(splitpath(normpath(relative))) == ".." ?
            throw(ArgumentError("invalid bundle path: $relative")) : nothing
            path=joinpath(root, relative)
            isfile(path) && !islink(path) && bytes2hex(open(sha256, path)) == digest ||
                throw(ArgumentError("bundle file is missing or changed: $relative"))
        end
    end
    manifest=TOML.parsefile(joinpath(root, "campaign.toml"))
    manifest["schema"] == 2 || throw(ArgumentError("unsupported campaign schema"))
    return map(manifest["benchmarks"]) do id
        path=joinpath(root, id)
        state=TOML.parsefile(joinpath(path, "state.toml"))
        state["state"] == "complete" || throw(ArgumentError("benchmark is incomplete: $id"))
        read_benchmark(path)
    end
end

export lock_campaign, read_campaign
