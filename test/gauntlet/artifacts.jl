module GauntletArtifacts

import TOML
using DataFrames
using JLD2
using LineCableModels: basis, domain, observe, observables
using LineCableModels.Engine: LineParameters, LineParametersBenchmark,
                              absolute_error, compare, frequencies,
                              relative_error, Z, Y
using Pkg.Artifacts
using SHA

export ARTIFACT_ROOT, ARTIFACTS_TOML, SNAPSHOT_SCHEMA_VERSION,
       artifact_name, benchmark_stage, bind_published_artifact, case_stage,
       cleanup_work, collection_archive_name, collection_release,
       collection_stage, finalize_staging, gauntlet_cleanup,
       gauntlet_instrumented, gauntlet_mode, gauntlet_stage_force,
       package_collection, prepare_staging, release_tag, report

const GAUNTLET_ROOT = @__DIR__
const ARTIFACT_ROOT = joinpath(GAUNTLET_ROOT, ".artifacts")
const ARTIFACTS_TOML = joinpath(GAUNTLET_ROOT, "Artifacts.toml")
const WORK_ROOT = joinpath(GAUNTLET_ROOT, "benchmarks", ".work")
const SNAPSHOT_SCHEMA_VERSION = 2
const VALID_MODES = (:snapshot, :live, :record)

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

case_stage(collection::Symbol, benchmark_id::Symbol; kwargs...) =
    benchmark_stage(collection, benchmark_id; kwargs...)

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

function gauntlet_mode()
    value = Symbol(get(ENV, "LINECABLEMODELS_GAUNTLET_MODE", "snapshot"))
    value in VALID_MODES || throw(ArgumentError(
        "LINECABLEMODELS_GAUNTLET_MODE must be snapshot, live, or record; got $(repr(value))",
    ))
    haskey(ENV, "CI") && value !== :snapshot &&
        throw(ArgumentError(
            "CI permits only LINECABLEMODELS_GAUNTLET_MODE=snapshot",
        ))
    return value
end

function _boolean_setting(name::AbstractString)
    value = lowercase(strip(get(ENV, name, "false")))
    value in ("true", "false") || throw(ArgumentError(
        "$name must be true or false; got $(repr(value))",
    ))
    return value == "true"
end

gauntlet_cleanup() = _boolean_setting("LINECABLEMODELS_GAUNTLET_CLEANUP")
gauntlet_stage_force() =
    _boolean_setting("LINECABLEMODELS_GAUNTLET_STAGE_FORCE")

include("reporting.jl")

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
    occupied && !force && throw(ArgumentError(
        "Gauntlet staging is not empty: $staging_root\n" *
        "Set LINECABLEMODELS_GAUNTLET_STAGE_FORCE=true to replace it.",
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
        report_files = _write_report(collection, stage)
        (
            collection,
            path = stage,
            schema_version = SNAPSHOT_SCHEMA_VERSION,
            report = report_files.frame,
            report_files
        )
    end
end

function _package_collision(path::AbstractString)
    return ArgumentError(
        "Gauntlet release package already exists: $path\n" *
        "Pass --force to the external packaging command to replace this local package.",
    )
end

function _validate_staged_report(collection::Symbol, stage::AbstractString)
    report_files = (
        jld2_path = joinpath(stage, "report.jld2"),
        tsv_path = joinpath(stage, "report.tsv"),
        digest_path = joinpath(stage, "report.sha256")
    )
    all(isfile, report_files) || throw(ArgumentError(
        "the staged $collection collection has not been finalized: $stage",
    ))
    expected = Dict{String, String}()
    for line in eachline(report_files.digest_path)
        fields = split(strip(line); limit = 2)
        length(fields) == 2 || throw(ArgumentError(
            "invalid staged report digest: $(report_files.digest_path)",
        ))
        expected[strip(fields[2])] = fields[1]
    end
    for path in (report_files.jld2_path, report_files.tsv_path)
        haskey(expected, basename(path)) || throw(ArgumentError(
            "staged report digest does not cover $(basename(path))",
        ))
        expected[basename(path)] == bytes2hex(sha256(read(path))) ||
            throw(ArgumentError("staged report digest does not match $path"))
    end
    report(stage; backend = collection)
    return report_files
end

function _write_toml(path::AbstractString, document::AbstractDict)
    temporary = path * ".new"
    isfile(temporary) && rm(temporary; force = true)
    open(temporary, "w") do io
        TOML.print(io, document; sorted = true)
    end
    mv(temporary, path; force = true)
    return path
end

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
    benchmarks = joinpath(stage, "benchmarks")
    isdir(benchmarks) && !isempty(readdir(benchmarks)) || throw(ArgumentError(
        "no staged $collection benchmarks exist: $benchmarks",
    ))
    report_files = _validate_staged_report(collection, stage)
    destination = collection_release(collection, release_version; artifact_root)
    ispath(destination) && !force && throw(_package_collision(destination))
    ispath(destination) && rm(destination; recursive = true, force = true)
    mkpath(destination)
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
        cp(benchmarks, joinpath(directory, "benchmarks"))
        for path in values(report_files)
            cp(path, joinpath(directory, basename(path)))
        end
        _write_toml(joinpath(directory, "release.toml"), release_document)
    end
    archive = joinpath(
        destination,
        collection_archive_name(collection, release_version)
    )
    archive_sha256 = archive_artifact(hash, archive)
    package_document = Dict(
        "artifact" => Dict(
            "archive" => basename(archive),
            "archive_sha256" => archive_sha256,
            "name" => artifact_name(collection),
            "tree_hash" => string(hash)
        ),
        "release" => release_document
    )
    package_path = _write_toml(
        joinpath(destination, "package.toml"),
        package_document
    )
    return (
        collection,
        version = release_version,
        tag,
        reason = release_reason,
        git_commit = commit,
        path = destination,
        archive,
        archive_sha256,
        tree_hash = string(hash),
        artifact = artifact_name(collection),
        package_path
    )
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

end
