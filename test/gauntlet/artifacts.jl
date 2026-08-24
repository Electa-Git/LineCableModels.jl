module GauntletArtifacts

import Pkg
using DataFrames
using JLD2
using LineCableModels: basis, domain, observables
using LineCableModels.Engine: LineParameters, LineParametersBenchmark, compare
using Pkg.Artifacts
using SHA

export ARTIFACT_ROOT, ARTIFACTS_TOML, GAUNTLET_VERSION,
       artifact_name, backend_archive_name, backend_stage, case_stage,
       cleanup_work, finalize_artifacts, gauntlet_cleanup, gauntlet_force,
       gauntlet_instrumented, gauntlet_mode, prepare_artifacts,
       publish_artifact, release_tag, report

const GAUNTLET_ROOT = @__DIR__
const ARTIFACT_ROOT = joinpath(GAUNTLET_ROOT, ".artifacts")
const ARTIFACTS_TOML = joinpath(GAUNTLET_ROOT, "Artifacts.toml")
const WORK_ROOT = joinpath(GAUNTLET_ROOT, "cases", ".work")
const GAUNTLET_VERSION = v"1.0.0"
const VALID_MODES = (:snapshot, :live, :record)

_version_label() = "v$(GAUNTLET_VERSION)"

release_tag() = "gauntlet-$(_version_label())"

function artifact_name(backend::Symbol)
    version = replace(_version_label(), "." => "_")
    return "gauntlet_$(backend)_$(version)"
end

backend_archive_name(backend::Symbol) = "benchmarks-$(backend)-$(_version_label()).tar.gz"

function backend_stage(
        backend::Symbol;
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    return joinpath(artifact_root, string(backend), _version_label())
end

function case_stage(
        backend::Symbol,
        name::Symbol;
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    return joinpath(backend_stage(backend; artifact_root), "cases", string(name))
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
gauntlet_force() = _boolean_setting("LINECABLEMODELS_GAUNTLET_FORCE")

include("reporting.jl")

function gauntlet_instrumented()
    options = Base.JLOptions()
    return !iszero(options.code_coverage) || !iszero(options.malloc_log)
end

function cleanup_work(; work_root::AbstractString = WORK_ROOT)
    ispath(work_root) && rm(work_root; recursive = true, force = true)
    return work_root
end

function _artifact_bindings(artifacts_toml::AbstractString)
    isfile(artifacts_toml) || return String[]
    document = Pkg.Artifacts.load_artifacts_toml(artifacts_toml)
    suffix = "_" * replace(_version_label(), "." => "_")
    return sort!(
        String[name
               for name in keys(document)
               if startswith(name, "gauntlet_") && endswith(name, suffix)]
    )
end

function _record_collisions(
        ; artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    collisions = String[]
    if isdir(artifact_root)
        for backend in readdir(artifact_root; join = true)
            isdir(backend) || continue
            current = joinpath(backend, _version_label())
            ispath(current) && push!(collisions, current)
        end
    end
    append!(
        collisions,
        ["$artifacts_toml::$name" for name in _artifact_bindings(artifacts_toml)]
    )
    return sort!(unique!(collisions))
end

function _collision_error(collisions)
    listed = join(["  - $collision" for collision in collisions], "\n")
    return ArgumentError(
        "Gauntlet artifacts already exist for $(_version_label()):\n" *
        "$listed\nSet LINECABLEMODELS_GAUNTLET_FORCE=true to replace them.",
    )
end

function prepare_artifacts(
        ; artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML,
        force::Bool = false
)
    collisions = _record_collisions(; artifact_root, artifacts_toml)
    isempty(collisions) || force || throw(_collision_error(collisions))
    if force && isdir(artifact_root)
        for backend in readdir(artifact_root; join = true)
            isdir(backend) || continue
            current = joinpath(backend, _version_label())
            ispath(current) && rm(current; recursive = true, force = true)
        end
    end
    return artifact_root
end

function _staged_backends(; artifact_root::AbstractString = ARTIFACT_ROOT)
    isdir(artifact_root) || return Symbol[]
    backends = Symbol[]
    for entry in readdir(artifact_root)
        cases = joinpath(artifact_root, entry, _version_label(), "cases")
        isdir(cases) && !isempty(readdir(cases)) && push!(backends, Symbol(entry))
    end
    return sort!(backends; by = string)
end

function _finalize_backend(
        backend::Symbol;
        artifact_root::AbstractString = ARTIFACT_ROOT,
        force::Bool = false
)
    stage = backend_stage(backend; artifact_root)
    cases = joinpath(stage, "cases")
    isdir(cases) && !isempty(readdir(cases)) || throw(ArgumentError(
        "no recorded $backend cases exist for $(_version_label())",
    ))
    archive = joinpath(stage, backend_archive_name(backend))
    isfile(archive) && !force && throw(_collision_error([archive]))
    report_files = _write_report(backend, stage)
    hash = create_artifact() do directory
        cp(cases, joinpath(directory, "cases"))
        for path in (
            report_files.jld2_path,
            report_files.tsv_path,
            report_files.digest_path
        )
            cp(path, joinpath(directory, basename(path)))
        end
        write(joinpath(directory, "backend.txt"), "$(backend)\n")
        write(
            joinpath(directory, "gauntlet-version.txt"),
            "$(GAUNTLET_VERSION)\n"
        )
    end
    archive_sha256 = archive_artifact(hash, archive)
    return (
        backend,
        path = stage,
        archive,
        archive_sha256,
        tree_hash = string(hash),
        artifact = artifact_name(backend),
        gauntlet_version = GAUNTLET_VERSION,
        report = report_files.frame
    )
end

function finalize_artifacts(
        ; artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML,
        force::Bool = false
)
    backends = _staged_backends(; artifact_root)
    isempty(backends) && throw(ArgumentError(
        "record mode produced no backend snapshots for $(_version_label())",
    ))
    bindings = Set(_artifact_bindings(artifacts_toml))
    existing = ["$artifacts_toml::$(artifact_name(backend))"
                for backend in backends if artifact_name(backend) in bindings]
    isempty(existing) || force || throw(_collision_error(existing))
    collections = [_finalize_backend(backend; artifact_root, force)
                   for backend in backends]
    mkpath(dirname(artifacts_toml))
    temporary = artifacts_toml * ".new"
    isfile(temporary) && rm(temporary; force = true)
    isfile(artifacts_toml) ? cp(artifacts_toml, temporary) : write(temporary, "")
    try
        for collection in collections
            bind_artifact!(
                temporary,
                collection.artifact,
                Base.SHA1(collection.tree_hash);
                lazy = true,
                force
            )
        end
        mv(temporary, artifacts_toml; force = true)
    catch
        isfile(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return collections
end

function publish_artifact(
        backend::Symbol,
        url::AbstractString;
        artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    isempty(strip(url)) && throw(ArgumentError("artifact download URL cannot be empty"))
    isfile(artifacts_toml) || throw(ArgumentError(
        "Gauntlet Artifacts.toml is missing: $artifacts_toml",
    ))
    artifact = artifact_name(backend)
    hash = artifact_hash(artifact, artifacts_toml)
    hash === nothing && throw(ArgumentError(
        "record the $backend collection for $(_version_label()) before publishing it",
    ))
    archive = joinpath(
        backend_stage(backend; artifact_root),
        backend_archive_name(backend)
    )
    isfile(archive) || throw(ArgumentError("Gauntlet archive is missing: $archive"))
    archive_sha256 = bytes2hex(sha256(read(archive)))
    bind_artifact!(
        artifacts_toml,
        artifact,
        hash;
        download_info = [(String(url), archive_sha256)],
        lazy = true,
        force = true
    )
    return (; backend, artifact, tree_hash = string(hash), archive_sha256)
end

end
