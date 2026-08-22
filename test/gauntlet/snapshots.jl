artifact_name(name::Symbol) = "pscad_gauntlet_$(name)"
artifact_name(case::GauntletCase) = artifact_name(case.name)

function _bound_hash(
        case::GauntletCase;
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    isfile(artifacts_toml) || return nothing
    return artifact_hash(artifact_name(case), artifacts_toml)
end

function snapshot_path(
        case::GauntletCase;
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    hash = _bound_hash(case; artifacts_toml)
    hash === nothing && throw(ArgumentError(
        "Gauntlet artifact $(artifact_name(case)) is not bound in $artifacts_toml. " *
        "Record and persist it explicitly before running snapshot mode.",
    ))
    if !artifact_exists(hash)
        try
            ensure_artifact_installed(artifact_name(case), artifacts_toml)
        catch error
            throw(ErrorException(
                "Gauntlet artifact $(artifact_name(case)) is not installed and has no " *
                "usable published download. Original error: $(sprint(showerror, error))",
            ))
        end
    end
    path = joinpath(artifact_path(hash), "snapshot.jld2")
    isfile(path) || throw(ArgumentError(
        "Gauntlet artifact $(artifact_name(case)) has no snapshot.jld2",
    ))
    return path
end

function _write_snapshot(
        path::AbstractString,
        case::GauntletCase,
        reference::LineParameters,
        candidate::LineParameters,
        reference_comparison::LineParametersComparison,
        julia_benchmark;
        pscad_version::AbstractString,
        pscad_elapsed_seconds::Real
)
    validate_structure(case, reference)
    validate_structure(case, candidate)
    mkpath(dirname(path))
    temporary = path * ".new"
    isfile(temporary) && rm(temporary; force = true)
    try
        JLD2.jldsave(
            temporary;
            format_version = SNAPSHOT_FORMAT_VERSION,
            case_name = string(case.name),
            case_sha256 = case_digest(case),
            problem = case.problem,
            formulation = formulation_record(case),
            port_order = case.port_order,
            frequencies = copy(case.problem.frequencies),
            pscad_version = String(pscad_version),
            pscad_formulation = string(case.pscad_formulation),
            reference,
            accepted = candidate,
            reference_comparison,
            pscad_elapsed_seconds = Float64(pscad_elapsed_seconds),
            julia_benchmark,
            recorded_at_utc = string(now(UTC))
        )
        mv(temporary, path; force = true)
    catch
        isfile(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return path
end

_snapshot_digest(path::AbstractString) = bytes2hex(sha256(read(path)))

function _atomic_stage!(temporary::AbstractString, destination::AbstractString)
    backup = destination * ".old"
    ispath(backup) && rm(backup; recursive = true, force = true)
    ispath(destination) && mv(destination, backup; force = true)
    try
        mv(temporary, destination; force = true)
        ispath(backup) && rm(backup; recursive = true, force = true)
    catch
        !ispath(destination) && ispath(backup) && mv(backup, destination; force = true)
        rethrow()
    end
    return destination
end

function persist_snapshot(
        case::GauntletCase,
        reference::LineParameters,
        candidate::LineParameters,
        reference_comparison::LineParametersComparison,
        julia_benchmark;
        pscad_version::AbstractString,
        pscad_elapsed_seconds::Real,
        artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    mkpath(artifact_root)
    destination = joinpath(artifact_root, string(case.name))
    temporary = destination * ".new"
    ispath(temporary) && rm(temporary; recursive = true, force = true)
    mkpath(temporary)
    snapshot = joinpath(temporary, "snapshot.jld2")
    _write_snapshot(
        snapshot,
        case,
        reference,
        candidate,
        reference_comparison,
        julia_benchmark;
        pscad_version,
        pscad_elapsed_seconds
    )
    digest = _snapshot_digest(snapshot)
    write(joinpath(temporary, "snapshot.sha256"), "$digest  snapshot.jld2\n")

    hash = create_artifact() do directory
        cp(snapshot, joinpath(directory, "snapshot.jld2"))
        cp(
            joinpath(temporary, "snapshot.sha256"),
            joinpath(directory, "snapshot.sha256")
        )
    end
    archive = joinpath(temporary, "$(case.name).tar.gz")
    archive_sha256 = archive_artifact(hash, archive)
    mkpath(dirname(artifacts_toml))
    bind_artifact!(
        artifacts_toml,
        artifact_name(case),
        hash;
        lazy = true,
        force = true
    )
    _atomic_stage!(temporary, destination)
    return (
        path = destination,
        snapshot_sha256 = digest,
        archive_sha256,
        tree_hash = string(hash),
        artifact = artifact_name(case)
    )
end

function publish_artifact(
        name::Symbol,
        url::AbstractString;
        artifact_root::AbstractString = ARTIFACT_ROOT,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    isempty(strip(url)) && throw(ArgumentError("artifact download URL cannot be empty"))
    isfile(artifacts_toml) || throw(ArgumentError(
        "Gauntlet Artifacts.toml is missing: $artifacts_toml",
    ))
    artifact = artifact_name(name)
    hash = artifact_hash(artifact, artifacts_toml)
    hash === nothing && throw(ArgumentError(
        "record $artifact before publishing it",
    ))
    archive = joinpath(artifact_root, string(name), "$name.tar.gz")
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
    return (; artifact, tree_hash = string(hash), archive_sha256)
end

function publish_artifact(case::GauntletCase, url::AbstractString; kwargs...)
    publish_artifact(case.name, url; kwargs...)
end

function _load_snapshot_path(case::GauntletCase, path::AbstractString)
    isfile(path) || throw(ArgumentError("Gauntlet snapshot is missing: $path"))
    digest_path = joinpath(dirname(path), "snapshot.sha256")
    isfile(digest_path) || throw(ArgumentError(
        "Gauntlet snapshot digest is missing: $digest_path",
    ))
    expected_digest = first(split(strip(read(digest_path, String))))
    observed_digest = _snapshot_digest(path)
    expected_digest == observed_digest || throw(ArgumentError(
        "Gauntlet snapshot SHA-256 does not match $digest_path",
    ))

    snapshot = JLD2.load(path)
    required = (
        "format_version", "case_name", "case_sha256", "problem", "formulation",
        "port_order", "frequencies", "pscad_version", "pscad_formulation",
        "reference", "accepted", "reference_comparison",
        "pscad_elapsed_seconds", "julia_benchmark", "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["format_version"] == SNAPSHOT_FORMAT_VERSION || throw(ArgumentError(
        "Gauntlet snapshot $path has unsupported format version $(snapshot["format_version"])",
    ))
    snapshot["case_name"] == string(case.name) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to case $(snapshot["case_name"]), not $(case.name)",
    ))
    snapshot["case_sha256"] == case_digest(case) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match $(case.source_file). " *
        "Review the case and record it explicitly.",
    ))
    snapshot["pscad_formulation"] == string(case.pscad_formulation) ||
        throw(ArgumentError(
            "Gauntlet snapshot PSCAD formulation does not match the case definition",
        ))
    snapshot["formulation"] == formulation_record(case) || throw(ArgumentError(
        "Gauntlet snapshot formulation does not match the case definition",
    ))
    snapshot["frequencies"] == case.problem.frequencies || throw(ArgumentError(
        "Gauntlet snapshot frequencies do not match the case definition",
    ))
    snapshot["port_order"] == case.port_order || throw(ArgumentError(
        "Gauntlet snapshot port order does not match the case definition",
    ))
    reference = snapshot["reference"]
    accepted = snapshot["accepted"]
    reference isa LineParameters || throw(ArgumentError(
        "Gauntlet snapshot $path has no reference Engine.LineParameters",
    ))
    accepted isa LineParameters || throw(ArgumentError(
        "Gauntlet snapshot $path has no accepted Engine.LineParameters",
    ))
    validate_structure(case, reference)
    validate_structure(case, accepted)
    return (; reference, accepted, metadata = snapshot)
end

function load_snapshot(
        case::GauntletCase;
        path::Union{Nothing, AbstractString} = nothing,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    resolved = path === nothing ? snapshot_path(case; artifacts_toml) : String(path)
    return _load_snapshot_path(case, resolved)
end

function load_prior_snapshot(
        case::GauntletCase;
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    _bound_hash(case; artifacts_toml) === nothing && return nothing
    return load_snapshot(case; artifacts_toml)
end

function run_snapshot(
        case::GauntletCase;
        path::Union{Nothing, AbstractString} = nothing,
        artifacts_toml::AbstractString = ARTIFACTS_TOML,
        benchmark_samples::Int = 10,
        benchmark_seconds::Real = 10
)
    loaded = load_snapshot(case; path, artifacts_toml)
    candidate = compute!(case.problem, case.formulation)
    validate_structure(case, candidate)
    reference_comparison = compare(loaded.reference, candidate)
    regression_comparison = compare(loaded.accepted, candidate)
    timing = benchmark_local(
        case;
        samples = benchmark_samples,
        seconds = benchmark_seconds
    )
    performance = performance_comparison(
        loaded.metadata["julia_benchmark"],
        timing,
        case.tolerances.performance
    )
    return (
        mode = :snapshot,
        reference = loaded.reference,
        accepted = loaded.accepted,
        candidate,
        comparison = reference_comparison,
        regression = regression_comparison,
        performance,
        metadata = loaded.metadata,
        pscad = nothing,
        timings = (pscad = loaded.metadata["pscad_elapsed_seconds"], julia = timing)
    )
end
