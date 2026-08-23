function _bound_hash(
        case::GauntletCase;
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    isfile(artifacts_toml) || return nothing
    return artifact_hash(artifact_name(case.backend), artifacts_toml)
end

function snapshot_path(
        case::GauntletCase;
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    artifact = artifact_name(case.backend)
    hash = _bound_hash(case; artifacts_toml)
    hash === nothing && throw(ArgumentError(
        "Gauntlet collection $artifact is not bound in $artifacts_toml. " *
        "Record and publish $(release_tag()) explicitly before running snapshot mode.",
    ))
    if !artifact_exists(hash)
        try
            ensure_artifact_installed(artifact, artifacts_toml)
        catch error
            throw(ErrorException(
                "Gauntlet collection $artifact is not installed and has no usable " *
                "download. Original error: $(sprint(showerror, error))",
            ))
        end
    end
    path = joinpath(
        artifact_path(hash),
        "cases",
        string(case.name),
        "snapshot.jld2"
    )
    isfile(path) || throw(ArgumentError(
        "Gauntlet collection $artifact has no snapshot for $(case.name)",
    ))
    return path
end

function _write_snapshot(
        path::AbstractString,
        case::GauntletCase,
        reference::LineParameters,
        candidate::LineParameters,
        reference_comparison::LineParametersBenchmark,
        julia_benchmark;
        reference_execution::NamedTuple
)
    validate_structure(case, reference)
    validate_structure(case, candidate)
    reference_execution.backend === case.backend || throw(ArgumentError(
        "reference execution backend $(reference_execution.backend) does not match " *
        "case backend $(case.backend)",
    ))
    mkpath(dirname(path))
    temporary = path * ".new"
    isfile(temporary) && rm(temporary; force = true)
    execution = _execution_record(reference_execution)
    try
        JLD2.jldsave(
            temporary;
            gauntlet_version = string(GAUNTLET_VERSION),
            case_name = string(case.name),
            backend = string(case.backend),
            case_sha256 = case_digest(case),
            problem = (
                reference = case.reference_problem,
                candidate = case.problem
            ),
            formulation = formulation_record(case),
            port_order = case.port_order,
            frequencies = copy(case.problem.frequencies),
            reference_execution = execution,
            reference,
            accepted = candidate,
            reference_comparison,
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

function _execution_record(execution::NamedTuple)
    required = (:backend, :elapsed_seconds, :exit_code)
    missing = filter(field -> !hasproperty(execution, field), required)
    isempty(missing) || throw(ArgumentError(
        "reference execution is missing fields: $(join(missing, ", "))",
    ))
    version = if hasproperty(execution, :version)
        execution.version
    elseif hasproperty(execution, :pscad_version)
        execution.pscad_version
    else
        throw(ArgumentError("reference execution is missing its backend version"))
    end
    return (
        backend = execution.backend,
        version = String(version),
        elapsed_seconds = execution.elapsed_seconds,
        exit_code = execution.exit_code
    )
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
        reference_comparison::LineParametersBenchmark,
        julia_benchmark;
        reference_execution::NamedTuple,
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    destination = case_stage(case.backend, case.name; artifact_root)
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
        reference_execution
    )
    digest = _snapshot_digest(snapshot)
    write(joinpath(temporary, "snapshot.sha256"), "$digest  snapshot.jld2\n")
    _atomic_stage!(temporary, destination)
    return (
        path = destination,
        snapshot_sha256 = digest,
        backend = case.backend,
        gauntlet_version = GAUNTLET_VERSION
    )
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
        "gauntlet_version", "case_name", "backend", "case_sha256", "problem",
        "formulation", "port_order", "frequencies", "reference_execution",
        "reference", "accepted", "reference_comparison", "julia_benchmark",
        "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["gauntlet_version"] == string(GAUNTLET_VERSION) || throw(ArgumentError(
        "Gauntlet snapshot $path uses version $(snapshot["gauntlet_version"]), " *
        "expected $(GAUNTLET_VERSION)",
    ))
    snapshot["case_name"] == string(case.name) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to case $(snapshot["case_name"]), not $(case.name)",
    ))
    snapshot["backend"] == string(case.backend) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to backend $(snapshot["backend"]), " *
        "not $(case.backend)",
    ))
    snapshot["case_sha256"] == case_digest(case) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match $(case.source_file). " *
        "Review the case and record it explicitly.",
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
        artifacts_toml::AbstractString = ARTIFACTS_TOML,
        force::Bool = gauntlet_force()
)
    force && return nothing
    hash = _bound_hash(case; artifacts_toml)
    hash === nothing && return nothing
    if !artifact_exists(hash)
        artifact = artifact_name(case.backend)
        try
            ensure_artifact_installed(artifact, artifacts_toml)
        catch error
            throw(ErrorException(
                "Gauntlet collection $artifact is not installed and has no usable " *
                "download. Original error: $(sprint(showerror, error))",
            ))
        end
    end
    path = joinpath(
        artifact_path(hash),
        "cases",
        string(case.name),
        "snapshot.jld2"
    )
    isfile(path) || return nothing
    return _load_snapshot_path(case, path)
end

function run_snapshot(
        case::GauntletCase;
        path::Union{Nothing, AbstractString} = nothing,
        artifacts_toml::AbstractString = ARTIFACTS_TOML,
        options::ComputeOptions = ComputeOptions(),
        benchmark_samples::Int = 10,
        benchmark_seconds::Real = 10
)
    loaded = load_snapshot(case; path, artifacts_toml)
    candidate = compute!(case.problem, case.formulation; options)
    validate_structure(case, candidate)
    reference_comparison = compare(loaded.reference, candidate)
    regression_comparison = compare(loaded.accepted, candidate)
    timing = benchmark_local(
        case;
        options,
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
        reference_execution = loaded.metadata["reference_execution"],
        timings = (
            reference = loaded.metadata["reference_execution"].elapsed_seconds,
            julia = timing
        )
    )
end
