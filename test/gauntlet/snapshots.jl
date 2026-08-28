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
        "Record, package, and publish the collection explicitly before running snapshot mode.",
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
        "benchmarks",
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
            schema_version = SNAPSHOT_SCHEMA_VERSION,
            benchmark_id = string(case.name),
            case_id = string(case.case_id),
            collection = string(case.backend),
            case_name = string(case.name),
            backend = string(case.backend),
            case_source_sha256 = case_digest(case),
            benchmark_source_sha256 = benchmark_digest(case),
            parameter_manifest = case.parameter_manifest,
            applied_variation = case.applied_variation,
            correlation = case.correlation,
            problem = (
                reference = case.reference_problem,
                candidate = case.problem
            ),
            formulation = formulation_record(case),
            calculations = (
                reference = (
                    id = case.backend,
                    owner = :external,
                    formulation = formulation_record(case).reference,
                    options = (;)
                ),
                candidate = (
                    id = :line_cable_models,
                    owner = :engine,
                    formulation = formulation_record(case).candidate,
                    options = (;)
                )
            ),
            comparison_policy = :line_parameters,
            tolerances = case.tolerances,
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
    destination = benchmark_stage(case.backend, case.name; artifact_root)
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
        benchmark_id = case.name,
        case_id = case.case_id,
        collection = case.backend,
        backend = case.backend,
        schema_version = SNAPSHOT_SCHEMA_VERSION
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
        "schema_version", "benchmark_id", "case_id",
        "collection", "case_source_sha256", "benchmark_source_sha256",
        "parameter_manifest", "applied_variation", "correlation", "problem",
        "formulation", "calculations", "comparison_policy", "tolerances", "port_order",
        "frequencies", "reference_execution",
        "reference", "accepted", "reference_comparison", "julia_benchmark",
        "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["schema_version"] == SNAPSHOT_SCHEMA_VERSION || throw(ArgumentError(
        "Gauntlet snapshot $path uses schema $(snapshot["schema_version"]), " *
        "expected $SNAPSHOT_SCHEMA_VERSION",
    ))
    snapshot["benchmark_id"] == string(case.name) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to benchmark $(snapshot["benchmark_id"]), not $(case.name)",
    ))
    snapshot["case_id"] == string(case.case_id) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to case $(snapshot["case_id"]), not $(case.case_id)",
    ))
    snapshot["collection"] == string(case.backend) || throw(ArgumentError(
        "Gauntlet snapshot $path belongs to collection $(snapshot["collection"]), " *
        "not $(case.backend)",
    ))
    snapshot["case_source_sha256"] == case_digest(case) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match the physical case source.",
    ))
    snapshot["benchmark_source_sha256"] == benchmark_digest(case) || throw(ArgumentError(
        "Gauntlet snapshot $path does not match the benchmark source.",
    ))
    isequal(snapshot["parameter_manifest"], case.parameter_manifest) ||
        throw(ArgumentError("Gauntlet snapshot parameter manifest does not match"))
    isequal(snapshot["applied_variation"], case.applied_variation) ||
        throw(ArgumentError("Gauntlet snapshot variation does not match"))
    isequal(snapshot["correlation"], case.correlation) ||
        throw(ArgumentError("Gauntlet snapshot correlation record does not match"))
    isequal(snapshot["tolerances"], case.tolerances) ||
        throw(ArgumentError("Gauntlet snapshot tolerances do not match"))
    _semantic_formulation_pair(snapshot["formulation"]) ==
    _semantic_formulation_pair(formulation_record(case)) || throw(ArgumentError(
        "Gauntlet snapshot formulation does not match the case definition",
    ))
    stored_problem = snapshot["problem"]
    _same_problem_structure(stored_problem.reference, case.reference_problem) ||
        throw(ArgumentError(
            "Gauntlet snapshot reference problem structure does not match the case definition",
        ))
    _same_problem_structure(stored_problem.candidate, case.problem) ||
        throw(ArgumentError(
            "Gauntlet snapshot candidate problem structure does not match the case definition",
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
        "benchmarks",
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
        options::NamedTuple = (;)
)
    run = computation_options(Val(GauntletCase), options)
    candidate_options = merge(
        run.candidate,
        (output_basis = run.output_basis,)
    )
    loaded = load_snapshot(case; path, artifacts_toml)
    candidate = compute(
        case.problem,
        case.formulation;
        options = candidate_options
    )
    validate_structure(case, candidate)
    reference_comparison = compare(loaded.reference, candidate)
    regression_comparison = compare(loaded.accepted, candidate)
    timing = benchmark_local(
        case;
        options = candidate_options,
        samples = run.benchmark.samples,
        seconds = run.benchmark.seconds
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
