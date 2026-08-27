struct BenchmarkCalculation{P, F, O <: NamedTuple}
    id::Symbol
    owner::Symbol
    problem::P
    formulation::F
    options::O

    function BenchmarkCalculation(
            id::Symbol,
            owner::Symbol,
            problem::P,
            formulation::F,
            options::O
    ) where {P, F, O <: NamedTuple}
        occursin(_CASE_IDENTIFIER, string(id)) || throw(ArgumentError(
            "calculation identifiers must be lowercase; got $(repr(id))",
        ))
        owner in (:engine, :uq, :external) || throw(ArgumentError(
            "calculation owner must be :engine, :uq, or :external",
        ))
        return new{P, F, O}(id, owner, problem, formulation, options)
    end
end

function benchmark_calculation(
        id::Symbol,
        owner::Symbol,
        problem,
        formulation;
        options::NamedTuple = (;)
)
    return BenchmarkCalculation(id, owner, problem, formulation, options)
end

struct LineParametersPolicy end
struct UQMomentPolicy end

struct OwnedBenchmark{M <: LoadedCase, R <: BenchmarkCalculation,
                      C <: BenchmarkCalculation, P, T}
    id::Symbol
    case_id::Symbol
    collection::Symbol
    source_file::String
    source_sha256::String
    model::M
    reference::R
    candidate::C
    comparison_policy::P
    tolerances::T

    function OwnedBenchmark(
            id::Symbol,
            case_id::Symbol,
            collection::Symbol,
            source_file::String,
            source_sha256::String,
            model::M,
            reference::R,
            candidate::C,
            comparison_policy::P,
            tolerances::T
    ) where {M <: LoadedCase, R <: BenchmarkCalculation,
             C <: BenchmarkCalculation, P, T}
        occursin(_CASE_IDENTIFIER, string(id)) || throw(ArgumentError(
            "benchmark identifiers must be lowercase; got $(repr(id))",
        ))
        occursin(_CASE_IDENTIFIER, string(collection)) || throw(ArgumentError(
            "benchmark collections must be lowercase identifiers",
        ))
        case_id === model.id || throw(ArgumentError(
            "benchmark :$id names case :$case_id but loaded :$(model.id)",
        ))
        reference.id == candidate.id && throw(ArgumentError(
            "benchmark calculations must have distinct identifiers",
        ))
        reference.owner === :external && throw(ArgumentError(
            "OwnedBenchmark cannot execute an external reference",
        ))
        candidate.owner === :external && throw(ArgumentError(
            "OwnedBenchmark cannot execute an external candidate",
        ))
        isfile(source_file) || throw(ArgumentError(
            "benchmark source file is missing: $source_file",
        ))
        return new{M, R, C, P, T}(
            id,
            case_id,
            collection,
            source_file,
            source_sha256,
            model,
            reference,
            candidate,
            comparison_policy,
            tolerances
        )
    end
end

function benchmark_definition(
        id::Symbol,
        case_id::Symbol,
        collection::Symbol,
        source_file::AbstractString,
        model::LoadedCase,
        reference::BenchmarkCalculation,
        candidate::BenchmarkCalculation,
        comparison_policy,
        tolerances
)
    path = realpath(source_file)
    return OwnedBenchmark(
        id,
        case_id,
        collection,
        path,
        bytes2hex(sha256(read(path))),
        model,
        reference,
        candidate,
        comparison_policy,
        tolerances
    )
end

function _compute_owned(calculation::BenchmarkCalculation)
    return isempty(calculation.options) ?
           compute(calculation.problem, calculation.formulation) :
           compute(
               calculation.problem,
               calculation.formulation;
               options = calculation.options
           )
end

function _execute(calculation::BenchmarkCalculation)
    started = time_ns()
    result = _compute_owned(calculation)
    elapsed = (time_ns() - started) * 1.0e-9
    return (; result, elapsed_seconds = elapsed)
end

function _owned_performance_settings(tolerances)
    hasproperty(tolerances, :performance) || return nothing
    settings = tolerances.performance
    keys(settings) == (:minimum_speedup, :samples, :seconds) || throw(ArgumentError(
        "owned performance settings must contain minimum_speedup, samples, and seconds",
    ))
    settings.minimum_speedup isa Real && isfinite(settings.minimum_speedup) &&
    settings.minimum_speedup > 1 || throw(ArgumentError(
        "owned minimum speedup must be finite and greater than one",
    ))
    settings.samples isa Integer && !(settings.samples isa Bool) &&
    settings.samples > 0 || throw(ArgumentError(
        "owned timing samples must be a positive integer",
    ))
    settings.seconds isa Real && isfinite(settings.seconds) &&
    settings.seconds > 0 || throw(ArgumentError(
        "owned timing duration must be positive and finite",
    ))
    return (
        minimum_speedup = Float64(settings.minimum_speedup),
        samples = Int(settings.samples),
        seconds = Float64(settings.seconds)
    )
end

function _benchmark_owned(calculation::BenchmarkCalculation, settings)
    benchmark = BenchmarkTools.@benchmarkable _compute_owned($calculation) samples=settings.samples seconds=settings.seconds evals=1
    trial = BenchmarkTools.run(benchmark)
    times = Float64.(trial.times) .* 1.0e-9
    return (
        minimum_seconds = minimum(times),
        median_seconds = median(times),
        bytes = trial.memory,
        allocations = trial.allocs,
        samples = length(times)
    )
end

function _owned_performance(benchmark::OwnedBenchmark)
    settings = _owned_performance_settings(benchmark.tolerances)
    settings === nothing && return nothing
    reference = _benchmark_owned(benchmark.reference, settings)
    candidate = _benchmark_owned(benchmark.candidate, settings)
    speedup = candidate.median_seconds / reference.median_seconds
    comparable = !gauntlet_instrumented()
    passes = comparable ? speedup >= settings.minimum_speedup : nothing
    return (; reference, candidate, speedup, comparable, passes, settings)
end

_normalize(::LineParametersPolicy, result, ::LoadedCase) = result
_normalize(::UQMomentPolicy, result, model::LoadedCase) =
    extract_moments(result, model.port_order)

_comparison_policy_record(::LineParametersPolicy) = :line_parameters
_comparison_policy_record(::UQMomentPolicy) = :uq_moments

_owned_formulation_record(formulation::AnalyticalFormulation) =
    formulation_record(formulation)
function _owned_formulation_record(formulation::LineCableModels.LinearError)
    return (
        kind = :linear_error,
        inner = _owned_formulation_record(formulation.inner),
        options = formulation.options
    )
end
function _owned_formulation_record(formulation::LineCableModels.MonteCarlo)
    distribution = formulation.distribution isa Symbol ?
                   formulation.distribution : string(typeof(formulation.distribution))
    return (
        kind = :monte_carlo,
        inner = _owned_formulation_record(formulation.inner),
        trials = formulation.trials,
        confidence = formulation.confidence,
        cdf_tolerance = formulation.cdf_tol,
        distribution,
        seed = formulation.seed,
        return_samples = formulation.return_samples,
        return_histograms = formulation.return_histograms,
        bins = formulation.bins,
        options = formulation.options
    )
end

function calculation_record(calculation::BenchmarkCalculation)
    return (
        id = calculation.id,
        owner = calculation.owner,
        formulation = _owned_formulation_record(calculation.formulation),
        execution_options = calculation.options
    )
end

function _semantic_owned_formulation_record(record::NamedTuple)
    kind = get(record, :kind, nothing)
    kind === :linear_error && return merge(
        record,
        (inner = _semantic_owned_formulation_record(record.inner),)
    )
    kind === :monte_carlo || return record
    options = record.options
    normalized_options = (
        retain_details = get(options, :retain_details, false),
        on_error = get(options, :on_error, :fail),
        max_failures = get(options, :max_failures, 100)
    )
    return merge(
        record,
        (
            inner = _semantic_owned_formulation_record(record.inner),
            options = normalized_options
        )
    )
end

function _semantic_calculation_records(records::NamedTuple)
    return map(records) do calculation
        merge(
            calculation,
            (
                formulation = _semantic_owned_formulation_record(
                    calculation.formulation
                ),
            )
        )
    end
end

function _owned_metadata(benchmark::OwnedBenchmark, reference_result, candidate_result)
    monte_carlo = candidate_result isa LineCableModels.MonteCarloResult ? (
        trials = only(candidate_result.trial_counts),
        seed = candidate_result.root_seed,
        distribution = benchmark.candidate.formulation.distribution
    ) : nothing
    return (
        benchmark_id = benchmark.id,
        case_id = benchmark.case_id,
        collection = benchmark.collection,
        case_source_sha256 = benchmark.model.source_sha256,
        benchmark_source_sha256 = benchmark.source_sha256,
        parameter_manifest = parameter_manifest(benchmark.model),
        applied_variation = variation_record(benchmark.model.variation),
        correlation = correlation_record(benchmark.model),
        calculations = (
            reference = calculation_record(benchmark.reference),
            candidate = calculation_record(benchmark.candidate)
        ),
        comparison_policy = _comparison_policy_record(benchmark.comparison_policy),
        monte_carlo
    )
end

function _moment_record(value::MomentResult)
    _validate_moment_result(value)
    return (
        values = value.values,
        frequencies = copy(value.frequencies),
        basis = value.basis,
        domain = nameof(value.domain),
        port_order = copy(value.port_order)
    )
end

function _moment_from_record(record)
    record.domain === :PhaseDomain || throw(ArgumentError(
        "only phase-domain UQ moment snapshots are supported",
    ))
    return MomentResult(
        record.values,
        record.frequencies,
        record.basis,
        PhaseDomain,
        record.port_order
    )
end

_moment_comparison_record(value::MomentBenchmark) = value.errors

function _same_moment_errors(observed::NamedTuple, stored::NamedTuple)
    keys(observed) == keys(stored) || return false
    return all(keys(observed)) do quantity
        left_product = getproperty(observed, quantity)
        right_product = getproperty(stored, quantity)
        all((:mean, :std)) do statistic
            left = getproperty(left_product, statistic)
            right = getproperty(right_product, statistic)
            isequal(left.absolute, right.absolute) &&
                isequal(left.relative, right.relative)
        end
    end
end

function _owned_comparison_passes(
        ::LineParametersPolicy,
        comparison::LineParametersBenchmark,
        tolerances
)
    return comparison_passes(comparison.Z, tolerances.Z) &&
           comparison_passes(comparison.Y, tolerances.Y)
end

function _owned_snapshot_path(
        benchmark::OwnedBenchmark;
        artifacts_toml::AbstractString = ARTIFACTS_TOML,
        required::Bool = true
)
    artifact = artifact_name(benchmark.collection)
    isfile(artifacts_toml) || (required ? throw(ArgumentError(
        "Gauntlet Artifacts.toml is missing: $artifacts_toml",
    )) : return nothing)
    hash = artifact_hash(artifact, artifacts_toml)
    hash === nothing && (required ? throw(ArgumentError(
        "Gauntlet collection $artifact is not bound in $artifacts_toml",
    )) : return nothing)
    if !artifact_exists(hash)
        try
            ensure_artifact_installed(artifact, artifacts_toml)
        catch error
            required || return nothing
            throw(ErrorException(
                "Gauntlet collection $artifact is unavailable. Original error: " *
                sprint(showerror, error),
            ))
        end
    end
    path = joinpath(
        artifact_path(hash),
        "benchmarks",
        string(benchmark.id),
        "snapshot.jld2"
    )
    isfile(path) || (required ? throw(ArgumentError(
        "Gauntlet collection $artifact has no snapshot for $(benchmark.id)",
    )) : return nothing)
    return path
end

snapshot_path(benchmark::OwnedBenchmark; artifacts_toml::AbstractString = ARTIFACTS_TOML) =
    _owned_snapshot_path(benchmark; artifacts_toml, required = true)

function _write_owned_snapshot(
        path::AbstractString,
        benchmark::OwnedBenchmark,
        reference::MomentResult,
        candidate::MomentResult,
        comparison::MomentBenchmark,
        timings,
        run_metadata
)
    mkpath(dirname(path))
    temporary = path * ".new"
    isfile(temporary) && rm(temporary; force = true)
    try
        JLD2.jldsave(
            temporary;
            schema_version = SNAPSHOT_SCHEMA_VERSION,
            benchmark_id = string(benchmark.id),
            case_id = string(benchmark.case_id),
            collection = string(benchmark.collection),
            case_source_sha256 = run_metadata.case_source_sha256,
            benchmark_source_sha256 = run_metadata.benchmark_source_sha256,
            parameter_manifest = run_metadata.parameter_manifest,
            applied_variation = run_metadata.applied_variation,
            correlation = run_metadata.correlation,
            calculations = run_metadata.calculations,
            comparison_policy = run_metadata.comparison_policy,
            tolerances = benchmark.tolerances,
            monte_carlo = run_metadata.monte_carlo,
            port_order = benchmark.model.port_order,
            frequencies = copy(reference.frequencies),
            reference = _moment_record(reference),
            accepted_reference = _moment_record(reference),
            accepted_candidate = _moment_record(candidate),
            reference_comparison = _moment_comparison_record(comparison),
            timings,
            environment = _performance_identity(),
            recorded_at_utc = string(now(UTC))
        )
        mv(temporary, path; force = true)
    catch
        isfile(temporary) && rm(temporary; force = true)
        rethrow()
    end
    return path
end

function persist_snapshot(
        benchmark::OwnedBenchmark,
        reference::MomentResult,
        candidate::MomentResult,
        comparison::MomentBenchmark,
        timings,
        run_metadata;
        artifact_root::AbstractString = ARTIFACT_ROOT
)
    destination = benchmark_stage(
        benchmark.collection,
        benchmark.id;
        artifact_root
    )
    temporary = destination * ".new"
    ispath(temporary) && rm(temporary; recursive = true, force = true)
    mkpath(temporary)
    snapshot = joinpath(temporary, "snapshot.jld2")
    _write_owned_snapshot(
        snapshot,
        benchmark,
        reference,
        candidate,
        comparison,
        timings,
        run_metadata
    )
    digest = _snapshot_digest(snapshot)
    write(joinpath(temporary, "snapshot.sha256"), "$digest  snapshot.jld2\n")
    _atomic_stage!(temporary, destination)
    return (
        path = destination,
        snapshot_sha256 = digest,
        benchmark_id = benchmark.id,
        case_id = benchmark.case_id,
        collection = benchmark.collection,
        schema_version = SNAPSHOT_SCHEMA_VERSION
    )
end

function _load_owned_snapshot_path(
        benchmark::OwnedBenchmark,
        path::AbstractString
)
    isfile(path) || throw(ArgumentError("Gauntlet snapshot is missing: $path"))
    digest_path = joinpath(dirname(path), "snapshot.sha256")
    isfile(digest_path) || throw(ArgumentError(
        "Gauntlet snapshot digest is missing: $digest_path",
    ))
    expected_digest = first(split(strip(read(digest_path, String))))
    expected_digest == _snapshot_digest(path) || throw(ArgumentError(
        "Gauntlet snapshot SHA-256 does not match $digest_path",
    ))
    snapshot = JLD2.load(path)
    required = (
        "schema_version", "benchmark_id", "case_id",
        "collection", "case_source_sha256", "benchmark_source_sha256",
        "parameter_manifest", "applied_variation", "correlation",
        "calculations", "comparison_policy", "tolerances", "monte_carlo", "port_order",
        "frequencies", "accepted_reference", "accepted_candidate",
        "reference_comparison", "timings", "environment", "recorded_at_utc"
    )
    missing = filter(key -> !haskey(snapshot, key), required)
    isempty(missing) || throw(ArgumentError(
        "Gauntlet snapshot $path is missing fields: $(join(missing, ", "))",
    ))
    snapshot["schema_version"] == SNAPSHOT_SCHEMA_VERSION || throw(ArgumentError(
        "Gauntlet snapshot $path does not use schema $SNAPSHOT_SCHEMA_VERSION",
    ))
    snapshot["benchmark_id"] == string(benchmark.id) ||
        throw(ArgumentError("Gauntlet snapshot benchmark ID does not match"))
    snapshot["case_id"] == string(benchmark.case_id) ||
        throw(ArgumentError("Gauntlet snapshot case ID does not match"))
    snapshot["collection"] == string(benchmark.collection) ||
        throw(ArgumentError("Gauntlet snapshot collection does not match"))
    current = _owned_metadata(benchmark, nothing, nothing)
    snapshot["case_source_sha256"] == bytes2hex(sha256(read(benchmark.model.source_file))) ||
        throw(ArgumentError("Gauntlet snapshot case source digest does not match"))
    snapshot["benchmark_source_sha256"] == bytes2hex(sha256(read(benchmark.source_file))) ||
        throw(ArgumentError("Gauntlet snapshot benchmark source digest does not match"))
    for field in ("parameter_manifest", "applied_variation", "correlation",
                  "comparison_policy")
        isequal(snapshot[field], getproperty(current, Symbol(field))) ||
            throw(ArgumentError("Gauntlet snapshot $field does not match"))
    end
    isequal(
        _semantic_calculation_records(snapshot["calculations"]),
        _semantic_calculation_records(current.calculations)
    ) || throw(ArgumentError("Gauntlet snapshot calculations does not match"))
    isequal(snapshot["tolerances"], benchmark.tolerances) ||
        throw(ArgumentError("Gauntlet snapshot tolerances do not match"))
    snapshot["port_order"] == benchmark.model.port_order ||
        throw(ArgumentError("Gauntlet snapshot terminal order does not match"))
    accepted_reference = _moment_from_record(snapshot["accepted_reference"])
    accepted_candidate = _moment_from_record(snapshot["accepted_candidate"])
    snapshot["frequencies"] == accepted_reference.frequencies ||
        throw(ArgumentError("Gauntlet snapshot frequency record does not match"))
    observed = compare(accepted_reference, accepted_candidate)
    _same_moment_errors(observed.errors, snapshot["reference_comparison"]) ||
        throw(ArgumentError("Gauntlet snapshot comparison does not match its moments"))
    return (; accepted_reference, accepted_candidate, metadata = snapshot)
end

function load_snapshot(
        benchmark::OwnedBenchmark;
        path::Union{Nothing, AbstractString} = nothing,
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    resolved = path === nothing ? snapshot_path(benchmark; artifacts_toml) : String(path)
    return _load_owned_snapshot_path(benchmark, resolved)
end

function load_prior_snapshot(
        benchmark::OwnedBenchmark;
        artifacts_toml::AbstractString = ARTIFACTS_TOML
)
    path = _owned_snapshot_path(benchmark; artifacts_toml, required = false)
    path === nothing && return nothing
    return _load_owned_snapshot_path(benchmark, path)
end

function _owned_comparison_passes(
        ::UQMomentPolicy,
        comparison::MomentBenchmark,
        tolerances
)
    return moment_comparison_passes(comparison, tolerances)
end

function run_benchmark(benchmark::OwnedBenchmark)
    mode = gauntlet_mode()
    benchmark.source_sha256 == bytes2hex(sha256(read(benchmark.source_file))) ||
        throw(ArgumentError("benchmark source changed after definition"))
    benchmark.model.source_sha256 == bytes2hex(sha256(read(benchmark.model.source_file))) ||
        throw(ArgumentError("case source changed after loading"))
    if mode === :record
        get(ENV, "LINECABLEMODELS_GAUNTLET_PERSIST", "false") == "true" ||
            throw(ArgumentError(
                "record mode requires LINECABLEMODELS_GAUNTLET_PERSIST=true",
            ))
        get(ENV, "LINECABLEMODELS_GAUNTLET_RUNNER", "false") == "true" ||
            throw(ArgumentError(
                "record mode must run through test/gauntlet/runtests.jl",
            ))
    end
    accepted = mode === :snapshot ? load_snapshot(benchmark) :
               mode === :record ? load_prior_snapshot(benchmark) : nothing
    reference_execution = _execute(benchmark.reference)
    candidate_execution = _execute(benchmark.candidate)
    reference = _normalize(
        benchmark.comparison_policy,
        reference_execution.result,
        benchmark.model
    )
    candidate = _normalize(
        benchmark.comparison_policy,
        candidate_execution.result,
        benchmark.model
    )
    comparison = compare(reference, candidate)
    passes = _owned_comparison_passes(
        benchmark.comparison_policy,
        comparison,
        benchmark.tolerances.reference
    )
    failure_detail = if benchmark.comparison_policy isa UQMomentPolicy
        summary = moment_error_summary(
            comparison,
            benchmark.tolerances.reference
        )
        "; errors=$(repr(summary))"
    else
        ""
    end
    passes || throw(ArgumentError(
        "owned benchmark $(benchmark.id) exceeds its reference tolerance" *
        failure_detail,
    ))
    regression = accepted === nothing ? nothing : (
        reference = compare(accepted.accepted_reference, reference),
        candidate = compare(accepted.accepted_candidate, candidate)
    )
    if regression !== nothing
        _owned_comparison_passes(
            benchmark.comparison_policy,
            regression.reference,
            benchmark.tolerances.regression
        ) || throw(ArgumentError(
            "owned benchmark $(benchmark.id) reference regression exceeds tolerance",
        ))
        _owned_comparison_passes(
            benchmark.comparison_policy,
            regression.candidate,
            benchmark.tolerances.regression
        ) || throw(ArgumentError(
            "owned benchmark $(benchmark.id) candidate regression exceeds tolerance",
        ))
    end
    performance = _owned_performance(benchmark)
    performance === nothing || performance.passes === nothing ||
        performance.passes || throw(ArgumentError(
            "owned benchmark $(benchmark.id) Monte Carlo/LEP speedup " *
            "$(performance.speedup) is below " *
            "$(performance.settings.minimum_speedup)",
        ))
    timings = (
        execution = (
            reference = reference_execution.elapsed_seconds,
            candidate = candidate_execution.elapsed_seconds
        ),
        benchmark = performance
    )
    metadata = _owned_metadata(
        benchmark,
        reference_execution.result,
        candidate_execution.result
    )
    artifact = mode === :record ? persist_snapshot(
        benchmark,
        reference,
        candidate,
        comparison,
        timings,
        metadata
    ) : nothing
    return (
        mode,
        reference_result = reference_execution.result,
        candidate_result = candidate_execution.result,
        reference,
        candidate,
        comparison,
        passes,
        regression,
        artifact,
        timings,
        performance,
        metadata
    )
end
