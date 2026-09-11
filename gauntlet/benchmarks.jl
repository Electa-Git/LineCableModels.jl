function _benchmark_performance_settings(tolerances)
    haskey(tolerances, :performance) || return nothing
    settings = tolerances.performance
    keys(settings) == (:minimum_speedup, :samples, :seconds) || throw(ArgumentError(
        "benchmark performance settings must contain minimum_speedup, samples, and seconds",
    ))
    settings.minimum_speedup isa Real && isfinite(settings.minimum_speedup) &&
    settings.minimum_speedup > 1 || throw(ArgumentError(
        "benchmark minimum speedup must be finite and greater than one",
    ))
    settings.samples isa Integer && !(settings.samples isa Bool) &&
    settings.samples > 0 || throw(ArgumentError(
        "benchmark timing samples must be a positive integer",
    ))
    settings.seconds isa Real && isfinite(settings.seconds) &&
    settings.seconds > 0 || throw(ArgumentError(
        "benchmark timing duration must be positive and finite",
    ))
    return (
        minimum_speedup = Float64(settings.minimum_speedup),
        samples = Int(settings.samples),
        seconds = Float64(settings.seconds)
    )
end

function _benchmark_owned(calculation::BenchmarkCalculation, settings)
    observations=NamedTuple[]
    started=time_ns()
    for _ in 1:settings.samples
        elapsed=@timed _execute(calculation)
        execution=get(details(elapsed.value.result), :execution, (;))
        reused=get(execution, :reused, false)
        push!(observations, (seconds = elapsed.time, bytes = elapsed.bytes, reused))
        (time_ns()-started)*1e-9 >= settings.seconds && break
    end
    return (
        scope = :compute_wall, median_seconds = median(row.seconds for row in observations),
        bytes = maximum(row.bytes for row in observations), samples = length(observations),
        observations, environment = _performance_identity())
end

function _benchmark_performance(benchmark::BenchmarkDefinition)
    settings=_benchmark_performance_settings(benchmark.tolerances)
    settings === nothing && return nothing
    reference=_benchmark_owned(benchmark.reference, settings)
    candidate=_benchmark_owned(benchmark.candidate, settings)
    # The candidate is the implementation under test, so its speedup over the
    # reference is the reference wall time divided by the candidate wall time.
    speedup=reference.median_seconds/candidate.median_seconds
    comparable=!gauntlet_instrumented() &&
               !any(row.reused for row in reference.observations) &&
               !any(row.reused for row in candidate.observations)
    passes=comparable ? speedup >= settings.minimum_speedup : nothing
    return (; reference, candidate, speedup, comparable, passes, settings)
end

_normalize(result::Union{AbstractCoreResult, AbstractParametricResult, MomentResult}, model) = result
function _normalize(result::AbstractUncertaintyResult, model)
    extract_moments(result, model.port_order)
end

"""
    validate(benchmark::BenchmarkDefinition)

Check comparison controls, declared terminal order and timing settings before
either calculation starts. Checks requiring calculated values remain in `compare`.
"""
function validate(benchmark::BenchmarkDefinition)
    settings = benchmark.comparison_settings
    validate(BenchmarkTableDefinition(; settings...))
    _benchmark_performance_settings(benchmark.tolerances)
    if haskey(benchmark.tolerances, :reference)
        settings.statistics == (:mean, :std) || throw(ArgumentError(
            "reference acceptance limits apply only to declared moment comparisons"))
        limits=benchmark.tolerances.reference
        keys(limits) == (:mean, :std) || throw(ArgumentError("moment tolerances need mean and std"))
        for group in limits
            keys(group) == (:R, :L, :C, :G) || throw(ArgumentError("moment tolerances need R, L, C, G"))
            for limit in group
                keys(limit) == (:absolute, :relative) && all(v -> v isa Real && isfinite(v) && v >= 0, limit) ||
                    throw(ArgumentError("moment limits must be finite nonnegative absolute and relative tolerances"))
            end
        end
    end
    if benchmark.reference.formulation isa Union{Gridspace,LineCableModels.Combinatorial}
        settings.pairing === nothing && throw(ArgumentError("a reference result space requires explicit pairing before execution"))
    end
    a, b = benchmark.reference.problem, benchmark.candidate.problem
    if a isa LineParametersProblem && b isa LineParametersProblem
        a.system.terminal_order == b.system.terminal_order &&
            a.system.connection_order == b.system.connection_order || throw(ArgumentError(
                "benchmark terminal identities or ordering differ"))
        a.frequencies == b.frequencies || throw(ArgumentError("benchmark frequency coordinates differ"))
    end
    return nothing
end

function formulation_record(formulation::PSCAD.PSCADFormulation)
    _selection_value(NamedTuple(formulation))
end
function formulation_record(formulation::AbstractFormulation)
    _selection_value(formulation)
end
function formulation_record(formulation::Engine.LineCableModelsFEM)
    return _selection_value(NamedTuple(formulation))
end
function formulation_record(formulation::LineCableModels.LinearError)
    return (
        kind = :linear_error,
        inner = formulation_record(formulation.inner),
        options = formulation.options
    )
end
function formulation_record(formulation::LineCableModels.MonteCarlo)
    distribution = formulation.distribution isa Symbol ?
                   formulation.distribution : string(typeof(formulation.distribution))
    return (
        kind = :monte_carlo,
        inner = formulation_record(formulation.inner),
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


"""
    run_benchmark(definition; directory=nothing, mode=:live)

Compute the two declared operands with their unchanged problems, formulations and
options. Comparison direction and RMS settings belong to the definition. Differences
between models are retained observations. Optional performance checks are separate.
When `directory` is supplied, completed calculations and analysis are recoverable.
Reuse depends on the numerical declaration and saved-file integrity, not live source
files. Execution-session metadata records provenance without freezing the environment.
Reports are saved before optional timing checks, so timing failures do not discard them.
`mode=:record` stages that same complete bundle for explicit artifact packaging.
"""
function run_benchmark(benchmark::BenchmarkDefinition; directory = nothing,
        implementation = (), session = nothing, mode::Symbol = :live,
        measure_performance::Bool = true, recover_solvers::Bool = false)
    validate(benchmark)
    directory === nothing || validate(Base.write,directory)
    calculations=(reference=calculation_record(benchmark.reference),
        candidate=calculation_record(benchmark.candidate))
    mode in (:live, :record) || throw(ArgumentError(
        "run_benchmark executes declarations; use read_collection for retained snapshots"))
    if directory === nothing && mode === :record
        directory=benchmark_stage(benchmark.collection, benchmark.id)
        ispath(directory) &&
            throw(ArgumentError("staged benchmark already exists: $directory"))
    end
    session === nothing && directory !== nothing && (session=execution_record())
    reference_execution=_execute(benchmark.reference;
        directory = directory === nothing ? nothing : joinpath(directory, "reference"),
        model = benchmark.model, implementation, session, recover_solvers)
    candidate_execution=_execute(benchmark.candidate;
        directory = directory === nothing ? nothing : joinpath(directory, "candidate"),
        model = benchmark.model, implementation, session, recover_solvers)
    # Result-space axes retain the actual resolved problems, including port identity.
    expected=benchmark.model.nominal_problem.system
    for (calculation, execution) in ((benchmark.reference, reference_execution),
        (benchmark.candidate, candidate_execution))
        problems=execution.result isa ParametricResult ? NamedTuple(execution.result).axes.problems :
                 (calculation.problem,)
        for problem in problems
            problem isa LineParametersProblem || continue
            problem.system.terminal_order == expected.terminal_order &&
            problem.system.connection_order == expected.connection_order ||
                throw(ArgumentError("benchmark result-space terminal identities or ordering differ"))
        end
    end
    reference=_normalize(reference_execution.result, benchmark.model)
    candidate=_normalize(candidate_execution.result, benchmark.model)
    definition=BenchmarkTableDefinition(;benchmark.comparison_settings...)
    reference_metadata=(port_order=get(details(reference_execution.result isa ParametricResult ? first(reference_execution.result) : reference_execution.result),:coordinates,benchmark.model.port_order),
        formulation=calculation_record(benchmark.reference).formulation, axes=reference isa ParametricResult ? reference.axes : nothing)
    candidate_metadata=(port_order=get(details(candidate_execution.result isa ParametricResult ? first(candidate_execution.result) : candidate_execution.result),:coordinates,benchmark.model.port_order),
        formulation=calculation_record(benchmark.candidate).formulation, axes=candidate isa ParametricResult ? candidate.axes : nothing)
    publication=report(definition,(
        reference=(result=reference,metadata=reference_metadata),
        candidate=(result=candidate,metadata=candidate_metadata),
        context=(id=benchmark.id,case_id=benchmark.case_id,collection=benchmark.collection)))
    comparison=publication.published.comparisons
    passes=haskey(benchmark.tolerances, :reference) ?
           moment_comparison_passes(comparison, benchmark.tolerances.reference) : nothing
    metadata=(
        benchmark_id = benchmark.id,
        case_id = benchmark.case_id,
        collection = benchmark.collection,
        case_source_sha256 = benchmark.model.source_sha256,
        benchmark_source_sha256 = benchmark.source_sha256,
        parameter_manifest = parameter_manifest(benchmark.model),
        applied_variation = variation_record(benchmark.model.variation),
        correlation = correlation_record(benchmark.model),
        session,
        calculations,
        comparison_settings = benchmark.comparison_settings
    )
    artifact=nothing
    if directory !== nothing
        operands=map((:reference, :candidate)) do role
            saved=read_calculation(joinpath(directory, string(role), "calculation.jld2"))
            BenchmarkCalculation(role, saved, saved.metadata.formulation)
        end
        retained=BenchmarkDefinition(
            benchmark.id, benchmark.case_id, benchmark.collection,
            benchmark.source_file, benchmark.source_sha256, benchmark.model, operands...,
            benchmark.comparison_settings, benchmark.tolerances)
        artifact=record_benchmark(retained, publication; directory = joinpath(directory, "analyses"))
    end
    performance_path=directory === nothing ? nothing : joinpath(directory,"performance.jld2")
    performance=if !measure_performance
        performance_path !== nothing && isfile(performance_path) ?
            JLD2.load(performance_path,"performance") : nothing
    else
        value=_benchmark_performance(benchmark)
        if performance_path !== nothing && value !== nothing
            temporary=tempname(directory)
            try
                JLD2.jldsave(temporary;performance=value,session)
                mv(temporary,performance_path;force=true)
            finally
                isfile(temporary) && rm(temporary)
            end
        end
        value
    end
    timings=(scope = :compute_wall,
        execution = (
            reference = (seconds = reference_execution.elapsed_seconds,
                reused = reference_execution.reused, session=reference_execution.session),
            candidate = (seconds = candidate_execution.elapsed_seconds,
                reused = candidate_execution.reused, session=candidate_execution.session)),
        benchmark = performance)
    return (; mode, reference_result = reference_execution.result,
        candidate_result = candidate_execution.result,
        reference, candidate, comparison,
        passes, performance, timings, metadata, artifact, report=publication)
end

"""
    benchmark_definition(model; id, source_file, reference, formulations, kwargs...)

Declare one case, one explicit reference and a scalar or Gridspace of candidate
formulations. The candidate follows the ordinary compute overload. ReportBuilder
owns quantities, bands and comparison controls.
"""
function benchmark_definition(model::LoadedCase; id::Symbol, source_file::AbstractString,
        reference, formulations, collection::Symbol=:manual,
        options::NamedTuple=(;), report=BenchmarkTableDefinition(), tolerances=(;))
    baseline=reference isa BenchmarkCalculation ? reference :
        BenchmarkCalculation(:reference,model.problem,reference)
    candidate=BenchmarkCalculation(:candidate,model.problem,formulations;options)
    return benchmark_definition(id,model.id,collection,source_file,model,baseline,candidate,
        report.settings,tolerances)
end
