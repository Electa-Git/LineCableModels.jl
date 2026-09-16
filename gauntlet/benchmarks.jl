function _benchmark_performance_settings(tolerances)
    haskey(tolerances, :performance) || return nothing
    settings = tolerances.performance
    # Saved work orders may carry a former speedup threshold. It has no role in
    # measurement or execution; retain only the requested sampling budget.
    all(key -> haskey(settings, key), (:samples, :seconds)) &&
    isempty(setdiff(keys(settings), (:samples, :seconds, :minimum_speedup))) || throw(ArgumentError(
        "benchmark performance settings must contain samples and seconds",
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
        samples = Int(settings.samples),
        seconds = Float64(settings.seconds)
    )
end

function _benchmark_owned(calculation::BenchmarkCalculation, settings; role=calculation.id)
    observations=NamedTuple[]
    prepared=_performance_calculation(calculation)
    external=_external_formulation(prepared.formulation)
    if !external
        _performance_span(;sample="warmup",samples=settings.samples,role,
                record=_-> _performance_observation!(role;sample="warmup")) do
            _compute_calculation(prepared)
        end
    end
    _performance_observation!(role;sample="warmup")
    started=time_ns()
    for sample in 1:settings.samples
        _performance_span(;sample,samples=settings.samples,role,record=elapsed->begin
            reused=_result_reused(elapsed.value)
            push!(observations,(seconds=elapsed.time,bytes=elapsed.bytes,reused,
                source_timings=_source_timings(elapsed.value)))
            _performance_observation!(role;sample,seconds=elapsed.time,reused)
        end) do
            @timed _compute_calculation(prepared)
        end
        (time_ns()-started)*1e-9 >= settings.seconds && break
    end
    _performance_observation!(role;finished=true)
    return (
        scope = :compute_call_wall, median_seconds = median(row.seconds for row in observations),
        bytes = maximum(row.bytes for row in observations), samples = length(observations),
        observations, environment = _performance_identity(),
        calculation = _numerical_record(calculation_record(prepared)),
        policy=(progress=false,diagnostics=:quiet,callbacks=false,
            warmup=external ? :native_not_repeated : :owned_call,
            allocation_scope=:julia, allocation_statistic=:maximum,settings=_selection_value(prepared.options)))
end

function _benchmark_performance(benchmark::BenchmarkDefinition)
    settings=_benchmark_performance_settings(benchmark.tolerances)
    settings === nothing && return nothing
    reference=_benchmark_owned(benchmark.reference, settings;role=:reference)
    candidate=_benchmark_owned(benchmark.candidate, settings;role=:candidate)
    # Keep the declared direction: reference wall time divided by candidate time.
    speedup=reference.median_seconds/candidate.median_seconds
    comparable=!gauntlet_instrumented() &&
               reference.scope === candidate.scope &&
               reference.environment == candidate.environment &&
               all(key->getproperty(reference.policy,key)==getproperty(candidate.policy,key),
                   (:progress,:diagnostics,:callbacks,:allocation_scope)) &&
               !any(row.reused for row in reference.observations) &&
               !any(row.reused for row in candidate.observations)
    return (; reference, candidate, speedup, comparable, settings)
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
        options = formulation.options.data
    )
end
function formulation_record(formulation::LineCableModels.MonteCarlo)
    distribution = formulation.options.data.distribution isa Symbol ?
                   formulation.options.data.distribution : string(typeof(formulation.options.data.distribution))
    return (
        kind = :monte_carlo,
        inner = formulation_record(formulation.inner),
        trials = formulation.options.data.trials,
        confidence = formulation.options.data.confidence,
        cdf_tolerance = formulation.options.data.cdf_tol,
        distribution,
        seed = formulation.options.data.seed,
        return_samples = formulation.options.data.return_samples,
        return_histograms = formulation.options.data.return_histograms,
        bins = formulation.options.data.bins,
        options = _selection_value(formulation.options)
    )
end


"""
    run_benchmark(definition; directory=nothing, mode=:live)

Compute the two declared operands with their unchanged problems, formulations and
options. Comparison direction and RMS settings belong to the definition. Differences
between models and optional timing ratios are retained observations, without
acceptance verdicts.
When `directory` is supplied, completed calculations and analysis are recoverable.
Reuse depends on the numerical declaration and saved-file integrity, not live source
files. Each execution session records the Julia and package versions and Git state.
A timing computation error still retains the comparison report before the error is propagated.
`mode=:record` stages that same complete bundle for explicit artifact packaging.
"""
function run_benchmark(benchmark::BenchmarkDefinition; directory = nothing,
        implementation = (), session = nothing, mode::Symbol = :live,
        measure_performance::Bool = true, recover_solvers::Bool = false,
        reuse_directories=String[])
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
    executions=map((:reference,:candidate)) do role
        receiver=LineCableModels.progress_receiver()
        receiver === nothing || LineCableModels.report_progress(receiver,(kind=:operand,role,state=:running,stage=:preparing,
            backend=_execution_backend(getproperty(benchmark,role).formulation),
            mode=_execution_mode(getproperty(benchmark,role).formulation)))
        execution=LineCableModels.with_progress_scope(;role) do
            _execute(getproperty(benchmark,role);
                directory=directory === nothing ? nothing : joinpath(directory,string(role)),
                model=benchmark.model,implementation,session,recover_solvers,
                reuse_directories=[joinpath(source,string(role)) for source in reuse_directories])
        end
        jobs_reused = execution.result isa ParametricResult ?
            count(value->_result_reused(value;partial=false),execution.result) :
            Int(_result_reused(execution.result;partial=false))
        jobs_reused = something(get(execution, :jobs_reused, nothing), jobs_reused)
        receiver === nothing || LineCableModels.report_progress(receiver,
            (kind=:operand,role,state=:complete,stage=:computed,
                reused=execution.reused,jobs_reused,
                seconds=execution.elapsed_seconds,
                compute_seconds=execution.reused ? nothing : get(execution.timing, :seconds, nothing),
                timing_reused=execution.reused || get(execution.timing, :reused_points, 0) > 0 ||
                    _result_reused(execution.result)))
        execution
    end
    reference_execution,candidate_execution=executions
    receiver=LineCableModels.progress_receiver()
    receiver === nothing || LineCableModels.report_progress(receiver,(stage=:validating,))
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
    reference=reference_execution.result
    candidate=candidate_execution.result
    definition=BenchmarkTableDefinition(;benchmark.comparison_settings...)
    reference_metadata=(port_order=get(details(reference_execution.result isa ParametricResult ? first(reference_execution.result) : reference_execution.result).data,:coordinates,benchmark.model.port_order),
        formulation=calculation_record(benchmark.reference).formulation, axes=reference isa ParametricResult ? reference.axes : nothing)
    candidate_metadata=(port_order=get(details(candidate_execution.result isa ParametricResult ? first(candidate_execution.result) : candidate_execution.result).data,:coordinates,benchmark.model.port_order),
        formulation=calculation_record(benchmark.candidate).formulation, axes=candidate isa ParametricResult ? candidate.axes : nothing)
    performance_path=directory === nothing ? nothing : joinpath(directory,"performance.jld2")
    performance_error=nothing
    performance_session=session
    performance_checks=(checksum_verified=missing,workload_verified=missing)
    performance=try
        if performance_path !== nothing && isfile(performance_path)
            retained=read_benchmark(performance_path,Val(:performance);calculations)
            if measure_performance && retained.performance !== nothing
                _benchmark_performance_settings((performance=retained.performance.settings,))==_benchmark_performance_settings(benchmark.tolerances) ||
                    throw(ArgumentError("retained performance settings differ; use a new benchmark attempt"))
            end
            performance_session=retained.session
            performance_checks=(checksum_verified=retained.checksum_verified,
                workload_verified=retained.workload_verified)
            retained.performance
        elseif !measure_performance
            nothing
    else
        value=nothing
        performance_session=session
        settings=_benchmark_performance_settings(benchmark.tolerances)
        if settings !== nothing && all(execution->execution.reused,executions)
            for source in reuse_directories
                saved_path=joinpath(source,"performance.jld2")
                isfile(saved_path) || continue
                saved=read_benchmark(saved_path,Val(:performance);
                    calculations=(reference=nothing,candidate=nothing))
                retained=saved.performance
                retained !== nothing && _benchmark_performance_settings((performance=retained.settings,))==settings || continue
                all(getproperty(retained,role).calculation==
                    _numerical_record(calculation_record(
                        _performance_calculation(getproperty(benchmark,role))))
                    for role in (:reference,:candidate)) || continue
                value=retained
                performance_session=saved.session
                performance_checks=(checksum_verified=saved.checksum_verified,workload_verified=true)
                break
            end
        end
        value === nothing && (value=_benchmark_performance(benchmark))
        value === nothing || (performance_checks=merge(performance_checks,(workload_verified=true,)))
        if performance_path !== nothing && value !== nothing
            temporary=tempname(directory)
            try
                JLD2.jldsave(temporary;schema_version=1,performance=value,session=performance_session,
                    calculations=(reference=value.reference.calculation,candidate=value.candidate.calculation))
                mv(temporary,performance_path)
                write(performance_path*".sha256",bytes2hex(open(sha256,performance_path))*"  performance.jld2\n")
            finally
                isfile(temporary) && rm(temporary)
            end
        end
        value
    end
    catch error
        error isa InterruptException && rethrow()
        performance_error=error
        @error "Performance measurement failed; continuing comparison report persistence" exception=(error,catch_backtrace())
        nothing
    end
    _performance_observation!(:reference;finished=true)
    _performance_observation!(:candidate;finished=true)
    measurements=(execution=(reference=(backend=string(_execution_label(benchmark.reference.formulation)),timing=reference_execution.timing,
            reused=reference_execution.reused,session=reference_execution.session,execution_wall_seconds=reference_execution.elapsed_seconds),
        candidate=(backend=string(_execution_label(benchmark.candidate.formulation)),timing=candidate_execution.timing,
            reused=candidate_execution.reused,session=candidate_execution.session,execution_wall_seconds=candidate_execution.elapsed_seconds)),
        performance,performance_checks...,session=performance_session)
    receiver === nothing || LineCableModels.report_progress(receiver,(stage=:reporting,))
    publication=report(definition,(
        reference=(result=reference,metadata=reference_metadata),
        candidate=(result=candidate,metadata=candidate_metadata),
        context=(id=benchmark.id,case_id=benchmark.case_id,collection=benchmark.collection),measurements))
    comparison=publication.published.comparisons
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
        receiver === nothing || LineCableModels.report_progress(receiver,(stage=:saving_report,))
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
    performance_error === nothing || throw(performance_error)
    timings=(scope = :execution_wall,
        execution = (
            reference = (seconds = reference_execution.elapsed_seconds,
                reused = reference_execution.reused || _result_reused(reference_execution.result),
                session=reference_execution.session,
                compute=reference_execution.timing),
            candidate = (seconds = candidate_execution.elapsed_seconds,
                reused = candidate_execution.reused || _result_reused(candidate_execution.result),
                session=candidate_execution.session,
                compute=candidate_execution.timing)),
        benchmark = performance)
    return (; mode, reference_result = reference_execution.result,
        candidate_result = candidate_execution.result,
        reference, candidate, comparison,
        performance, timings, metadata, artifact, report=publication)
end

"""
    benchmark_definition(model; id, source_file, reference, formulations, kwargs...)

Declare one case, one explicit reference and a scalar or Gridspace of candidate
formulations. The candidate follows the ordinary compute overload. ReportBuilder
owns quantities, bands and comparison controls.
"""
function benchmark_definition(model::LoadedCase; id::Symbol, source_file::AbstractString,
        reference, formulations, collection::Symbol=:manual,
        options::Union{NamedTuple,Grammar.ComputationOptions}=Grammar.ComputationOptions(), report=BenchmarkTableDefinition(), tolerances=(;))
    baseline=reference isa BenchmarkCalculation ? reference :
        BenchmarkCalculation(:reference,model.problem,reference)
    candidate=BenchmarkCalculation(:candidate,model.problem,formulations;options)
    return benchmark_definition(id,model.id,collection,source_file,model,baseline,candidate,
        report.settings,tolerances)
end
