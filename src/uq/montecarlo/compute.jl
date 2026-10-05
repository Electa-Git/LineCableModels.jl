function _dkw_trials(observable_count::Integer, confidence::Real, cdf_tol::Real)
    observable_count > 0 || throw(ArgumentError("observable count must be positive"))
    alpha = 1 - confidence
    return ceil(Int, log(2 * observable_count / alpha) / (2 * cdf_tol^2))
end

_observable_count(value::Engine.CableConstants) = 4length(value)

function _observable_count(value::Engine.LineParameters)
    n = size(observe(value, Engine.Z), 1)
    return 2 * n * (n + 1) * length(observe(value, Engine.frequencies))
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with the `R`, `L`, `C` and `G` samples of `plan.trials` realizations,
one row per assembly of `input`, the first accepted result. The samples keep the scalar
type of that result.
"""
function initialize_buffers(::MonteCarlo, ::Type, input::Engine.CableConstants, plan,
        buffers)
    T = eltype(observe(input, R))
    return merge(buffers, (
        R = Matrix{T}(undef, length(input), plan.trials),
        L = Matrix{T}(undef, length(input), plan.trials),
        C = Matrix{T}(undef, length(input), plan.trials),
        G = Matrix{T}(undef, length(input), plan.trials)
    ))
end

function _record_sample!(
        buffers::NamedTuple{(:R, :L, :C, :G)},
        value::Engine.CableConstants,
        trial::Int,
        expected_axis
)
    value.cores == expected_axis.cores || throw(DimensionMismatch(
        "Monte Carlo cable-constant realizations produced incompatible assembly labels",
    ))
    value.frequency == expected_axis.frequency || throw(DimensionMismatch(
        "Monte Carlo cable-constant realizations produced incompatible frequencies",
    ))
    buffers.R[:, trial] .= observe(value, R)
    buffers.L[:, trial] .= observe(value, L)
    buffers.C[:, trial] .= observe(value, C)
    buffers.G[:, trial] .= observe(value, Engine.G)
    return buffers
end

function _sample_axis(value::Engine.CableConstants)
    (
        cores = copy(value.cores),
        frequency = value.frequency
    )
end
_sample_axis(value::Engine.LineParameters) = observe(value, Engine.frequencies)

function _cable_summaries(values::AbstractMatrix)
    return [SampleSummary(collect(@view values[assembly, :]))
            for assembly in axes(values, 1)]
end

function _sample_shunt_model(value::Union{Engine.LineParameters,Engine.CableConstants})
    model=get(details(value).data,:shunt_model,nothing)
    model===nothing && return ComputationDetails()
    return ComputationDetails(; shunt_model=(requested=model.requested,effective=model.effective,
        domains=[(design=d.design,terminals=collect(d.terminals),effective=d.effective)
            for d in model.domains]))
end

function _aggregate(
        sample_values::NamedTuple{(:R, :L, :C, :G)},
        first_result::Engine.CableConstants,
        formulation::MonteCarlo
)
    summaries = (
        R = _cable_summaries(sample_values.R),
        L = _cable_summaries(sample_values.L),
        C = _cable_summaries(sample_values.C),
        G = _cable_summaries(sample_values.G)
    )
    representation = materialize(first_result, summaries; details=_sample_shunt_model(first_result))
    retained = formulation.options.data.return_samples ? sample_values : nothing
    hist = formulation.options.data.return_histograms ?
           (
        R = [HistogramDensity(collect(@view sample_values.R[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.R, 1)],
        L = [HistogramDensity(collect(@view sample_values.L[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.L, 1)],
        C = [HistogramDensity(collect(@view sample_values.C[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.C, 1)],
        G = [HistogramDensity(collect(@view sample_values.G[assembly, :]);
                 bins = formulation.options.data.bins) for assembly in axes(sample_values.G, 1)]
    ) : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

"""
$(TYPEDSIGNATURES)

Extend `buffers` with the `R`, `L`, `C` and `G` samples of `plan.trials` realizations, in
the phase and frequency layout of `input`, the first accepted result, with scalar type `T`.
"""
function initialize_buffers(::MonteCarlo, ::Type{T}, input::Engine.LineParameters, plan,
        buffers) where {T}
    dimensions = (size(observe(input, Engine.Z))..., plan.trials)
    Rs = Array{T}(undef, dimensions)
    Ls = similar(Rs)
    Cs = similar(Rs)
    Gs = similar(Rs)
    return merge(buffers, (; R = Rs, L = Ls, C = Cs, G = Gs))
end

function _record_sample!(
        buffers::NamedTuple{(:R, :L, :C, :G)},
        value::Engine.LineParameters,
        trial::Int,
        expected_frequencies
)
    impedance = observe(value, Engine.Z)
    size(impedance) == size(buffers.R)[1:3] || throw(DimensionMismatch(
        "Monte Carlo realizations produced incompatible impedance dimensions",
    ))
    observe(value, Engine.frequencies) == expected_frequencies || throw(DimensionMismatch(
        "Monte Carlo realizations produced incompatible frequency axes",
    ))
    for index in CartesianIndices(size(impedance))
        i, j, k = index.I
        buffers.R[i, j, k, trial] = observe(value, R, i, j, k)
        buffers.L[i, j, k, trial] = observe(value, L, i, j, k)
        buffers.G[i, j, k, trial] = observe(value, Engine.G, i, j, k)
        buffers.C[i, j, k, trial] = observe(value, C, i, j, k)
    end
    return buffers
end

function _map_samples(function_value, sample_values::Array{<:Real, 4})
    output_size = size(sample_values)[1:3]
    indices = CartesianIndices(output_size)
    first_index = first(indices)
    first_value = function_value(collect(@view sample_values[first_index.I..., :]))
    output = Array{typeof(first_value)}(undef, output_size)
    output[first_index] = first_value
    for index in Iterators.drop(indices, 1)
        output[index] = function_value(collect(@view sample_values[index.I..., :]))
    end
    return output
end

function _aggregate(
        sample_values::NamedTuple{(:R, :L, :C, :G)},
        first_result::Engine.LineParameters,
        formulation::MonteCarlo
)
    summary_values = Tuple(
        _map_samples(SampleSummary, values)
    for values in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )
    summaries = NamedTuple{(:R, :L, :C, :G)}(summary_values)
    retained_details=_sample_shunt_model(first_result)
    source_details=details(first_result).data
    haskey(source_details,:coordinates) && (retained_details=ComputationDetails(
        merge(retained_details.data,(coordinates=source_details.coordinates,))))
    haskey(source_details,:comparison_unsupported) && (retained_details=ComputationDetails(
        merge(retained_details.data,(comparison_unsupported=source_details.comparison_unsupported,))))
    representation = materialize(first_result, summaries; details=retained_details)
    hist = formulation.options.data.return_histograms ?
           NamedTuple{(:R, :L, :C, :G)}(Tuple(
        _map_samples(
            values -> HistogramDensity(values; bins = formulation.options.data.bins),
            samples
        )
    for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )) : nothing
    retained = formulation.options.data.return_samples ? sample_values : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

const _MONTE_CARLO_STACKTRACE_LIMIT = 8

function _failure_stack(backtrace)
    summaries = @NamedTuple{
        function_name::String,
        file::String,
        line::Int
    }[]
    for frame in Iterators.take(
        Base.StackTraces.stacktrace(backtrace),
        _MONTE_CARLO_STACKTRACE_LIMIT
    )
        push!(summaries,
            (
                function_name = string(frame.func),
                file = string(frame.file),
                line = Int(frame.line)
            ))
    end
    return summaries
end

function _failure_record(
        attempt::Int,
        target_trial::Int,
        stage::Symbol,
        sample,
        exception,
        backtrace
)
    return (
        attempt,
        target_trial,
        stage,
        sample,
        error = (
            type = string(typeof(exception)),
            message = sprint(showerror, exception),
            stack = _failure_stack(backtrace)
        )
    )
end

function _failure_summary(failures, accepted::Int, attempts::Int)
    by_type = Dict{String, Int}()
    by_stage = Dict{Symbol, Int}()
    for failure in failures
        error_type = failure.error.type
        by_type[error_type] = get(by_type, error_type, 0) + 1
        by_stage[failure.stage] = get(by_stage, failure.stage, 0) + 1
    end
    type_counts = [(type = error_type, count = by_type[error_type])
                   for error_type in sort!(collect(keys(by_type)))]
    stage_counts = [(stage, count = by_stage[stage])
                    for stage in sort!(collect(keys(by_stage)))]
    return (
        attempts,
        accepted,
        failed = length(failures),
        acceptance_rate = accepted / attempts,
        by_type = type_counts,
        by_stage = stage_counts
    )
end

function _resampling_limit_error(failures, accepted::Int, attempts::Int, maximum::Int)
    summary = _failure_summary(failures, accepted, attempts)
    final_failure = last(failures)
    stack = final_failure.error.stack
    location = isempty(stack) ? "unknown location" :
               "$(first(stack).file):$(first(stack).line) in $(first(stack).function_name)"
    throw(ErrorException(
        "Monte Carlo resampling limit of $maximum failures was exhausted after " *
        "$(summary.attempts) attempts ($(summary.accepted) accepted); final " *
        "$(final_failure.error.type) during $(final_failure.stage) at $location: " *
        final_failure.error.message,
    ))
end

function realize(rng::Random.AbstractRNG, point::Gridpoint{Engine.LineParametersProblem}, distribution)
    return DataModel.realize_clearance(rng, point, distribution)
end

function _uncertain_arguments(point::Gridpoint)
    records=NamedTuple[]
    function collect!(value,path)
        if value isa UncertainValue
            push!(records,(argument_path=path,nominal=value.nominal,standard_deviation=value.sigma))
        elseif value isa Gridpoint
            for (index,argument) in enumerate(value.args)
                collect!(argument,(path...,index))
            end
        elseif value isa Union{Tuple,AbstractArray}
            for (index,argument) in enumerate(value)
                collect!(argument,(path...,index))
            end
        end
    end
    collect!(point,())
    return records
end

function _monte_carlo(point, formulation::MonteCarlo, options, seed, details_owner, child_options)
    clearance = point isa Gridpoint{Engine.LineParametersProblem} ?
                DataModel.collect_clearance_requirements(point) : nothing
    physical_inputs=clearance===nothing ? nothing : Engine.completed_inputs(clearance.declaration[])
    clearance===nothing || (clearance.declaration[]=nothing)
    try
        aggregate = _monte_carlo(point, formulation, options, seed, details_owner, child_options, clearance)
        if physical_inputs===nothing
            # Capture the nominal physical declaration once after successful
            # sampling. A Monte Carlo builder need not accept first-order numbers.
            declaration=realize(Random.Xoshiro(0),point,(_rng,mean,_sigma) -> mean)
            physical_inputs=merge(Engine.completed_inputs(declaration),(interpretation=:nominal_declaration,))
        end
        physical_inputs=merge(physical_inputs,(uncertain_arguments=_uncertain_arguments(point),))
        representation=Engine.retain_gridpoint(aggregate.representation,Commons.gridpoint_id();
            fields=merge(Engine.completed_formulation(formulation.inner),(inputs=physical_inputs,
                uncertainty=(estimator=:empirical,representation=:marginal_mean_std),
                uncertainty_descriptions=_uncertainty_descriptions(MonteCarlo))))
        return merge(aggregate,(;representation))
    finally
        DataModel.warn_clearance_summary(clearance)
    end
end

# Draw the arguments, build the problem inside the clearance scope and compute it, once.
# After a failure that the study resamples, add its record to `failures` and return
# `nothing`. `stage` and `drawn` keep the stage and the drawn arguments for that record.
function _attempt(rng, sampler, attempts::Int, accepted::Int)
    (; point, formulation, clearance, child_options, failures, stage, drawn) = sampler
    stage[] = :sample
    drawn[] = nothing
    try
        realization = DataModel.with_clearance(clearance) do
            arguments = realize_arguments(rng, point, formulation.options.data.distribution)
            drawn[] = arguments
            stage[] = :build
            realize(point, arguments)
        end
        stage[] = :compute
        return compute(realization, formulation.inner; options=child_options)
    catch exception
        backtrace = catch_backtrace()
        formulation.options.data.on_error === :resample && exception isa DomainError ||
            rethrow()
        push!(failures, _failure_record(
            attempts,
            accepted + 1,
            stage[],
            drawn[],
            exception,
            backtrace
        ))
        length(failures) < formulation.options.data.max_failures ||
            _resampling_limit_error(
                failures,
                accepted,
                attempts,
                formulation.options.data.max_failures
            )
        return nothing
    end
end

function _monte_carlo(point, formulation::MonteCarlo, options, seed, details_owner, child_options, clearance)
    rng = Random.Xoshiro(seed)
    failures = NamedTuple[]
    timing = get(options.data, :timing, false)
    timing isa Bool || throw(ArgumentError("timing must be Bool"))
    trial_timings = timing ? NamedTuple[] : nothing
    progress = point isa Gridpoint{<:Engine.LineParametersProblem} && verbosity(options, :progress) > 0
    started = progress ? time_ns() : UInt64(0)
    sampler = (; point, formulation, clearance, child_options, failures,
        stage = Ref(:sample), drawn = Ref{Any}(nothing))

    # Pilot: draw until the first accepted realization.
    attempts = 0
    pilot = nothing
    while pilot === nothing
        attempts += 1
        pilot = _attempt(rng, sampler, attempts, 0)
    end
    first_result = something(pilot)

    # The storage of every trial, allocated once from the pilot result.
    record = formulation.options.data.retain_details ?
             computation_details(details_owner, first_result) : nothing
    ntrials = something(
        formulation.options.data.trials,
        _dkw_trials(
            _observable_count(first_result),
            formulation.options.data.confidence,
            formulation.options.data.cdf_tol
        )
    )
    sample_values = initialize_buffers(formulation, Float64, first_result,
        (; trials = ntrials), (;))
    sample_axis = _sample_axis(first_result)
    retained = formulation.options.data.retain_details ?
               Vector{typeof(record)}(undef, ntrials) : nothing

    trials = _trials!(rng, first_result, sample_values, retained, record, sampler,
        details_owner, sample_axis, trial_timings, attempts, ntrials, progress, started)
    attempts = trials.attempts
    accepted = ntrials

    failure_summary = _failure_summary(failures, accepted, attempts)
    retained_details = retained === nothing ? nothing :
                       (
        trials = retained,
        failures,
        failure_summary,
        clearance = DataModel.clearance_summary(clearance)
    )
    aggregation_started = progress ? time_ns() : UInt64(0)
    aggregate = _aggregate(sample_values, first_result, formulation)
    if progress
        now = time_ns()
        if now - trials.last_log >= 5_000_000_000
            @info "Monte Carlo aggregation completed" _group=:progress accepted trials=ntrials attempts rejected=length(failures) aggregation_seconds=(now-aggregation_started)*1e-9 elapsed_seconds=(now-started)*1e-9
        end
    end
    return merge(aggregate,
        (; trials = ntrials, seed, details = retained_details, timing=trial_timings))
end

# Record the pilot result in `buffers`, then draw and record further trials until the study
# accepts `ntrials` of them. Return `attempts` and the time of the last progress message.
function _trials!(rng, first_result::R, buffers, retained, record, sampler, details_owner,
        sample_axis, trial_timings, attempts::Int, ntrials::Int, progress::Bool,
        started::UInt64) where {R}
    formulation = sampler.formulation
    failures = sampler.failures
    previous = started
    last_log = started
    average_seconds = 0.0
    accepted = 0
    value = first_result
    while true
        accepted += 1
        _record_sample!(buffers, value, accepted, sample_axis)
        retained === nothing || (retained[accepted] = record)
        trial_timings === nothing || push!(trial_timings, details(value).data.timing)
        if progress
            now = time_ns()
            interval = (now - previous) * 1e-9
            average_seconds = accepted == 1 ? interval : 0.2 * interval + 0.8 * average_seconds
            previous = now
            if now - last_log >= 5_000_000_000
                @info "Monte Carlo sampling progress" _group=:progress accepted trials=ntrials attempts rejected=length(failures) elapsed_seconds=(now-started)*1e-9 eta_hours=(ntrials-accepted)*average_seconds/3600
                last_log = now
            end
        end
        accepted < ntrials || break
        result = nothing
        while result === nothing
            attempts += 1
            result = _attempt(rng, sampler, attempts, accepted)
        end
        record = formulation.options.data.retain_details ?
                 computation_details(details_owner, result) : nothing
        _sample_shunt_model(result) == _sample_shunt_model(first_result) || throw(ArgumentError(
            "Monte Carlo realizations changed shunt-model coverage; select a fixed geometry model for this study"))
        result isa R || throw(ArgumentError(
            "Monte Carlo realizations produced incompatible result types",
        ))
        if retained !== nothing
            typeof(record) === eltype(retained) || throw(ArgumentError(
                "Monte Carlo trials produced incompatible details record types",
            ))
        end
        value = result
    end
    return (; attempts, last_log)
end

function compute(problem::ParametricProblem, formulation::MonteCarlo)
    levels = verbosity(get(problem.options.data, :verbosity, (default=0,)))
    progress = problem.space isa LineCableModels.Gridspace{<:Engine.LineParametersProblem} && get(levels, :progress, levels.default) > 0
    child_options = progress ? ComputationOptions(merge(problem.options.data,
        (verbosity=merge(problem.options.data.verbosity, (progress=0,)),))) : problem.options
    logger = haskey(problem.options.data, :verbosity) ? VerbosityLogger(Logging.current_logger(), levels) : Logging.current_logger()
    return Logging.with_logger(logger) do
    started = progress ? time_ns() : UInt64(0)
    previous = started
    last_log = started
    average_seconds = 0.0
    source_id=Commons.gridpoint_id().source_id
    Base.get_extension(LineCableModels, :LineCableModelsMeasurementsExt) === nothing &&
        throw(ArgumentError("MonteCarlo requires the Measurements extension to construct its result; " *
            "load it with `using Measurements` before compute"))
    point_count = length(problem.space)
    point_count > 0 || throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))
    progress && @info "Monte Carlo computation started" _group=:progress populations=point_count trials=formulation.options.data.trials
    root_seed = formulation.options.data.seed === nothing ? rand(Random.RandomDevice(), UInt64) :
                formulation.options.data.seed
    details_owner = formulation.options.data.retain_details ?
                    typeof(formulation.inner) : nothing
    point_source = points(problem.space)
    first_item = iterate(point_source)
    first_item === nothing && throw(DimensionMismatch(
        "problem-space iteration ended before its declared cardinality",
    ))
    first_point, state = first_item
    first_seed = root_seed
    first_aggregate = _monte_carlo(
        first_point,
        formulation,
        problem.options,
        first_seed,
        details_owner,
        child_options
    )
    validate(typeof(first_aggregate.representation), MonteCarloResult)

    values = Vector{typeof(first_aggregate.representation)}(undef, point_count)
    stats_values = Vector{typeof(first_aggregate.statistics)}(undef, point_count)
    sample_values = formulation.options.data.return_samples ?
                    Vector{typeof(first_aggregate.samples)}(undef, point_count) : nothing
    histogram_values = formulation.options.data.return_histograms ?
                       Vector{typeof(first_aggregate.histograms)}(undef, point_count) :
                       nothing
    seeds = Vector{UInt64}(undef, point_count)
    trial_counts = Vector{Int}(undef, point_count)
    retained = formulation.options.data.retain_details ?
               Vector{typeof(first_aggregate.details)}(undef, point_count) : nothing

    values[1] = Engine.retain_gridpoint(first_aggregate.representation,Commons.gridpoint_id(;source_id))
    stats_values[1] = first_aggregate.statistics
    sample_values === nothing || (sample_values[1] = first_aggregate.samples)
    histogram_values === nothing ||
        (histogram_values[1] = first_aggregate.histograms)
    seeds[1] = first_aggregate.seed
    trial_counts[1] = first_aggregate.trials
    retained === nothing || (retained[1] = first_aggregate.details)
    timing = get(problem.options.data, :timing, false)
    population_timings = timing ? Vector{Vector{NamedTuple}}(undef, point_count) : nothing
    timing && (population_timings[1]=first_aggregate.timing)
    if progress
        now = time_ns()
        average_seconds = (now - previous) * 1e-9
        previous = now
        if now - last_log >= 5_000_000_000
            @info "Monte Carlo population progress" _group=:progress completed=1 populations=point_count elapsed_seconds=(now-started)*1e-9 eta_hours=(point_count-1)*average_seconds/3600
            last_log = now
        end
    end

    for index in 2:point_count
        item = iterate(point_source, state)
        item === nothing && throw(DimensionMismatch(
            "problem-space iteration ended before its declared cardinality",
        ))
        point, state = item
        point_seed = root_seed ⊻ (UInt64(index - 1) * 0x9e3779b97f4a7c15)
        aggregate = _monte_carlo(
            point,
            formulation,
            problem.options,
            point_seed,
            details_owner,
            child_options
        )
        typeof(aggregate.representation) === eltype(values) || throw(ArgumentError(
            "Monte Carlo points produced incompatible core result types",
        ))
        typeof(aggregate.statistics) === eltype(stats_values) || throw(ArgumentError(
            "Monte Carlo points produced incompatible statistics product types",
        ))
        values[index] = Engine.retain_gridpoint(aggregate.representation,
            Commons.gridpoint_id(;source_id,problem_index=index))
        stats_values[index] = aggregate.statistics
        if sample_values !== nothing
            typeof(aggregate.samples) === eltype(sample_values) || throw(ArgumentError(
                "Monte Carlo points produced incompatible sample product types",
            ))
            sample_values[index] = aggregate.samples
        end
        if histogram_values !== nothing
            typeof(aggregate.histograms) === eltype(histogram_values) ||
                throw(ArgumentError(
                    "Monte Carlo points produced incompatible histogram product types",
                ))
            histogram_values[index] = aggregate.histograms
        end
        seeds[index] = aggregate.seed
        trial_counts[index] = aggregate.trials
        if retained !== nothing
            typeof(aggregate.details) === eltype(retained) || throw(ArgumentError(
                "Monte Carlo points produced incompatible details product types",
            ))
            retained[index] = aggregate.details
        end
        timing && (population_timings[index]=aggregate.timing)
        if progress
            now = time_ns()
            interval = (now - previous) * 1e-9
            average_seconds = 0.2 * interval + 0.8 * average_seconds
            previous = now
            if now - last_log >= 5_000_000_000
                @info "Monte Carlo population progress" _group=:progress completed=index populations=point_count elapsed_seconds=(now-started)*1e-9 eta_hours=(point_count-index)*average_seconds/3600
                last_log = now
            end
        end
    end
    iterate(point_source, state) === nothing || throw(DimensionMismatch(
        "problem-space iteration exceeded its declared cardinality",
    ))

    retained_details = retained === nothing ? (;) :
                       (
        trials = getproperty.(retained, :trials),
        failures = getproperty.(retained, :failures),
        failure_summary = getproperty.(retained, :failure_summary),
        clearance = getproperty.(retained, :clearance)
    )
    timing && (retained_details=merge(retained_details, (timing=population_timings,)))
    result = MonteCarloResult(
        formulation,
        values,
        stats_values,
        sample_values,
        histogram_values,
        root_seed,
        seeds,
        trial_counts,
        ComputationDetails(retained_details)
    )
    progress && @info "Monte Carlo computation completed successfully" _group=:progress completed=point_count populations=point_count elapsed_seconds=(time_ns()-started)*1e-9
    result
    end
end
