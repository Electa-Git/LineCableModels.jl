function _dkw_trials(observable_count::Integer, confidence::Real, cdf_tol::Real)
    observable_count > 0 || throw(ArgumentError("observable count must be positive"))
    alpha = 1 - confidence
    return ceil(Int, log(2 * observable_count / alpha) / (2 * cdf_tol^2))
end

_observable_count(::DataModel.CableConstants) = 3

function _observable_count(value::Engine.LineParameters)
    n = size(observe(value, Engine.Z), 1)
    return 2 * n * (n + 1) * length(observe(value, Engine.frequencies))
end

function _sample_storage(first_result::DataModel.CableConstants, trials::Int)
    T = typeof(observe(first_result, R))
    return (
        R = Vector{T}(undef, trials),
        L = Vector{T}(undef, trials),
        C = Vector{T}(undef, trials)
    )
end

function _record_sample!(
        storage::NamedTuple{(:R, :L, :C)},
        value::DataModel.CableConstants,
        trial::Int,
        ::Nothing
)
    storage.R[trial] = observe(value, R)
    storage.L[trial] = observe(value, L)
    storage.C[trial] = observe(value, C)
    return storage
end

_sample_axis(::DataModel.CableConstants) = nothing
_sample_axis(value::Engine.LineParameters) = observe(value, Engine.frequencies)

function _aggregate(
        sample_values::NamedTuple{(:R, :L, :C)},
        ::DataModel.CableConstants,
        formulation::MonteCarlo
)
    Rs = sample_values.R
    Ls = sample_values.L
    Cs = sample_values.C
    summaries = (
        R = SampleSummary(Rs),
        L = SampleSummary(Ls),
        C = SampleSummary(Cs)
    )
    representation = DataModel.CableConstants(
        Statistics.mean(summaries.R),
        Statistics.mean(summaries.L),
        Statistics.mean(summaries.C)
    )
    retained = formulation.return_samples ? sample_values : nothing
    hist = formulation.return_histograms ?
           (
        R = HistogramDensity(Rs; bins = formulation.bins),
        L = HistogramDensity(Ls; bins = formulation.bins),
        C = HistogramDensity(Cs; bins = formulation.bins)
    ) : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

function _sample_storage(first_result::Engine.LineParameters, trials::Int)
    first_impedance = observe(first_result, Engine.Z)
    dimensions = (size(first_impedance)..., trials)
    Rs = Array{Float64}(undef, dimensions)
    Ls = similar(Rs)
    Cs = similar(Rs)
    Gs = similar(Rs)
    return (; R = Rs, L = Ls, C = Cs, G = Gs)
end

function _record_sample!(
        storage::NamedTuple{(:R, :L, :C, :G)},
        value::Engine.LineParameters,
        trial::Int,
        expected_frequencies
)
    impedance = observe(value, Engine.Z)
    size(impedance) == size(storage.R)[1:3] || throw(DimensionMismatch(
        "Monte Carlo realisations produced incompatible impedance dimensions",
    ))
    observe(value, Engine.frequencies) == expected_frequencies || throw(DimensionMismatch(
        "Monte Carlo realisations produced incompatible frequency axes",
    ))
    for index in CartesianIndices(size(impedance))
        i, j, k = index.I
        storage.R[i, j, k, trial] = observe(value, R, i, j, k)
        storage.L[i, j, k, trial] = observe(value, L, i, j, k)
        storage.G[i, j, k, trial] = observe(value, Engine.G, i, j, k)
        storage.C[i, j, k, trial] = observe(value, C, i, j, k)
    end
    return storage
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
    means = (
        R = Statistics.mean.(summaries.R),
        L = Statistics.mean.(summaries.L),
        C = Statistics.mean.(summaries.C),
        G = Statistics.mean.(summaries.G)
    )
    angular = reshape(2π .* observe(first_result, Engine.frequencies), 1, 1, :)
    representation = Engine.LineParameters(
        domain(first_result),
        complex.(means.R, means.L .* angular),
        complex.(means.G, means.C .* angular),
        observe(first_result, Engine.frequencies);
        basis = basis(first_result)
    )
    hist = formulation.return_histograms ?
           NamedTuple{(:R, :L, :C, :G)}(Tuple(
            _map_samples(
                values -> HistogramDensity(values; bins = formulation.bins),
                samples
            )
            for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
        )) : nothing
    retained = formulation.return_samples ? sample_values : nothing
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
        push!(summaries, (
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
    type_counts = [
        (type = error_type, count = by_type[error_type])
        for error_type in sort!(collect(keys(by_type)))
    ]
    stage_counts = [
        (stage, count = by_stage[stage])
        for stage in sort!(collect(keys(by_stage)))
    ]
    return (
        attempts,
        accepted,
        failed = length(failures),
        acceptance_rate = accepted / attempts,
        by_type = type_counts,
        by_stage = stage_counts
    )
end

function _retry_limit_error(failures, accepted::Int, attempts::Int, maximum::Int)
    summary = _failure_summary(failures, accepted, attempts)
    final_failure = last(failures)
    stack = final_failure.error.stack
    location = isempty(stack) ? "unknown location" :
               "$(first(stack).file):$(first(stack).line) in $(first(stack).function_name)"
    throw(ErrorException(
        "Monte Carlo retry limit of $maximum failures was exhausted after " *
        "$(summary.attempts) attempts ($(summary.accepted) accepted); final " *
        "$(final_failure.error.type) during $(final_failure.stage) at $location: " *
        final_failure.error.message,
    ))
end

function _monte_carlo(point, formulation::MonteCarlo, options, seed, details_owner)
    rng = Random.Xoshiro(seed)
    failures = NamedTuple[]
    attempts = 0
    accepted = 0
    ntrials = formulation.trials
    first_result = nothing
    sample_values = nothing
    sample_axis = nothing
    retained = nothing

    while ntrials === nothing || accepted < ntrials
        attempts += 1
        target_trial = accepted + 1
        sample = nothing
        stage = :sample
        value = nothing
        succeeded = false
        try
            sample = ParametricBuilder.realize_arguments(
                rng,
                point,
                formulation.distribution
            )
            stage = :build
            realization = ParametricBuilder.realize(point, sample)
            stage = :compute
            value = compute(realization, formulation.inner; options)
            succeeded = true
        catch exception
            backtrace = catch_backtrace()
            formulation.options.on_error === :retry && exception isa DomainError ||
                rethrow()
            push!(failures, _failure_record(
                attempts,
                target_trial,
                stage,
                sample,
                exception,
                backtrace
            ))
            length(failures) < formulation.options.max_failures ||
                _retry_limit_error(
                    failures,
                    accepted,
                    attempts,
                    formulation.options.max_failures
                )
        end
        succeeded || continue

        record = formulation.options.retain_details ?
                 computation_details(Val(details_owner), value) : nothing
        if accepted == 0
            first_result = value
            ntrials = something(
                ntrials,
                _dkw_trials(
                    _observable_count(first_result),
                    formulation.confidence,
                    formulation.cdf_tol
                )
            )
            sample_values = _sample_storage(first_result, ntrials)
            sample_axis = _sample_axis(first_result)
            if formulation.options.retain_details
                retained = Vector{typeof(record)}(undef, ntrials)
            end
        else
            typeof(value) === typeof(first_result) || throw(ArgumentError(
                "Monte Carlo realisations produced incompatible result types",
            ))
            if retained !== nothing
                typeof(record) === eltype(retained) || throw(ArgumentError(
                    "Monte Carlo trials produced incompatible details record types",
                ))
            end
        end

        accepted = target_trial
        _record_sample!(sample_values, value, accepted, sample_axis)
        retained === nothing || (retained[accepted] = record)
    end

    failure_summary = _failure_summary(failures, accepted, attempts)
    retained_details = retained === nothing ? nothing : (
        trials = retained,
        failures,
        failure_summary
    )
    return merge(
        _aggregate(sample_values, first_result, formulation),
        (; trials = ntrials, seed, details = retained_details)
    )
end

function compute(problem::ParametricProblem, formulation::MonteCarlo)
    point_count = length(problem.space)
    point_count > 0 || throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))
    root_seed = formulation.seed === nothing ? rand(Random.RandomDevice(), UInt64) :
                formulation.seed
    details_owner = formulation.options.retain_details ?
                    computation_owner(formulation.inner) : nothing
    point_source = ParametricBuilder.points(problem.space)
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
        details_owner
    )
    check_core_result(typeof(first_aggregate.representation))

    values = Vector{typeof(first_aggregate.representation)}(undef, point_count)
    stats_values = Vector{typeof(first_aggregate.statistics)}(undef, point_count)
    sample_values = formulation.return_samples ?
                    Vector{typeof(first_aggregate.samples)}(undef, point_count) : nothing
    histogram_values = formulation.return_histograms ?
                       Vector{typeof(first_aggregate.histograms)}(undef, point_count) :
                       nothing
    seeds = Vector{UInt64}(undef, point_count)
    trial_counts = Vector{Int}(undef, point_count)
    retained = formulation.options.retain_details ?
               Vector{typeof(first_aggregate.details)}(undef, point_count) : nothing

    values[1] = first_aggregate.representation
    stats_values[1] = first_aggregate.statistics
    sample_values === nothing || (sample_values[1] = first_aggregate.samples)
    histogram_values === nothing ||
        (histogram_values[1] = first_aggregate.histograms)
    seeds[1] = first_aggregate.seed
    trial_counts[1] = first_aggregate.trials
    retained === nothing || (retained[1] = first_aggregate.details)

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
            details_owner
        )
        typeof(aggregate.representation) === eltype(values) || throw(ArgumentError(
            "Monte Carlo points produced incompatible core result types",
        ))
        typeof(aggregate.statistics) === eltype(stats_values) || throw(ArgumentError(
            "Monte Carlo points produced incompatible statistics product types",
        ))
        values[index] = aggregate.representation
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
    end
    iterate(point_source, state) === nothing || throw(DimensionMismatch(
        "problem-space iteration exceeded its declared cardinality",
    ))

    retained_details = retained === nothing ? (;) : (
        trials = getproperty.(retained, :trials),
        failures = getproperty.(retained, :failures),
        failure_summary = getproperty.(retained, :failure_summary)
    )
    return MonteCarloResult(
        formulation,
        values,
        stats_values,
        sample_values,
        histogram_values,
        root_seed,
        seeds,
        trial_counts,
        retained_details
    )
end
