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

function _monte_carlo(point, formulation::MonteCarlo, options, seed, details_owner)
    rng = Random.Xoshiro(seed)
    first_problem = ParametricBuilder.realize(rng, point, formulation.distribution)
    first_result = compute(first_problem, formulation.inner; options)
    ntrials = something(
        formulation.trials,
        _dkw_trials(_observable_count(first_result), formulation.confidence, formulation.cdf_tol)
    )
    retained = if formulation.options.retain_details
        first_record = computation_details(Val(details_owner), first_result)
        records = Vector{typeof(first_record)}(undef, ntrials)
        records[1] = first_record
        records
    else
        nothing
    end
    sample_values = _sample_storage(first_result, ntrials)
    sample_axis = _sample_axis(first_result)
    _record_sample!(sample_values, first_result, 1, sample_axis)
    for trial in 2:ntrials
        realization = ParametricBuilder.realize(rng, point, formulation.distribution)
        value = compute(realization, formulation.inner; options)
        typeof(value) === typeof(first_result) || throw(ArgumentError(
            "Monte Carlo realisations produced incompatible result types",
        ))
        _record_sample!(sample_values, value, trial, sample_axis)
        if retained !== nothing
            record = computation_details(Val(details_owner), value)
            typeof(record) === eltype(retained) || throw(ArgumentError(
                "Monte Carlo trials produced incompatible details record types",
            ))
            retained[trial] = record
        end
    end
    return merge(
        _aggregate(sample_values, first_result, formulation),
        (; trials = ntrials, seed, details = retained)
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

    return MonteCarloResult(
        formulation,
        values,
        stats_values,
        sample_values,
        histogram_values,
        root_seed,
        seeds,
        trial_counts,
        retained === nothing ? (;) : (trials = retained,)
    )
end
