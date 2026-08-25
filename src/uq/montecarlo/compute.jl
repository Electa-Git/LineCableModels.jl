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

function _histogram(values::AbstractVector{<:Real}, bins::Union{Nothing, Int})
    lo, hi = extrema(values)
    if lo == hi
        width = max(abs(float(lo)) * sqrt(eps(Float64)), sqrt(eps(Float64)))
        edges = [float(lo) - width, float(hi) + width]
    else
        count = something(bins, max(1, ceil(Int, sqrt(length(values)))))
        edges = collect(range(float(lo), float(hi); length = count + 1))
    end
    counts = zeros(Float64, length(edges) - 1)
    for value in values
        index = value == last(edges) ? length(counts) :
                clamp(searchsortedlast(edges, value), 1, length(counts))
        counts[index] += 1
    end
    return HistogramDensity(edges, counts ./ (length(values) .* diff(edges)))
end

function _sample_storage(first_result::DataModel.CableConstants, trials::Int)
    T = typeof(observe(first_result, R))
    return DataModel.CableConstants(
        Vector{T}(undef, trials),
        Vector{T}(undef, trials),
        Vector{T}(undef, trials)
    )
end

function _record_sample!(
        storage::DataModel.CableConstants,
        value::DataModel.CableConstants,
        trial::Int,
        ::Nothing
)
    observe(storage, R)[trial] = observe(value, R)
    observe(storage, L)[trial] = observe(value, L)
    observe(storage, C)[trial] = observe(value, C)
    return storage
end

_sample_axis(::DataModel.CableConstants) = nothing
_sample_axis(value::Engine.LineParameters) = observe(value, Engine.frequencies)

function _aggregate(
        sample_values::DataModel.CableConstants,
        ::DataModel.CableConstants,
        formulation::MonteCarlo
)
    Rs = observe(sample_values, R)
    Ls = observe(sample_values, L)
    Cs = observe(sample_values, C)
    summaries = DataModel.CableConstants(SampleSummary(Rs), SampleSummary(Ls), SampleSummary(Cs))
    representation = DataModel.CableConstants(
        observe(summaries, R, Statistics.mean),
        observe(summaries, L, Statistics.mean),
        observe(summaries, C, Statistics.mean)
    )
    retained = formulation.return_samples ? sample_values : nothing
    hist = formulation.return_histograms ?
           DataModel.CableConstants(
        _histogram(Rs, formulation.bins),
        _histogram(Ls, formulation.bins),
        _histogram(Cs, formulation.bins)
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
    return RLCG(Rs, Ls, Cs, Gs; basis = basis(first_result))
end

function _record_sample!(
        storage::RLCG,
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
    basis(value) === basis(storage) || throw(ArgumentError(
        "Monte Carlo realisations produced incompatible result bases",
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
        sample_values::RLCG,
        first_result::Engine.LineParameters,
        formulation::MonteCarlo
)
    summary_values = (
        _map_samples(SampleSummary, values)
    for values in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )
    summaries = RLCG(summary_values...; basis = basis(sample_values))
    means = RLCG(
        observe(summaries, R, Statistics.mean),
        observe(summaries, L, Statistics.mean),
        observe(summaries, C, Statistics.mean),
        observe(summaries, Engine.G, Statistics.mean);
        basis = basis(summaries)
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
           RLCG((
        _map_samples(values -> _histogram(values, formulation.bins), samples)
    for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...; basis = basis(sample_values)) : nothing
    retained = formulation.return_samples ? sample_values : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

function _monte_carlo(point, formulation::MonteCarlo, options, seed)
    rng = Random.Xoshiro(seed)
    first_problem = ParametricBuilder.realize(rng, point, formulation.distribution)
    first_result = compute(first_problem, formulation.inner; options)
    retained = if formulation.options.retain_details
        [computation_details(
            Val(ParametricBuilder._computation_owner(formulation.inner)),
            first_result
        )]
    else
        nothing
    end
    ntrials = something(
        formulation.trials,
        _dkw_trials(_observable_count(first_result), formulation.confidence, formulation.cdf_tol)
    )
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
        formulation.options.retain_details && push!(
            retained,
            computation_details(
                Val(ParametricBuilder._computation_owner(formulation.inner)),
                value
            )
        )
    end
    return merge(
        _aggregate(sample_values, first_result, formulation),
        (; trials = ntrials, seed, details = retained)
    )
end

function compute(problem::ParametricProblem, formulation::MonteCarlo)
    values = nothing
    stats_values = nothing
    sample_values = nothing
    histogram_values = nothing
    seeds = UInt64[]
    trial_counts = Int[]
    retained = nothing
    root_seed = formulation.seed === nothing ? rand(Random.RandomDevice(), UInt64) :
                formulation.seed
    for (index, point) in enumerate(ParametricBuilder.points(problem.space))
        point_seed = root_seed ⊻ (UInt64(index - 1) * 0x9e3779b97f4a7c15)
        aggregate = _monte_carlo(
            point,
            formulation,
            problem.options,
            point_seed
        )
        values = ParametricBuilder._append_result(values, aggregate.representation)
        stats_values = ParametricBuilder._append_result(stats_values, aggregate.statistics)
        formulation.return_samples &&
            (sample_values = ParametricBuilder._append_result(
                sample_values,
                aggregate.samples
            ))
        formulation.return_histograms &&
            (histogram_values = ParametricBuilder._append_result(
                histogram_values,
                aggregate.histograms
            ))
        push!(seeds, aggregate.seed)
        push!(trial_counts, aggregate.trials)
        formulation.options.retain_details &&
            (retained = ParametricBuilder._append_result(retained, aggregate.details))
    end
    typed_values = values === nothing ? AbstractProblemResult[] : values
    typed_stats = stats_values === nothing ? Any[] : stats_values
    return MonteCarloResult(
        formulation,
        typed_values,
        typed_stats,
        formulation.return_samples ? something(sample_values, Any[]) : nothing,
        formulation.return_histograms ? something(histogram_values, Any[]) : nothing,
        root_seed,
        seeds,
        trial_counts,
        formulation.options.retain_details ?
        (trials = something(retained, Vector{NamedTuple}[]),) : (;)
    )
end
