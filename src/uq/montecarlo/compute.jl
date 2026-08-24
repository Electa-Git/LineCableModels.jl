function _dkw_trials(observable_count::Integer, confidence::Real, cdf_tol::Real)
    observable_count > 0 || throw(ArgumentError("observable count must be positive"))
    alpha = 1 - confidence
    return ceil(Int, log(2 * observable_count / alpha) / (2 * cdf_tol^2))
end

_observable_count(value::DataModel.CableConstants) = length(observables(value))

function _observable_count(value::Engine.LineParameters)
    observed = observables(value)
    n = size(observed.series_impedance, 1)
    return 2 * n * (n + 1) * length(observed.frequency)
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

function _aggregate(draws::Vector{<:DataModel.CableConstants}, formulation::MonteCarlo)
    Rs = [value.R for value in draws]
    Ls = [value.L for value in draws]
    Cs = [value.C for value in draws]
    summaries = DataModel.CableConstants(SampleSummary(Rs), SampleSummary(Ls), SampleSummary(Cs))
    representation = DataModel.CableConstants(
        summaries.R.mean,
        summaries.L.mean,
        summaries.C.mean
    )
    retained = formulation.return_samples ? DataModel.CableConstants(Rs, Ls, Cs) : nothing
    hist = formulation.return_histograms ?
           DataModel.CableConstants(
        _histogram(Rs, formulation.bins),
        _histogram(Ls, formulation.bins),
        _histogram(Cs, formulation.bins)
    ) : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

function _rlcg_samples(draws::Vector{<:Engine.LineParameters})
    first_result = first(draws)
    first_observed = observables(first_result)
    dimensions = (size(first_observed.series_impedance)..., length(draws))
    Rs = Array{Float64}(undef, dimensions)
    Ls = similar(Rs)
    Cs = similar(Rs)
    Gs = similar(Rs)
    angular = reshape(2π .* first_observed.frequency, 1, 1, :)
    for (trial, value) in enumerate(draws)
        observed = observables(value)
        size(observed.series_impedance) == size(first_observed.series_impedance) ||
            throw(DimensionMismatch(
                "Monte Carlo realisations produced incompatible impedance dimensions",
            ))
        observed.frequency == first_observed.frequency || throw(DimensionMismatch(
            "Monte Carlo realisations produced incompatible frequency axes",
        ))
        @views Rs[:, :, :, trial] .= real.(observed.series_impedance.values)
        @views Ls[:, :, :, trial] .= imag.(observed.series_impedance.values) ./ angular
        @views Gs[:, :, :, trial] .= real.(observed.shunt_admittance.values)
        @views Cs[:, :, :, trial] .= imag.(observed.shunt_admittance.values) ./ angular
    end
    return RLCG(Rs, Ls, Cs, Gs)
end

function _map_observables(function_value, sample_values::Array{<:Real, 4})
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

function _aggregate(draws::Vector{<:Engine.LineParameters}, formulation::MonteCarlo)
    sample_values = _rlcg_samples(draws)
    summaries = RLCG((
        _map_observables(SampleSummary, values)
    for values in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...)
    means = RLCG((
        map(summary -> summary.mean, values)
    for values in (summaries.R, summaries.L, summaries.C, summaries.G)
    )...)
    first_result = first(draws)
    angular = reshape(2π .* first_result.f, 1, 1, :)
    representation = Engine.LineParameters(
        domain(first_result),
        complex.(means.R, means.L .* angular),
        complex.(means.G, means.C .* angular),
        first_result.f;
        basis = basis(first_result)
    )
    hist = formulation.return_histograms ?
           RLCG((
        _map_observables(values -> _histogram(values, formulation.bins), samples)
    for samples in (sample_values.R, sample_values.L, sample_values.C, sample_values.G)
    )...) : nothing
    retained = formulation.return_samples ? sample_values : nothing
    return (; representation, statistics = summaries, samples = retained, histograms = hist)
end

function _monte_carlo(point, formulation::MonteCarlo, options, seed)
    rng = Random.Xoshiro(seed)
    first_problem = ParametricBuilder.realize(rng, point, formulation.distribution)
    first_result = compute(first_problem, formulation.inner; options)
    ntrials = something(
        formulation.trials,
        _dkw_trials(_observable_count(first_result), formulation.confidence, formulation.cdf_tol)
    )
    draws = Vector{typeof(first_result)}(undef, ntrials)
    draws[1] = first_result
    for trial in 2:ntrials
        realization = ParametricBuilder.realize(rng, point, formulation.distribution)
        draws[trial] = compute(realization, formulation.inner; options)
    end
    return merge(_aggregate(draws, formulation), (; trials = ntrials, seed))
end

function compute(problem::ParametricProblem, formulation::MonteCarlo)
    values = nothing
    stats_values = nothing
    sample_values = nothing
    histogram_values = nothing
    seeds = UInt64[]
    trial_counts = Int[]
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
        trial_counts
    )
end
