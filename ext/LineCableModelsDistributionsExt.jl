module LineCableModelsDistributionsExt

using Distributions
using Random
using Statistics

import LineCableModels
const PB = LineCableModels.ParametricBuilder
const Computation = LineCableModels.Computation

function PB._sample_uncertainty(
    rng::Random.AbstractRNG,
    value::PB.UncertainValue{<:Real},
    distribution::Distributions.UnivariateDistribution,
)
    standardized = rand(rng, distribution)
    distribution_mean = Statistics.mean(distribution)
    distribution_std = Statistics.std(distribution)
    isfinite(distribution_mean) && isfinite(distribution_std) &&
        distribution_std > zero(distribution_std) || throw(ArgumentError(
        "Monte Carlo distributions must have finite mean and positive finite standard deviation",
    ))
    isfinite(standardized) || throw(ArgumentError(
        "Monte Carlo distribution produced a non-finite realization",
    ))
    return value.nominal + value.sigma *
        (standardized - distribution_mean) / distribution_std
end

function Distributions.pdf(distribution::Computation.HistogramPDF, value::Real)
    value < first(distribution.edges) && return 0.0
    value > last(distribution.edges) && return 0.0
    index = value == last(distribution.edges) ? length(distribution.density) :
        searchsortedlast(distribution.edges, value)
    return distribution.density[index]
end

(distribution::Computation.HistogramPDF)(value::Real) =
    Distributions.pdf(distribution, value)
Distributions.minimum(distribution::Computation.HistogramPDF) =
    first(distribution.edges)
Distributions.maximum(distribution::Computation.HistogramPDF) =
    last(distribution.edges)
Distributions.insupport(
    distribution::Computation.HistogramPDF,
    value::Real,
) = minimum(distribution) <= value <= maximum(distribution)

function Distributions.logpdf(
    distribution::Computation.HistogramPDF,
    value::Real,
)
    probability = Distributions.pdf(distribution, value)
    return probability > 0 ? log(probability) : -Inf
end

function Distributions.cdf(distribution::Computation.HistogramPDF, value::Real)
    value <= first(distribution.edges) && return 0.0
    value >= last(distribution.edges) && return 1.0
    index = searchsortedlast(distribution.edges, value)
    widths = diff(distribution.edges)
    prior = index == 1 ? 0.0 :
        sum(distribution.density[1:(index - 1)] .* widths[1:(index - 1)])
    return prior + distribution.density[index] *
        (value - distribution.edges[index])
end

struct HistogramPDFSampler{H,T} <:
       Distributions.Sampleable{Distributions.Univariate,Distributions.Continuous}
    distribution::H
    cumulative_probability::Vector{T}
end

function Distributions.sampler(distribution::Computation.HistogramPDF)
    probabilities = distribution.density .* diff(distribution.edges)
    cumulative = cumsum(probabilities)
    cumulative[end] = one(eltype(cumulative))
    return HistogramPDFSampler(distribution, cumulative)
end

function Distributions.quantile(sampler::HistogramPDFSampler, probability::Real)
    probability <= 0 && return minimum(sampler.distribution)
    probability >= 1 && return maximum(sampler.distribution)
    index = searchsortedfirst(sampler.cumulative_probability, probability)
    prior = index == 1 ? zero(eltype(sampler.cumulative_probability)) :
        sampler.cumulative_probability[index - 1]
    density = sampler.distribution.density[index]
    iszero(density) && return sampler.distribution.edges[index]
    return sampler.distribution.edges[index] + (probability - prior) / density
end

function Distributions.quantile(
    distribution::Computation.HistogramPDF,
    probability::Real,
)
    0 <= probability <= 1 || throw(DomainError(
        probability,
        "probability must lie in [0, 1]",
    ))
    return Distributions.quantile(Distributions.sampler(distribution), probability)
end

function Base.rand(rng::Random.AbstractRNG, sampler::HistogramPDFSampler)
    probability = rand(rng)
    index = searchsortedfirst(sampler.cumulative_probability, probability)
    prior = index == 1 ? zero(probability) :
        sampler.cumulative_probability[index - 1]
    mass = sampler.cumulative_probability[index] - prior
    fraction = iszero(mass) ? zero(probability) : (probability - prior) / mass
    left = sampler.distribution.edges[index]
    right = sampler.distribution.edges[index + 1]
    return left + fraction * (right - left)
end

Base.rand(rng::Random.AbstractRNG, distribution::Computation.HistogramPDF) =
    rand(rng, Distributions.sampler(distribution))

function _raw_moment(distribution::Computation.HistogramPDF, order::Integer)
    order >= 0 || throw(ArgumentError("moment order must be nonnegative"))
    T = eltype(distribution.density)
    total = zero(T)
    exponent = order + 1
    for index in eachindex(distribution.density)
        left = distribution.edges[index]
        right = distribution.edges[index + 1]
        total += distribution.density[index] *
            (right^exponent - left^exponent) / exponent
    end
    return total
end

Distributions.moment(distribution::Computation.HistogramPDF, order::Integer) =
    _raw_moment(distribution, order)
Statistics.mean(distribution::Computation.HistogramPDF) =
    _raw_moment(distribution, 1)
function Statistics.var(distribution::Computation.HistogramPDF)
    variance = _raw_moment(distribution, 2) - Statistics.mean(distribution)^2
    return max(variance, zero(variance))
end
Statistics.std(distribution::Computation.HistogramPDF) =
    sqrt(Statistics.var(distribution))

function Distributions.mode(distribution::Computation.HistogramPDF)
    _, index = findmax(distribution.density)
    return (distribution.edges[index] + distribution.edges[index + 1]) / 2
end

function Distributions.modes(distribution::Computation.HistogramPDF)
    maximum_density = maximum(distribution.density)
    indices = findall(==(maximum_density), distribution.density)
    return [
        (distribution.edges[index] + distribution.edges[index + 1]) / 2
        for index in indices
    ]
end

end
