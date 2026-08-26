"""
    LineCableModelsDistributionsExt

Sample `UncertainValue` with Distributions.jl and expose `HistogramDensity`
through the Distributions API.
"""
module LineCableModelsDistributionsExt

using Distributions
using Random
using Statistics

import LineCableModels
const PB = LineCableModels.ParametricBuilder
const UQ = LineCableModels.UQ
import LineCableModels.ParametricBuilder: sample_uncertainty

function sample_uncertainty(
        rng::Random.AbstractRNG,
        value::PB.UncertainValue{<:Real},
        distribution::Distributions.UnivariateDistribution
)
    standardized = rand(rng, distribution)
    distribution_mean = Statistics.mean(distribution)
    distribution_std = Statistics.std(distribution)
    isfinite(distribution_mean) && isfinite(distribution_std) &&
    distribution_std > zero(distribution_std) || throw(ArgumentError(
        "Monte Carlo distributions must have finite mean and positive finite standard deviation",
    ))
    isfinite(standardized) || throw(ArgumentError(
        "Monte Carlo distribution produced a non-finite realisation",
    ))
    return value.nominal +
           value.sigma *
           (standardized - distribution_mean) / distribution_std
end

function Distributions.pdf(distribution::UQ.HistogramDensity, value::Real)
    value < first(distribution.edges) && return 0.0
    value > last(distribution.edges) && return 0.0
    index = value == last(distribution.edges) ? length(distribution.density) :
            searchsortedlast(distribution.edges, value)
    return distribution.density[index]
end

function (distribution::UQ.HistogramDensity)(value::Real)
    Distributions.pdf(distribution, value)
end
function Distributions.minimum(distribution::UQ.HistogramDensity)
    first(distribution.edges)
end
Distributions.maximum(distribution::UQ.HistogramDensity) = last(distribution.edges)
function Distributions.insupport(
        distribution::UQ.HistogramDensity,
        value::Real
)
    minimum(distribution) <= value <= maximum(distribution)
end

function Distributions.logpdf(
        distribution::UQ.HistogramDensity,
        value::Real
)
    probability = Distributions.pdf(distribution, value)
    return probability > 0 ? log(probability) : -Inf
end

function Distributions.cdf(distribution::UQ.HistogramDensity, value::Real)
    return UQ.cumulative_probability(distribution, value)
end

struct HistogramDensitySampler{H, T} <:
       Distributions.Sampleable{Distributions.Univariate, Distributions.Continuous}
    distribution::H
    cumulative_probability::Vector{T}
end

function Distributions.sampler(distribution::UQ.HistogramDensity)
    probabilities = distribution.density .* diff(distribution.edges)
    cumulative = cumsum(probabilities)
    cumulative[end] = one(eltype(cumulative))
    return HistogramDensitySampler(distribution, cumulative)
end

function Distributions.quantile(sampler::HistogramDensitySampler, probability::Real)
    probability <= 0 && return minimum(sampler.distribution)
    probability >= 1 && return maximum(sampler.distribution)
    index = searchsortedfirst(sampler.cumulative_probability, probability)
    prior = index == 1 ? zero(eltype(sampler.cumulative_probability)) :
            sampler.cumulative_probability[index - 1]
    density = sampler.distribution.density[index]
    iszero(density) && return sampler.distribution.edges[index]
    return sampler.distribution.edges[index] + (probability - prior) / density
end

function Base.rand(rng::Random.AbstractRNG, sampler::HistogramDensitySampler)
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

function Base.rand(rng::Random.AbstractRNG, distribution::UQ.HistogramDensity)
    rand(rng, Distributions.sampler(distribution))
end

function _raw_moment(distribution::UQ.HistogramDensity, order::Integer)
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

function Distributions.moment(distribution::UQ.HistogramDensity, order::Integer)
    _raw_moment(distribution, order)
end
Statistics.mean(distribution::UQ.HistogramDensity) = _raw_moment(distribution, 1)
function Statistics.var(distribution::UQ.HistogramDensity)
    variance = _raw_moment(distribution, 2) - Statistics.mean(distribution)^2
    return max(variance, zero(variance))
end
function Statistics.std(distribution::UQ.HistogramDensity)
    sqrt(Statistics.var(distribution))
end

function Distributions.mode(distribution::UQ.HistogramDensity)
    _, index = findmax(distribution.density)
    return (distribution.edges[index] + distribution.edges[index + 1]) / 2
end

function Distributions.modes(distribution::UQ.HistogramDensity)
    maximum_density = maximum(distribution.density)
    indices = findall(==(maximum_density), distribution.density)
    return [(distribution.edges[index] + distribution.edges[index + 1]) / 2
            for index in indices]
end

end
