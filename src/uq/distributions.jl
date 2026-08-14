function HistogramPDF(
        edges::AbstractVector{TE},
        density::AbstractVector{TD}
) where {TE <: Real, TD <: Real}
    T = float(promote_type(TE, TD))
    return HistogramPDF{T}(Vector{T}(edges), Vector{T}(density))
end

function _auto_nbins(
        values::AbstractVector{<:Real};
        nbins_min::Int = 10,
        nbins_max::Int = 200
)
    isempty(values) && throw(ArgumentError("cannot bin an empty sample"))
    sorted_values = sort(float.(values))
    span = last(sorted_values) - first(sorted_values)
    if !(isfinite(span) && span > 0)
        return nbins_min
    end
    interquartile_range = StatsBase.quantile(sorted_values, 0.75) -
                          StatsBase.quantile(sorted_values, 0.25)
    if !(isfinite(interquartile_range) && interquartile_range > 0)
        return clamp(ceil(Int, sqrt(length(values))), nbins_min, nbins_max)
    end
    width = 2 * interquartile_range / cbrt(length(values))
    raw_count = span / width
    return isfinite(raw_count) ?
           clamp(ceil(Int, raw_count), nbins_min, nbins_max) : nbins_max
end

function _pdf_from_hist(
        values::AbstractVector{<:Real};
        nbins::Union{Int, Nothing} = nothing
)
    isempty(values) && throw(ArgumentError("cannot fit a histogram to an empty sample"))
    count = isnothing(nbins) ? _auto_nbins(values) : nbins
    count > 0 || throw(ArgumentError("nbins must be positive"))
    histogram = fit(Histogram, float.(values); nbins = count, closed = :left)
    edges = collect(histogram.edges[1])
    density = histogram.weights ./ (length(values) .* diff(edges))
    return HistogramPDF(edges, density)
end

function _binindex(distribution::HistogramPDF, value::Real)
    value < first(distribution.edges) && return 0
    value > last(distribution.edges) && return 0
    value == last(distribution.edges) && return length(distribution.density)
    return searchsortedlast(distribution.edges, value)
end

(distribution::HistogramPDF)(value::Real) = Distributions.pdf(distribution, value)
Distributions.minimum(distribution::HistogramPDF) = first(distribution.edges)
Distributions.maximum(distribution::HistogramPDF) = last(distribution.edges)
function Distributions.insupport(distribution::HistogramPDF, value::Real)
    minimum(distribution) <= value <= maximum(distribution)
end

function Distributions.pdf(distribution::HistogramPDF{T}, value::Real) where {T}
    index = _binindex(distribution, value)
    return index == 0 ? zero(T) : distribution.density[index]
end

function Distributions.logpdf(distribution::HistogramPDF, value::Real)
    probability = Distributions.pdf(distribution, value)
    return probability > 0 ? log(probability) : -Inf
end

function Distributions.cdf(distribution::HistogramPDF{T}, value::Real) where {T}
    value < minimum(distribution) && return zero(T)
    value >= maximum(distribution) && return one(T)
    index = _binindex(distribution, value)
    widths = diff(distribution.edges)
    prior = sum(
        (distribution.density[j] * widths[j] for j in 1:(index - 1));
        init = zero(T)
    )
    return prior + distribution.density[index] * (value - distribution.edges[index])
end

struct HistogramPDFSampler{T <: Real} <:
       Distributions.Sampleable{Distributions.Univariate, Distributions.Continuous}
    distribution::HistogramPDF{T}
    cumulative_probability::Vector{T}
end

function Distributions.sampler(distribution::HistogramPDF)
    probabilities = distribution.density .* diff(distribution.edges)
    cumulative = cumsum(probabilities)
    cumulative[end] = one(eltype(cumulative))
    return HistogramPDFSampler(distribution, cumulative)
end

function Distributions.quantile(sampler::HistogramPDFSampler{T}, probability::Real) where {T}
    probability <= 0 && return minimum(sampler.distribution)
    probability >= 1 && return maximum(sampler.distribution)
    index = searchsortedfirst(sampler.cumulative_probability, probability)
    prior = index == 1 ? zero(T) : sampler.cumulative_probability[index - 1]
    density = sampler.distribution.density[index]
    density == 0 && return sampler.distribution.edges[index]
    return sampler.distribution.edges[index] + (probability - prior) / density
end

function Distributions.quantile(distribution::HistogramPDF, probability::Real)
    0 <= probability <= 1 || throw(
        DomainError(probability, "probability must lie in [0, 1]"),
    )
    Distributions.quantile(Distributions.sampler(distribution), probability)
end

function Base.rand(rng::AbstractRNG, sampler::HistogramPDFSampler)
    probability = rand(rng)
    index = searchsortedfirst(sampler.cumulative_probability, probability)
    prior = index == 1 ? zero(probability) : sampler.cumulative_probability[index - 1]
    mass = sampler.cumulative_probability[index] - prior
    fraction = iszero(mass) ? zero(probability) : (probability - prior) / mass
    left = sampler.distribution.edges[index]
    right = sampler.distribution.edges[index + 1]
    return left + fraction * (right - left)
end

function Base.rand(rng::AbstractRNG, distribution::HistogramPDF)
    rand(rng, Distributions.sampler(distribution))
end

function _raw_moment(distribution::HistogramPDF{T}, order::Integer) where {T}
    order >= 0 || throw(ArgumentError("moment order must be nonnegative"))
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

function Distributions.moment(distribution::HistogramPDF, order::Integer)
    _raw_moment(distribution, order)
end
Distributions.mean(distribution::HistogramPDF) = _raw_moment(distribution, 1)
function Distributions.var(distribution::HistogramPDF)
    value = _raw_moment(distribution, 2) - Distributions.mean(distribution)^2
    return max(value, zero(value))
end
Distributions.std(distribution::HistogramPDF) = sqrt(Distributions.var(distribution))

function Distributions.mode(distribution::HistogramPDF)
    _, index = findmax(distribution.density)
    return (distribution.edges[index] + distribution.edges[index + 1]) / 2
end

function Distributions.modes(distribution::HistogramPDF)
    maximum_density = maximum(distribution.density)
    indices = findall(value -> value ≈ maximum_density, distribution.density)
    return [(distribution.edges[index] + distribution.edges[index + 1]) / 2
            for
            index in indices]
end

"""Return one retained scalar cable-constant trial."""
function trial(result::CableConstantsMC, index::Integer)
    resistance_values = samples(result, :R)
    index in eachindex(resistance_values) || throw(BoundsError(result, index))
    return CableConstants(
        resistance_values[index],
        samples(result, :L)[index],
        samples(result, :C)[index]
    )
end

"""Reconstruct one retained joint line-parameter trial."""
function trial(
        result::LineParametersMC{S, Samples, Distributions, Surrogate},
        index::Integer
) where {
        S,
        Samples,
        Distributions,
        T,
        U,
        D,
        Basis,
        Surrogate <: LineParameters{T, U, D, Basis}
}
    resistance_values = samples(result, :R)
    inductance_values = samples(result, :L)
    capacitance_values = samples(result, :C)
    conductance_values = samples(result, :G)
    sample_size = size(resistance_values)
    all(
        size(values) == sample_size
    for
    values in (inductance_values, capacitance_values, conductance_values)
    ) ||
        throw(DimensionMismatch("stored R, L, C, and G samples must have equal dimensions"))
    sample_size[1] == sample_size[2] || throw(
        DimensionMismatch("stored line-parameter matrices must be square"),
    )
    sample_size[3] == nfrequencies(result) || throw(
        DimensionMismatch("stored samples and frequencies must agree"),
    )
    index in axes(resistance_values, 4) || throw(BoundsError(result, index))

    matrix_count = sample_size[1]
    frequency_count = sample_size[3]
    scalar_type = eltype(resistance_values)
    impedance = Array{Complex{scalar_type}, 3}(
        undef,
        matrix_count,
        matrix_count,
        frequency_count
    )
    admittance = similar(impedance)
    frequency_values = frequencies(result)
    for k in 1:frequency_count, j in 1:matrix_count, i in 1:matrix_count
        omega = 2π * frequency_values[k]
        impedance[i, j, k] = resistance_values[i, j, k, index] +
                             im * omega * inductance_values[i, j, k, index]
        admittance[i, j, k] = conductance_values[i, j, k, index] +
                              im * omega * capacitance_values[i, j, k, index]
    end
    return LineParameters(
        D,
        SeriesImpedance{eltype(impedance), Basis}(impedance),
        ShuntAdmittance{eltype(admittance), Basis}(admittance),
        frequency_values
    )
end

function Base.rand(rng::AbstractRNG, result::Union{CableConstantsMC, LineParametersMC})
    has_samples(result) || throw(
        ArgumentError("joint sampling requires mc(...; return_samples=true)"),
    )
    return trial(result, rand(rng, 1:ntrials(result)))
end

function Base.rand(result::Union{CableConstantsMC, LineParametersMC})
    rand(Random.default_rng(), result)
end
