"Return core results in Gridspace traversal order."
result(value::Union{LinearErrorResult, MonteCarloResult}) = value.values

Units.label(::Units.Quantity{:sample_count}) = "Count"
Units.symbol(::Units.Quantity{:sample_count}) = "n"
Units.label(::Units.Quantity{:probability}) = "Probability"
Units.symbol(::Units.Quantity{:probability}) = "p"
Units.label(::Units.Quantity{:cumulative_probability}) = "Cumulative probability"
Units.symbol(::Units.Quantity{:cumulative_probability}) = "F"
Units.label(::Units.Quantity{:probability_density}) = "Probability density"
Units.symbol(::Units.Quantity{:probability_density}) = "p"

const _DimensionlessStatisticalQuantity = Union{
    Units.Quantity{:sample_count},
    Units.Quantity{:probability},
    Units.Quantity{:cumulative_probability}
}
Units.native_unit(::_DimensionlessStatisticalQuantity) = Units.units(:base, :dimensionless)
Units.display_unit(::_DimensionlessStatisticalQuantity) = Units.units(:base, :dimensionless)

"Return the sample-summary product for each Monte Carlo point."
statistics(value::MonteCarloResult) = value.stats

"Return retained sample products, or `nothing` when retention was disabled."
samples(value::MonteCarloResult) = value.sample_values

"Return retained histogram products, or `nothing` when retention was disabled."
histograms(value::MonteCarloResult) = value.histogram_values

basis(value::MonteCarloResult) = basis(first(value.values))

"Return the uncertainty-bearing core results of a linear propagation."
uncertain_value(value::LinearErrorResult) = value.values

"""
$(TYPEDSIGNATURES)

Return the resolved root random seed of a Monte Carlo calculation.
"""
root_seed(value::MonteCarloResult) = value.root_seed

"""
$(TYPEDSIGNATURES)

Return the random seed used for the Gridspace point at `index`.
"""
point_seed(value::MonteCarloResult, index::Integer) = value.point_seeds[index]

"""
$(TYPEDSIGNATURES)

Return the number of trials evaluated for the Gridspace point at `index`.
"""
trial_count(value::MonteCarloResult, index::Integer) = value.trial_counts[index]

"""
$(TYPEDSIGNATURES)

Return the simultaneous empirical-CDF confidence of a Monte Carlo calculation.
"""
confidence(value::MonteCarloResult) = value.formulation.confidence

"""
$(TYPEDSIGNATURES)

Return the empirical-CDF tolerance used to size a Monte Carlo calculation.
"""
cdf_tolerance(value::MonteCarloResult) = value.formulation.cdf_tol

"""
$(TYPEDSIGNATURES)

Return the sampling distribution of a Monte Carlo calculation.
"""
sampling_distribution(value::MonteCarloResult) = value.formulation.distribution

const _MonteCarloProductSelector = Union{
    typeof(statistics),
    typeof(samples),
    typeof(histograms)
}

const _MonteCarloScientificSelector = Union{
    typeof(R),
    typeof(L),
    typeof(C),
    typeof(Engine.G)
}

const _StatisticSelector = Union{
    typeof(Statistics.mean),
    typeof(Statistics.std),
    typeof(Statistics.median),
    typeof(minimum),
    typeof(maximum)
}

function Units.quantity(
        ::_MonteCarloProductSelector,
        selector::_MonteCarloScientificSelector
)
    return Units.quantity(selector)
end

function _monte_carlo_product(value::MonteCarloResult, ::typeof(statistics), point::Integer)
    return value.stats[point]
end

function _monte_carlo_product(value::MonteCarloResult, ::typeof(samples), point::Integer)
    retained = value.sample_values
    retained === nothing && throw(ArgumentError("Monte Carlo samples were not retained"))
    return retained[point]
end

function _monte_carlo_product(value::MonteCarloResult, ::typeof(histograms), point::Integer)
    retained = value.histogram_values
    retained === nothing && throw(ArgumentError(
        "Monte Carlo histograms were not retained or derived",
    ))
    return retained[point]
end

function observe(
        value::MonteCarloResult,
        selector::_MonteCarloScientificSelector,
        point::Integer,
        indices...
)
    return observe(value.values[point], selector, indices...)
end

function observe(
        value::MonteCarloResult{<:Engine.LineParameters},
        ::typeof(frequencies),
        point::Integer,
        indices...
)
    return observe(value.values[point], frequencies, indices...)
end

function observe(
        value::MonteCarloResult,
        product::_MonteCarloProductSelector,
        selector::_MonteCarloScientificSelector,
        point::Integer,
        indices...
)
    stored = _monte_carlo_field(
        _monte_carlo_product(value, product, point),
        selector
    )
    return _product_value(stored, indices)
end

function observe(
        value::MonteCarloResult,
        ::typeof(statistics),
        selector::_MonteCarloScientificSelector,
        transform::_StatisticSelector,
        point::Integer,
        indices...
)
    stored = _monte_carlo_field(
        _monte_carlo_product(value, statistics, point),
        selector
    )
    return _statistic(transform, _product_value(stored, indices))
end

function _monte_carlo_observables(selectors::Tuple)
    product_selectors = (statistics, samples, histograms)
    products = Tuple((product, selector)
                     for product in product_selectors for selector in selectors)
    return (selectors..., products...)
end

function observables(
        ::Type{<:MonteCarloResult{T}}
) where {T <: DataModel.CableConstants}
    return _monte_carlo_observables((R, L, C))
end

function observables(
        ::Type{<:MonteCarloResult{T}}
) where {T <: Engine.LineParameters}
    return (frequencies, _monte_carlo_observables((R, L, C, Engine.G))...)
end

@inline _product_value(value, ::Tuple{}) = value
@inline _product_value(value, indices::Tuple) = getindex(value, indices...)

@inline _monte_carlo_field(product, ::typeof(R)) = product.R
@inline _monte_carlo_field(product, ::typeof(L)) = product.L
@inline _monte_carlo_field(product, ::typeof(C)) = product.C
@inline _monte_carlo_field(product, ::typeof(Engine.G)) = product.G

@inline _statistic(transform, value::AbstractArray) = map(transform, value)
@inline _statistic(transform, value) = transform(value)

function detach(summary::SampleSummary, factor)
    return SampleSummary(
        summary.mean * factor,
        summary.std * abs(factor),
        summary.min * factor,
        summary.q05 * factor,
        summary.median * factor,
        summary.q95 * factor,
        summary.max * factor,
        summary.n
    )
end

function detach(
        summaries::AbstractArray{<:SampleSummary},
        factor
)
    return map(summary -> detach(summary, factor), summaries)
end

function detach(histogram::HistogramDensity, factor)
    factor > zero(factor) || throw(ArgumentError("histogram conversion must be positive"))
    return HistogramDensity(
        histogram.edges .* factor,
        histogram.density ./ factor
    )
end
