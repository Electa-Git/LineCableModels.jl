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
        value::MonteCarloResult,
        product::_MonteCarloProductSelector,
        selector::_MonteCarloScientificSelector,
        point::Integer,
        indices...
)
    return observe(_monte_carlo_product(value, product, point), selector, indices...)
end

function observe(
        value::MonteCarloResult,
        ::typeof(statistics),
        selector::_MonteCarloScientificSelector,
        transform::_StatisticSelector,
        point::Integer,
        indices...
)
    return observe(
        _monte_carlo_product(value, statistics, point),
        selector,
        transform,
        indices...
    )
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
    return _monte_carlo_observables((R, L, C, Engine.G))
end

@inline _product_value(value, ::Tuple{}) = value
@inline _product_value(value, indices::Tuple) = getindex(value, indices...)

observe(product::RLCG, ::typeof(R), indices...) = _product_value(product.R, indices)
observe(product::RLCG, ::typeof(L), indices...) = _product_value(product.L, indices)
observe(product::RLCG, ::typeof(C), indices...) = _product_value(product.C, indices)
observe(product::RLCG, ::typeof(Engine.G), indices...) = _product_value(product.G, indices)
observe(product::RLC, ::typeof(R), indices...) = _product_value(product.R, indices)
observe(product::RLC, ::typeof(L), indices...) = _product_value(product.L, indices)
observe(product::RLC, ::typeof(C), indices...) = _product_value(product.C, indices)

@inline _statistic(transform, value::AbstractArray) = map(transform, value)
@inline _statistic(transform, value) = transform(value)

function observe(product::RLCG, ::typeof(R), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, R, indices...))
end
function observe(product::RLCG, ::typeof(L), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, L, indices...))
end
function observe(product::RLCG, ::typeof(C), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, C, indices...))
end
function observe(product::RLCG, ::typeof(Engine.G), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, Engine.G, indices...))
end
function observe(product::RLC, ::typeof(R), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, R, indices...))
end
function observe(product::RLC, ::typeof(L), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, L, indices...))
end
function observe(product::RLC, ::typeof(C), transform::_StatisticSelector, indices...)
    _statistic(transform, observe(product, C, indices...))
end

observables(::Type{<:RLC}) = (R, L, C)
observables(::Type{<:RLCG}) = (R, L, C, Engine.G)

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
