"Return primitive results in Gridspace traversal order."
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

"Return the uncertainty-bearing primitive results of a linear propagation."
uncertain_value(value::LinearErrorResult) = value.values

@inline _product_value(value, ::Tuple{}) = value
@inline _product_value(value, indices::Tuple) = getindex(value, indices...)

observe(product::RLCG, ::typeof(R), indices...) = _product_value(product.R, indices)
observe(product::RLCG, ::typeof(L), indices...) = _product_value(product.L, indices)
observe(product::RLCG, ::typeof(C), indices...) = _product_value(product.C, indices)
observe(product::RLCG, ::typeof(Engine.G), indices...) = _product_value(product.G, indices)

@inline _statistic(transform, value::AbstractArray) = map(transform, value)
@inline _statistic(transform, value) = transform(value)

const _StatisticSelector = Union{
    typeof(Statistics.mean),
    typeof(Statistics.std),
    typeof(Statistics.median),
    typeof(minimum),
    typeof(maximum)
}

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

const _CableUQValue = Union{SampleSummary, HistogramDensity, AbstractArray}

function observe(
        product::DataModel.CableConstants{<:_CableUQValue},
        ::typeof(R),
        indices...
)
    _product_value(product.R, indices)
end
function observe(
        product::DataModel.CableConstants{<:_CableUQValue},
        ::typeof(L),
        indices...
)
    _product_value(product.L, indices)
end
function observe(
        product::DataModel.CableConstants{<:_CableUQValue},
        ::typeof(C),
        indices...
)
    _product_value(product.C, indices)
end

function observe(
        product::DataModel.CableConstants{<:_CableUQValue},
        ::typeof(R),
        transform::_StatisticSelector,
        indices...
)
    _statistic(transform, observe(product, R, indices...))
end
function observe(
        product::DataModel.CableConstants{<:_CableUQValue},
        ::typeof(L),
        transform::_StatisticSelector,
        indices...
)
    _statistic(transform, observe(product, L, indices...))
end
function observe(
        product::DataModel.CableConstants{<:_CableUQValue},
        ::typeof(C),
        transform::_StatisticSelector,
        indices...
)
    _statistic(transform, observe(product, C, indices...))
end

observables(::Type{<:RLCG}) = (R, L, C, Engine.G)

function Grammar._detach_and_scale(summary::SampleSummary, factor)
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

function Grammar._detach_and_scale(
        summaries::AbstractArray{<:SampleSummary},
        factor
)
    return map(summary -> Grammar._detach_and_scale(summary, factor), summaries)
end

function Grammar._detach_and_scale(histogram::HistogramDensity, factor)
    factor > zero(factor) || throw(ArgumentError("histogram conversion must be positive"))
    return HistogramDensity(
        histogram.edges .* factor,
        histogram.density ./ factor
    )
end
