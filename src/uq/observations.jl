"Return primitive results in Gridspace traversal order."
result(value::Union{LinearErrorResult, MonteCarloResult}) = value.values

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

@inline function _observe_product(product::RLCG, field::Symbol, indices...)
    return _product_value(getfield(product, field), indices)
end

observe(product::RLCG, ::typeof(R), indices...) =
    _observe_product(product, :R, indices...)
observe(product::RLCG, ::typeof(L), indices...) =
    _observe_product(product, :L, indices...)
observe(product::RLCG, ::typeof(C), indices...) =
    _observe_product(product, :C, indices...)
observe(product::RLCG, ::typeof(Engine.G), indices...) =
    _observe_product(product, :G, indices...)

@inline _statistic(transform, value::AbstractArray) = map(transform, value)
@inline _statistic(transform, value) = transform(value)

const _StatisticSelector = Union{
    typeof(Statistics.mean),
    typeof(Statistics.std),
    typeof(Statistics.median),
    typeof(minimum),
    typeof(maximum)
}

for selector in (R, L, C, Engine.G)
    @eval function observe(
            product::RLCG,
            ::typeof($selector),
            transform::_StatisticSelector,
            indices...
    )
        return _statistic(transform, observe(product, $selector, indices...))
    end
end

const _CableUQValue = Union{SampleSummary, HistogramDensity, AbstractArray}

@inline function _observe_cable_product(product, field::Symbol, indices...)
    return _product_value(getfield(product, field), indices)
end

for (selector, field) in ((R, :R), (L, :L), (C, :C))
    @eval function observe(
            product::DataModel.CableConstants{<:_CableUQValue},
            ::typeof($selector),
            indices...
    )
        return _observe_cable_product(product, $(QuoteNode(field)), indices...)
    end
    @eval function observe(
            product::DataModel.CableConstants{<:_CableUQValue},
            ::typeof($selector),
            transform::_StatisticSelector,
            indices...
    )
        return _statistic(transform, observe(product, $selector, indices...))
    end
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
