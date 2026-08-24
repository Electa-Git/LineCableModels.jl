result(value::Union{LinearErrorResult, MonteCarloResult}) = value.values
statistics(value::MonteCarloResult) = value.stats
samples(value::MonteCarloResult) = value.sample_values
histograms(value::MonteCarloResult) = value.histogram_values
uncertain_value(value::LinearErrorResult) = value.values
observables(value::LinearErrorResult) = (result = result(value),)

function observables(
        value::MonteCarloResult{<:Any, <:Any, <:AbstractVector, Nothing, Nothing},
)
    return (result = result(value), statistics = statistics(value))
end
function observables(
        value::MonteCarloResult{<:Any, <:Any, <:AbstractVector, <:Any, Nothing},
)
    return (
        result = result(value),
        statistics = statistics(value),
        samples = samples(value)
    )
end
function observables(
        value::MonteCarloResult{<:Any, <:Any, <:AbstractVector, Nothing, <:Any},
)
    return (
        result = result(value),
        statistics = statistics(value),
        histograms = histograms(value)
    )
end
function observables(value::MonteCarloResult)
    return (
        result = result(value),
        statistics = statistics(value),
        samples = samples(value),
        histograms = histograms(value)
    )
end
