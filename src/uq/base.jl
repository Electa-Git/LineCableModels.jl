Base.length(value::LinearErrorResult) = length(value.values)
Base.size(value::LinearErrorResult) = (length(value),)
Base.getindex(value::LinearErrorResult, index::Integer) = value.values[index]
Base.iterate(value::LinearErrorResult, state...) = iterate(value.values, state...)
Base.firstindex(value::LinearErrorResult) = firstindex(value.values)
Base.lastindex(value::LinearErrorResult) = lastindex(value.values)

Base.length(value::MonteCarloResult) = length(value.values)
Base.size(value::MonteCarloResult) = (length(value),)
Base.getindex(value::MonteCarloResult, index::Integer) = value.values[index]
Base.iterate(value::MonteCarloResult, state...) = iterate(value.values, state...)
Base.firstindex(value::MonteCarloResult) = firstindex(value.values)
Base.lastindex(value::MonteCarloResult) = lastindex(value.values)
