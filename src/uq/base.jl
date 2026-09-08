Base.length(value::LinearErrorResult) = length(value.values)
Base.size(value::LinearErrorResult) = (length(value),)
Base.getindex(value::LinearErrorResult, index::Integer) = value.values[index]
Base.iterate(value::LinearErrorResult, state...) = iterate(value.values, state...)
Base.firstindex(value::LinearErrorResult) = firstindex(value.values)
Base.lastindex(value::LinearErrorResult) = lastindex(value.values)

Base.length(value::PolynomialChaosResult) = length(value.values)
Base.size(value::PolynomialChaosResult) = (length(value),)
Base.getindex(value::PolynomialChaosResult, index::Integer) = value.values[index]
Base.iterate(value::PolynomialChaosResult, state...) = iterate(value.values, state...)
Base.firstindex(value::PolynomialChaosResult) = firstindex(value.values)
Base.lastindex(value::PolynomialChaosResult) = lastindex(value.values)

Base.length(value::MonteCarloResult) = length(value.values)
Base.size(value::MonteCarloResult) = (length(value),)
Base.getindex(value::MonteCarloResult, index::Integer) = value.values[index]
Base.iterate(value::MonteCarloResult, state...) = iterate(value.values, state...)
Base.firstindex(value::MonteCarloResult) = firstindex(value.values)
Base.lastindex(value::MonteCarloResult) = lastindex(value.values)
