Base.size(value::LinearErrorResult) = size(value.values)
Base.length(value::LinearErrorResult) = length(value.values)
Base.getindex(value::LinearErrorResult, index::Integer) = value.values[index]
Base.iterate(value::LinearErrorResult, state...) = iterate(value.values, state...)
Base.IndexStyle(::Type{<:LinearErrorResult}) = IndexLinear()

Base.size(value::MonteCarloResult) = size(value.values)
Base.length(value::MonteCarloResult) = length(value.values)
Base.getindex(value::MonteCarloResult, index::Integer) = value.values[index]
Base.iterate(value::MonteCarloResult, state...) = iterate(value.values, state...)
Base.IndexStyle(::Type{<:MonteCarloResult}) = IndexLinear()
