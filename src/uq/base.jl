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

#! explicit-imports: off
# Base's iterator trait protocol exposes these values without public bindings.
Base.IteratorSize(::Type{<:SensitivityResult}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:SensitivityResult}) = Base.HasEltype()
#! explicit-imports: on
Base.eltype(::Type{<:SensitivityResult{<:Any, V}}) where {V} = eltype(V)
Base.length(value::SensitivityResult) = length(value.values)
Base.size(value::SensitivityResult) = (length(value),)
Base.getindex(value::SensitivityResult, index::Integer) = value.values[index]
Base.iterate(value::SensitivityResult, state...) = iterate(value.values, state...)
Base.firstindex(value::SensitivityResult) = firstindex(value.values)
Base.lastindex(value::SensitivityResult) = lastindex(value.values)
