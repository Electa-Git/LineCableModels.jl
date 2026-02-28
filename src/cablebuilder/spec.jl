# The Target parameter tells the compiler what strict physics struct to spit out.
abstract type AbstractSpec{Target} end

# The Introspection Kernel. 
# ntuple + fieldcount forces the compiler to unroll the fields at compile-time. 
# This preserves 100% type stability across mixed grid types.
# Wrapping in Val() forces the compiler to unroll this at compile-time.
@inline grid_args(spec::AbstractSpec) =
	ntuple(i -> getfield(spec, i), Val(fieldcount(typeof(spec))))

# The Generic Generator.
# Splats the unrolled tuple of grids into Iterators.product, 
# then splats the resulting combination directly into the Target constructor.
@inline generator(spec::AbstractSpec{T}) where {T} = (
	T(args...) for args in Iterators.product(grid_args(spec)...)
)

# ---------------------------------------------------------
# The Generic Iteration Protocol
# ---------------------------------------------------------
Base.iterate(spec::AbstractSpec) = iterate(generator(spec))
Base.iterate(spec::AbstractSpec, state) = iterate(generator(spec), state)

# Calculates the total combinatorial gangbang size natively.
# init=1 prevents it from shitting the bed if you ever define a struct with zero fields.
Base.length(spec::AbstractSpec) =
	prod((length(g) for g in grid_args(spec)), init = 1)

# Let Julia's native Generator interface infer the eltype dynamically
# Base.IteratorEltype(::Type{<:AbstractSpec}) = Base.EltypeUnknown()
Base.IteratorEltype(::Type{<:AbstractSpec}) = Base.HasEltype()
Base.eltype(::Type{<:AbstractSpec{Target}}) where {Target} = Target


# The Generic Stochastic Realizer
import Base: rand
using Distributions
@inline grid_args(spec::AbstractSpec) =
	ntuple(i -> getfield(spec, i), Val(fieldcount(typeof(spec))))

@inline function Base.rand(
	spec::AbstractSpec{Target};
	dist::Type{D} = Normal,
) where {Target, D <: ContinuousUnivariateDistribution}
	# Pass the Type positionally. No kwargs in the closure.
	samples = map(g -> rand(g, D), grid_args(spec))

	return Target(samples...)
end