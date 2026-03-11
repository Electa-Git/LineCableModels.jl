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
# The Generic Iteration Protocol (allocation-thin)
# ---------------------------------------------------------
@inline function Base.iterate(spec::AbstractSpec)
	gen = generator(spec)
	y = iterate(gen)
	y === nothing && return nothing
	return (y[1], (gen, y[2]))
end

@inline function Base.iterate(spec::AbstractSpec, state)
	gen, st = state
	y = iterate(gen, st)
	y === nothing && return nothing
	return (y[1], (gen, y[2]))
end

# Calculates the total combinatorial gangbang size natively.
# init=1 prevents it from shitting the bed if you ever define a struct with zero fields.
Base.length(spec::AbstractSpec) = prod((length(g) for g in grid_args(spec)), init = 1)
Base.IteratorSize(::Type{<:AbstractSpec}) = Base.HasLength()

# Forces Julia to infer the state tuple bottom-up from the materialized objects
Base.IteratorEltype(::Type{<:AbstractSpec}) = Base.EltypeUnknown()


# The Generic Stochastic Realizer
import Base: rand
using Distributions

# The positional fast-path for nested specs
@inline function Base.rand(
	spec::AbstractSpec{Target},
	::Type{D},
) where {Target, D <: ContinuousUnivariateDistribution}
	samples = map(g -> rand(g, D), grid_args(spec))
	return Target(samples...)
end

# The public API (kwargs)
Base.rand(
	spec::AbstractSpec;
	dist::Type{D} = Normal,
) where {D <: ContinuousUnivariateDistribution} = rand(spec, dist)
