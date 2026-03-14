# ==============================================================================
# THE GRIDSPACE ENGINE (Replaces AbstractSpec)
# ==============================================================================

# The Universal Staging Area. 
# Target is the strict `<: Real` physics struct. Args is the tuple of Grids/Scalars.
struct Gridspace{Target, Args <: Tuple}
	grids::Args
end

Gridspace{Target}(grids::Args) where {Target, Args <: Tuple} =
	Gridspace{Target, Args}(grids)

Grid(g::Gridspace) = g

# ==============================================================================
# THE THIN-ALLOCATION ITERATOR PROTOCOL
# ==============================================================================

# 1. The Initializer
@inline function Base.iterate(g::Gridspace{Target}) where {Target}
	# Iterators.product is natively stack-allocated. No closures.
	iter = Iterators.product(g.grids...)
	next = iterate(iter)

	next === nothing && return nothing

	args, state = next
	# Native splat into Target. Completely type-stable.
	return Target(args...), state
end

# 2. The Advancer
@inline function Base.iterate(g::Gridspace{Target}, state) where {Target}
	iter = Iterators.product(g.grids...)
	next = iterate(iter, state)

	next === nothing && return nothing

	args, new_state = next
	return Target(args...), new_state
end

# 3. Utilities (So the compiler knows exactly how big the loop is)
Base.IteratorSize(::Type{<:Gridspace}) = Base.HasShape{1}()
Base.length(g::Gridspace) = prod(length, g.grids)
Base.size(g::Gridspace) = (length(g),)

# ---------------------------------------------------------
# The Stochastic Sampler
# ---------------------------------------------------------
import Base: rand
using Distributions

@inline function Base.rand(
	g::Gridspace{Target},
	::Type{D},
) where {Target, D <: ContinuousUnivariateDistribution}
	# Map distributes your existing grid.jl rand() over the tuple
	samples = map(grid -> rand(grid, D), grid_args(g))
	return Target(samples...)
end

Base.rand(
	g::Gridspace;
	dist::Type{D} = Normal,
) where {D <: ContinuousUnivariateDistribution} = rand(g, dist)