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

# No more ntuple/fieldcount bullshit. We just return the tuple.
@inline grid_args(g::Gridspace) = g.grids

# The Generator. 
# Splats the unrolled iterators directly into the Target constructor.
@inline generator(g::Gridspace{Target}) where {Target} = (
	Target(args...) for args in Iterators.product(grid_args(g)...)
)

# ---------------------------------------------------------
# The Iteration Protocol
# ---------------------------------------------------------
@inline function Base.iterate(g::Gridspace)
	gen = generator(g)
	y = iterate(gen)
	y === nothing && return nothing
	return (y[1], (gen, y[2]))
end

@inline function Base.iterate(g::Gridspace, state)
	gen, st = state
	y = iterate(gen, st)
	y === nothing && return nothing
	return (y[1], (gen, y[2]))
end

Base.length(g::Gridspace) = prod((length(x) for x in grid_args(g)), init = 1)
Base.IteratorSize(::Type{<:Gridspace}) = Base.HasLength()
Base.IteratorEltype(::Type{<:Gridspace}) = Base.EltypeUnknown()

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