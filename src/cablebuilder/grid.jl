import Base: iterate, length, eltype, extrema
using Measurements

# ---------------------------------------------------------
# The Vaults (Do not expose these to the smooth-brains)
# ---------------------------------------------------------
struct DeterministicGrid{V}
	vals::V
end

struct RelativeGrid{V, P}
	vals::V
	pcts::P
end

struct AbsoluteGrid{V, P}
	vals::V
	abs_err::P
end

# ---------------------------------------------------------
# The Bouncers & Tags
# ---------------------------------------------------------
# If it's already a collection, pass it through
_tuplify(x::AbstractArray) = x
_tuplify(x::Tuple) = x

# If it's a scalar (Real, Symbol, Type, Material, etc.), wrap it in a 1-element Tuple so it can iterate
_tuplify(x) = (x,)

# The explicit tag for absolute standard deviations. 
# If someone asks what this does, fire them.
struct AbsoluteError{T}
	vals::T
end
AbsoluteError(x) = AbsoluteError{typeof(_tuplify(x))}(_tuplify(x))

# ---------------------------------------------------------
# Surface API 
# ---------------------------------------------------------
Grid(v) = DeterministicGrid(_tuplify(v))
Grid(v, p) = RelativeGrid(_tuplify(v), _tuplify(p))
Grid(v, a::AbsoluteError) = AbsoluteGrid(_tuplify(v), a.vals)

# If it's already a Grid vault, pass it through untouched.
Grid(g::DeterministicGrid) = g
Grid(g::RelativeGrid) = g
Grid(g::AbsoluteGrid) = g

# ---------------------------------------------------------
# The Iteration Protocol (Measurements Gangbang)
# ---------------------------------------------------------

# Deterministic
Base.iterate(g::DeterministicGrid) = iterate(g.vals)
Base.iterate(g::DeterministicGrid, state) = iterate(g.vals, state)
Base.length(g::DeterministicGrid) = length(g.vals)
Base.eltype(::Type{DeterministicGrid{V}}) where {V} = eltype(V)

# Relative
@inline _measurify_rel(::Nothing) = nothing
@inline _measurify_rel(res::Tuple) = (
	measurement(res[1][1], abs(res[1][1]) * (res[1][2] / 100.0)),
	res[2],
)
Base.iterate(g::RelativeGrid) = _measurify_rel(iterate(Iterators.product(g.vals, g.pcts)))
Base.iterate(g::RelativeGrid, state) =
	_measurify_rel(iterate(Iterators.product(g.vals, g.pcts), state))
Base.length(g::RelativeGrid) = length(g.vals) * length(g.pcts)
Base.eltype(::Type{<:RelativeGrid{V, P}}) where {V, P} =
	Measurement{promote_type(eltype(V), eltype(P))}

# Absolute
@inline _measurify_abs(::Nothing) = nothing
@inline _measurify_abs(res::Tuple) = (
	measurement(res[1][1], abs(res[1][2])),
	res[2],
)
Base.iterate(g::AbsoluteGrid) =
	_measurify_abs(iterate(Iterators.product(g.vals, g.abs_err)))
Base.iterate(g::AbsoluteGrid, state) =
	_measurify_abs(iterate(Iterators.product(g.vals, g.abs_err), state))
Base.length(g::AbsoluteGrid) = length(g.vals) * length(g.abs_err)
Base.eltype(::Type{<:AbsoluteGrid{V, P}}) where {V, P} =
	Measurement{promote_type(eltype(V), eltype(P))}

# ---------------------------------------------------------
# Boundaries (For Lemonparty Solvers)
# ---------------------------------------------------------
Base.extrema(g::DeterministicGrid) = (minimum(g.vals), maximum(g.vals))

function Base.extrema(g::RelativeGrid)
	v_min, v_max = minimum(g.vals), maximum(g.vals)
	p_max = maximum(abs, g.pcts) / 100.0
	return (v_min * (1.0 - p_max), v_max * (1.0 + p_max))
end

function Base.extrema(g::AbsoluteGrid)
	v_min, v_max = minimum(g.vals), maximum(g.vals)
	err_max = maximum(abs, g.abs_err)
	return (v_min - err_max, v_max + err_max)
end

# ---------------------------------------------------------
# The Stochastic Sampler
# ---------------------------------------------------------
import Base: rand
using Distributions

# 1. Positional Fast-Paths (Zero kwargs)
@inline Base.rand(g::DeterministicGrid, ::Type{D}) where {D} = rand(g.vals)

@inline function Base.rand(g::RelativeGrid, ::Type{D}) where {D}
	v, p = rand(g.vals), rand(g.pcts)
	σ = abs(v) * (p / 100.0)
	# The compiler deletes the unused branch at compile-time because D is known
	return D === Normal ? rand(Normal(v, σ)) : rand(Uniform(v - √3*σ, v + √3*σ))
end

@inline function Base.rand(g::AbsoluteGrid, ::Type{D}) where {D}
	v, σ = rand(g.vals), rand(g.abs_err)
	return D === Normal ? rand(Normal(v, σ)) : rand(Uniform(v - √3*σ, v + √3*σ))
end

# 2. Public API (Forwards kwargs to the positional fast-path)
Base.rand(
	g::Union{DeterministicGrid, RelativeGrid, AbsoluteGrid};
	dist::Type{D} = Normal,
) where {D} = rand(g, D)
