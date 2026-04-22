import Base: iterate, length, eltype, extrema
using Measurements

# ---------------------------------------------------------
# The Vaults (Strictly type-constrained to Tuples)
# ---------------------------------------------------------
struct DeterministicGrid{V <: Tuple}
	vals::V
end

struct RelativeGrid{V <: Tuple, P <: Tuple}
	vals::V
	rel_err::P
end

struct AbsoluteGrid{V <: Tuple, P <: Tuple}
	vals::V
	abs_err::P
end

# The explicit tag for absolute standard deviations. 
# If someone asks what this does, fire them.
struct AbsoluteError{T <: Tuple}
	vals::T
end

AbsoluteError(x::AbstractArray) = AbsoluteError(Tuple(x))
AbsoluteError(x::Real) = AbsoluteError((x,))


# ---------------------------------------------------------
# Surface API: The Formal Normalization Grammar
# Tuples and Arrays are collections. Everything else is a scalar.
# ---------------------------------------------------------

# --- 1. Deterministic ---
Grid(v::Tuple) = DeterministicGrid(v)
Grid(v::AbstractArray) = DeterministicGrid(Tuple(v))
Grid(v::Any) = DeterministicGrid((v,))

# --- 2. Relative (v, p) ---
Grid(v::Tuple, p::Tuple) = RelativeGrid(v, p)
Grid(v::Tuple, p::AbstractArray) = RelativeGrid(v, Tuple(p))
Grid(v::Tuple, p::Any) = RelativeGrid(v, (p,))

Grid(v::AbstractArray, p::Tuple) = RelativeGrid(Tuple(v), p)
Grid(v::AbstractArray, p::AbstractArray) = RelativeGrid(Tuple(v), Tuple(p))
Grid(v::AbstractArray, p::Any) = RelativeGrid(Tuple(v), (p,))

Grid(v::Any, p::Tuple) = RelativeGrid((v,), p)
Grid(v::Any, p::AbstractArray) = RelativeGrid((v,), Tuple(p))
Grid(v::Any, p::Any) = RelativeGrid((v,), (p,))

# --- 3. Absolute (v, a) ---
Grid(v::Tuple, a::AbsoluteError) = AbsoluteGrid(v, a.vals)
Grid(v::AbstractArray, a::AbsoluteError) = AbsoluteGrid(Tuple(v), a.vals)
Grid(v::Any, a::AbsoluteError) = AbsoluteGrid((v,), a.vals)

# Pass-through for already built vaults
Grid(g::Union{DeterministicGrid, RelativeGrid, AbsoluteGrid}) = g

# ---------------------------------------------------------
# The Iteration Protocol (Measurements Gangbang)
# ---------------------------------------------------------

# Deterministic
@inline Base.iterate(g::DeterministicGrid, state...) = iterate(g.vals, state...)
@inline Base.length(g::DeterministicGrid) = length(g.vals)
Base.eltype(::Type{<:DeterministicGrid{V}}) where {V} = eltype(V)

# Relative
@inline function Base.iterate(g::RelativeGrid, state...)
	res = iterate(Iterators.product(g.vals, g.rel_err), state...)
	res === nothing && return nothing
	((v, p), next_state) = res
	return measurement(v, abs(v) * (p / 100.0)), next_state
end
@inline Base.length(g::RelativeGrid) = length(g.vals) * length(g.rel_err)
Base.eltype(::Type{<:RelativeGrid{V, P}}) where {V, P} =
	Measurement{promote_type(eltype(V), eltype(P))}

# Absolute
@inline function Base.iterate(g::AbsoluteGrid, state...)
	res = iterate(Iterators.product(g.vals, g.abs_err), state...)
	res === nothing && return nothing
	((v, err), next_state) = res
	return measurement(v, abs(err)), next_state
end
@inline Base.length(g::AbsoluteGrid) = length(g.vals) * length(g.abs_err)
Base.eltype(::Type{<:AbsoluteGrid{V, P}}) where {V, P} =
	Measurement{promote_type(eltype(V), eltype(P))}

# ---------------------------------------------------------
# Boundaries (For Lemonparty Solvers)
# ---------------------------------------------------------
@inline Base.extrema(g::DeterministicGrid) = (minimum(g.vals), maximum(g.vals))

@inline function Base.extrema(g::RelativeGrid)
	v_min, v_max = minimum(g.vals), maximum(g.vals)
	p_max = maximum(abs, g.rel_err) / 100.0
	return (v_min * (1.0 - p_max), v_max * (1.0 + p_max))
end

@inline function Base.extrema(g::AbsoluteGrid)
	v_min, v_max = minimum(g.vals), maximum(g.vals)
	err_max = maximum(abs, g.abs_err)
	return (v_min - err_max, v_max + err_max)
end

# ---------------------------------------------------------
# The Stochastic Sampler
# ---------------------------------------------------------
import Base: rand
using Distributions

@inline Base.rand(g::DeterministicGrid, ::Type{D}) where {D} = rand(g.vals)

@inline function Base.rand(g::RelativeGrid, ::Type{D}) where {D}
	v, p = rand(g.vals), rand(g.rel_err)
	σ = abs(v) * (p / 100.0)
	σ == 0 && return float(v)
	return D <: Normal ? rand(Normal(v, σ)) : rand(Uniform(v - √3*σ, v + √3*σ))
end

@inline function Base.rand(g::AbsoluteGrid, ::Type{D}) where {D}
	v, σ = rand(g.vals), rand(g.abs_err)
	σ == 0 && return float(v)
	return D <: Normal ? rand(Normal(v, σ)) : rand(Uniform(v - √3*σ, v + √3*σ))
end

@inline Base.rand(
	g::Union{DeterministicGrid, RelativeGrid, AbsoluteGrid};
	dist::Type{D} = Normal,
) where {D} = rand(g, D)
