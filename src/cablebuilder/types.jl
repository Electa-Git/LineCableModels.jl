abstract type AbstractShape{T <: Real} end
abstract type AbstractCablePart end

@inline r_ex(p::AbstractCablePart) = r_ex(p.shape)
@inline r_in(p::AbstractCablePart) = r_in(p.shape)

struct ConductorPart{T, S <: AbstractShape{T}} <: AbstractCablePart
	cmp::Symbol
	shape::S
	material::Material{T}
end

@inline function ConductorPart(
	cmp::Symbol,
	shape::AbstractShape{S},
	mat::Material{M},
) where {S <: Real, M <: Real}
	T = promote_type(S, M)
	s = convert(AbstractShape{T}, shape)
	m = convert(Material{T}, mat)

	return ConductorPart{T, typeof(s)}(cmp, s, m)
end

struct InsulatorPart{T, S <: AbstractShape{T}} <: AbstractCablePart
	cmp::Symbol
	shape::S
	material::Material{T}
end

@inline function InsulatorPart(
	cmp::Symbol,
	shape::AbstractShape{S},
	mat::Material{M},
) where {S <: Real, M <: Real}
	T = promote_type(S, M)
	s = convert(AbstractShape{T}, shape)
	m = convert(Material{T}, mat)

	return InsulatorPart{T, typeof(s)}(cmp, s, m)
end

# ==========================================
# THE TOPOLOGICAL VAULT
# ==========================================
struct PartGroup{T <: Real, P <: Tuple} <: AbstractCablePart
	r_in::T
	r_ex::T
	origin::Tuple{T, T}
	n::Int
	m::Int
	parts::P
end

@inline r_ex(g::PartGroup) = g.r_ex
@inline r_in(g::PartGroup) = g.r_in

# ---------------------------------------------------------
# THE GLOBAL RECAST FALLBACKS
# ---------------------------------------------------------
# 1. Reals get standard numeric conversion
@inline recast(::Type{T}, x::Real) where {T} = convert(T, x)

# 2. Everything else (Symbols, Bools, Strings) is ignored and passed through safely
@inline recast(::Type{T}, x) where {T} = x