abstract type AbstractShape{T <: Real} end

abstract type AbstractCablePart end

@inline r_ex(p::AbstractCablePart) = r_ex(p.shape)

struct ConductorPart{T, S <: AbstractShape{T}} <: AbstractCablePart
	grp::Symbol
	shape::S
	material::Material{T}
end

@inline function ConductorPart(
	grp::Symbol,
	shape::AbstractShape{S},
	mat::Material{M},
) where {S <: Real, M <: Real}
	# 1. Ask the compiler what the safest common type is.
	T = promote_type(S, M)

	# 2. Let native dispatch handle the translation.
	s = convert(AbstractShape{T}, shape)
	m = convert(Material{T}, mat)

	return ConductorPart{T, typeof(s)}(grp, s, m)
end

struct InsulatorPart{T, S <: AbstractShape{T}} <: AbstractCablePart
	grp::Symbol
	shape::S
	material::Material{T}
end

@inline function InsulatorPart(
	grp::Symbol,
	shape::AbstractShape{S},
	mat::Material{M},
) where {S <: Real, M <: Real}
	# 1. Ask the compiler what the safest common type is.
	T = promote_type(S, M)

	# 2. Let native dispatch handle the translation.
	s = convert(AbstractShape{T}, shape)
	m = convert(Material{T}, mat)

	return InsulatorPart{T, typeof(s)}(grp, s, m)
end

# ---------------------------------------------------------
# THE GLOBAL RECAST FALLBACKS
# ---------------------------------------------------------
# 1. Reals get standard numeric conversion
@inline recast(::Type{T}, x::Real) where {T} = convert(T, x)

# 2. Everything else (Symbols, Bools, Strings) is ignored and passed through safely
@inline recast(::Type{T}, x) where {T} = x