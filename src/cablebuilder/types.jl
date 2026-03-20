abstract type AbstractShape{T <: Real} end

abstract type AbstractCablePart end

@inline r_ex(p::AbstractCablePart) = r_ex(p.shape)

struct ConductorPart{T, S <: AbstractShape{T}} <: AbstractCablePart
	grp::Symbol
	origin::Tuple{T, T}
	shape::S
	material::Material{T}
end

@inline function ConductorPart(
	grp::Symbol,
	origin::Tuple{O1, O2},
	shape::AbstractShape{S},
	mat::Material{M},
) where {O1 <: Real, O2 <: Real, S <: Real, M <: Real}
	# Ask the compiler what the safest common type is.
	T = promote_type(O1, O2, S, M)

	# Let native dispatch handle the translation.
	o = (convert(T, origin[1]), convert(T, origin[2]))
	s = convert(AbstractShape{T}, shape)
	m = convert(Material{T}, mat)

	return ConductorPart{T, typeof(s)}(grp, o, s, m)
end

struct InsulatorPart{T, S <: AbstractShape{T}} <: AbstractCablePart
	grp::Symbol
	origin::Tuple{T, T}
	shape::S
	material::Material{T}
end

@inline function InsulatorPart(
	grp::Symbol,
	origin::Tuple{O1, O2},
	shape::AbstractShape{S},
	mat::Material{M},
) where {O1 <: Real, O2 <: Real, S <: Real, M <: Real}
	# Ask the compiler what the safest common type is.
	T = promote_type(O1, O2, S, M)

	# Let native dispatch handle the translation.
	o = (convert(T, origin[1]), convert(T, origin[2]))
	s = convert(AbstractShape{T}, shape)
	m = convert(Material{T}, mat)

	return InsulatorPart{T, typeof(s)}(grp, o, s, m)
end

# ---------------------------------------------------------
# THE GLOBAL RECAST FALLBACKS
# ---------------------------------------------------------
# 1. Reals get standard numeric conversion
@inline recast(::Type{T}, x::Real) where {T} = convert(T, x)

# 2. Everything else (Symbols, Bools, Strings) is ignored and passed through safely
@inline recast(::Type{T}, x) where {T} = x