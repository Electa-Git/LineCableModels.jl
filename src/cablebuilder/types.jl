abstract type AbstractLayout end
struct Concentric <: AbstractLayout end
struct SectorShaped <: AbstractLayout end
abstract type AbstractShape{L <: AbstractLayout, T <: Real} end

abstract type AbstractCablePart end

@inline r_ex(p::AbstractCablePart) = r_ex(p.shape)

struct ConductorPart{L, T, S <: AbstractShape{L, T}} <: AbstractCablePart
	grp::Symbol
	shape::S
	material::Material{T}
end

@inline function ConductorPart(
	grp::Symbol,
	shape::AbstractShape{L, S},
	mat::Material{M},
) where {L, S <: Real, M <: Real}
	# 1. Ask the compiler what the safest common type is.
	T = promote_type(S, M)

	# 2. Let native dispatch handle the translation.
	s = convert(AbstractShape{L, T}, shape)
	m = convert(Material{T}, mat)

	return ConductorPart{L, T, typeof(s)}(grp, s, m)
end

struct InsulatorPart{L, T, S <: AbstractShape{L, T}} <: AbstractCablePart
	grp::Symbol
	shape::S
	material::Material{T}
end

@inline function InsulatorPart(
	grp::Symbol,
	shape::AbstractShape{L, S},
	mat::Material{M},
) where {L, S <: Real, M <: Real}
	# 1. Ask the compiler what the safest common type is.
	T = promote_type(S, M)

	# 2. Let native dispatch handle the translation.
	s = convert(AbstractShape{L, T}, shape)
	m = convert(Material{T}, mat)

	return InsulatorPart{L, T, typeof(s)}(grp, s, m)
end
