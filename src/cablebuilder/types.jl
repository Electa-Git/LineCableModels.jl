abstract type AbstractCablePart end

@inline r_ex(p::AbstractCablePart) = r_ex(p.shape)

struct ConductorPart{L, T, S <: AbstractShape{L, T}} <: AbstractCablePart
	grp::Symbol
	shape::S
	material::Material{T}
end

function ConductorPart(
	grp::Symbol,
	shape::AbstractShape{L, T_shape},
	mat::Material{T_mat},
) where {L, T_shape <: Real, T_mat <: Real}
	# 1. Ask the compiler what the safest common type is.
	T = promote_type(T_shape, T_mat)

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

function InsulatorPart(
	grp::Symbol,
	shape::AbstractShape{L, T_shape},
	mat::Material{T_mat},
) where {L, T_shape <: Real, T_mat <: Real}
	# 1. Ask the compiler what the safest common type is.
	T = promote_type(T_shape, T_mat)

	# 2. Let native dispatch handle the translation.
	s = convert(AbstractShape{L, T}, shape)
	m = convert(Material{T}, mat)

	return InsulatorPart{L, T, typeof(s)}(grp, s, m)
end