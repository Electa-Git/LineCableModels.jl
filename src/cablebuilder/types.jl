abstract type AbstractCablePart end

struct ConductorPart{L, T, S <: AbstractShape{L, T}} <: AbstractCablePart
	tag::Symbol
	shape::S
	material::Material{T}
end

function ConductorPart(
	tag::Symbol,
	shape::AbstractShape{L, T_shape},
	mat::Material{T_mat},
) where {L, T_shape, T_mat}
	# 1. Ask the compiler what the safest common type is.
	T_common = promote_type(T_shape, T_mat)

	# 2. Let native dispatch handle the translation.
	shape_promoted = convert(AbstractShape{L, T_common}, shape)
	mat_promoted = convert(Material{T_common}, mat)

	return ConductorPart{L, T_common, typeof(shape_promoted)}(
		tag,
		shape_promoted,
		mat_promoted,
	)
end

struct InsulatorPart{L, T, S <: AbstractShape{L, T}} <: AbstractCablePart
	tag::Symbol
	shape::S
	material::Material{T}
end

function InsulatorPart(
	tag::Symbol,
	shape::AbstractShape{L, T_shape},
	mat::Material{T_mat},
) where {L, T_shape, T_mat}
	# 1. Ask the compiler what the safest common type is.
	T_common = promote_type(T_shape, T_mat)

	# 2. Let native dispatch handle the translation.
	shape_promoted = convert(AbstractShape{L, T_common}, shape)
	mat_promoted = convert(Material{T_common}, mat)

	return InsulatorPart{L, T_common, typeof(shape_promoted)}(
		tag,
		shape_promoted,
		mat_promoted,
	)
end