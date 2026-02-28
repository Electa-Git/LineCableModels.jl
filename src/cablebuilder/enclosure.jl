struct Enclosure{L, T <: Real, S <: AbstractShape{L, T}} <: AbstractShape{L, T}
	base_shape::S
	filler_material::Material{T}
end

function Enclosure(
	base_shape::AbstractShape{L, T_shape},
	filler::Material{T_mat},
) where {L, T_shape <: Real, T_mat <: Real}
	# Find the mathematical truce between the shape and the fluid
	T_common = promote_type(T_shape, T_mat)

	# Force the inner shape and the material to upgrade if necessary
	s_promoted = convert(AbstractShape{L, T_common}, base_shape)
	f_promoted = convert(Material{T_common}, filler)

	# Shove the perfectly aligned structs into the auto-generated Vault
	return Enclosure{L, T_common, typeof(s_promoted)}(s_promoted, f_promoted)
end

function Base.convert(::Type{<:AbstractShape{L, T}}, e::Enclosure{L}) where {L, T <: Real}
	# Recursively upgrade the payload
	s_converted = convert(AbstractShape{L, T}, e.base_shape)
	f_converted = convert(Material{T}, e.filler_material)

	# Lock them in a new Vault
	return Enclosure{L, T, typeof(s_converted)}(s_converted, f_converted)
end

r_in(e::Enclosure) = r_in(e.base_shape)
r_ex(e::Enclosure) = r_ex(e.base_shape)
