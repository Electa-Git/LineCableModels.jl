struct Enclosure{T <: Real, S <: AbstractShape{T}} <: AbstractShape{T}
	base_shape::S
	filler_material::Material{T}
end

function Enclosure(
	base_shape::AbstractShape{T_shape},
	filler::Material{T_mat},
) where {T_shape <: Real, T_mat <: Real}

	T = promote_type(T_shape, T_mat)

	s = convert(AbstractShape{T}, base_shape)   # must return concrete
	f = convert(Material{T}, filler)

	return Enclosure{T, typeof(s)}(s, f)
end

function Base.convert(::Type{<:AbstractShape{T}}, e::Enclosure) where {T <: Real}
	# Recursively upgrade the payload
	s_converted = convert(AbstractShape{T}, e.base_shape)
	f_converted = convert(Material{T}, e.filler_material)

	# Lock them in a new Vault
	return Enclosure{T, typeof(s_converted)}(s_converted, f_converted)
end

# Override the global accessors because Enclosure is a diva
r_in(e::Enclosure) = r_in(e.base_shape)
r_ex(e::Enclosure) = r_ex(e.base_shape)

struct EnclosureBuilder{P, S, O, F}
	inner::S
	offset::O
	filler::Material{F}
end

@inline function EnclosureBuilder{P}(
	inner::S,
	offset::O,
	filler::Material{F},
) where {P, S, O, F}
	return EnclosureBuilder{P, S, O, F}(inner, offset, filler)
end

@inline function (b::EnclosureBuilder{P})(current_r::T) where {P, T <: Real}
	r0 = current_r + b.offset
	part = b.inner(r0)
	newshape = Enclosure(part.shape, b.filler)
	return P(part.grp, newshape, part.material)
end

struct EnclosureSpec{P, S, O, F} <: AbstractSpec{EnclosureBuilder{P}}
	inner::S
	offset::O
	filler::F
end

# # Make this diva explicit about what is iterable, and what is not. 
# @inline grid_args(spec::EnclosureSpec) = (spec.inner, spec.offset, spec.filler)

@inline function EnclosureSpec(::Type{P}, inner::S, offset::O, filler::F) where {P, S, O, F}
	return EnclosureSpec{P, S, O, F}(inner, offset, filler)
end

@inline function EnclosureSpec(::Type{P}, inner_spec, filler; offset = 0.0) where {P}
	filler_spec = convert(AbstractSpec{Material}, filler)
	return EnclosureSpec(P, inner_spec, Grid(offset), filler_spec)
end
