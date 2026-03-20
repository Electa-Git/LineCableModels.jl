# Shape params define the shape of a primitive, but not its material, group or location/layout.
abstract type AbstractShapeParams{T <: Real} end

@gridspace @relax struct Circular{T <: Real} <: AbstractShapeParams{T}
	r::T
end

@inline function validate(p::Circular)
	p.r > zero(p.r) || throw(DomainError(p.r, "Circular radius must be strictly positive."))
	return p
end

@gridspace @relax struct Rectangular{T <: Real} <: AbstractShapeParams{T}
	w::T
	h::T
end

@inline function validate(p::Rectangular)
	p.w > zero(p.w) ||
		throw(DomainError(p.w, "Rectangular width must be strictly positive."))
	p.h > zero(p.h) ||
		throw(DomainError(p.h, "Rectangular height must be strictly positive."))
	return p
end

# If a specific payload vault doesn't define intrinsic rules, it passes.
@inline validate(params::AbstractShapeParams) = params