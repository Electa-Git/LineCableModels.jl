# ==========================================
# 1. THE VAULT
# ==========================================
# Just strictly holds the universal boundaries.
@relax struct SolidCore{T <: Real, P <: AbstractShapeParams{T}} <: AbstractShape{T}
	r_in::T
	r_ex::T
	params::P
end

# # Outer constructor: Extract T3 directly from the abstract signature and promote
# @inline function SolidCore(
# 	r_in::T1,
# 	r_ex::T2,
# 	params::AbstractShapeParams{T3},
# ) where {T1 <: Real, T2 <: Real, T3 <: Real}

# 	T = promote_type(T1, T2, T3)

# 	# Idiomatic Julia conversion. Routes to the 1-liners in primitives.jl
# 	p_cast = convert(AbstractShapeParams{T}, params)

# 	return SolidCore{T, typeof(p_cast)}(convert(T, r_in), convert(T, r_ex), p_cast)
# end

# # Convert method: Propagates conversion cleanly down the tree
# @inline function Base.convert(
# 	::Type{<:AbstractShape{T}},
# 	s::SolidCore,
# ) where {T <: Real}

# 	p_cast = convert(AbstractShapeParams{T}, s.params)
# 	return SolidCore{T, typeof(p_cast)}(convert(T, s.r_in), convert(T, s.r_ex), p_cast)
# end

@inline function validate(part::ConductorPart{T, <:SolidCore}) where {T}
	shape = part.shape

	# Cascade to intrinsic validation
	validate(shape.params)

	# Topological bounds checks
	shape.r_in == zero(T) || throw(
		DomainError(
			shape.r_in, "Topological violation: SolidCore MUST start exactly at r=0.",
		),
	)

	shape.r_ex > shape.r_in || throw(
		DomainError(
			shape.r_ex,
			"Physics violation: Outer radius ($(shape.r_ex)) must be > inner radius ($(shape.r_in)).",
		),
	)

	return part
end

# ==========================================
# 2. THE FUNCTOR SPECIALIZATION (Spatial Collapse)
# ==========================================
@inline function build_part(
	::Type{Target},
	::Type{SolidCore},
	grp::Symbol,
	current_r::T,
	payload::Tuple{M, C},
) where {Target, T <: Real, M <: Material, C <: Circular}

	mat, params = payload

	current_r <= eps(T) || throw(DomainError(
		current_r, "Topological violation: SolidCore must be at r=0.",
	))

	shape = SolidCore(current_r, params.r, params)

	return Target(grp, shape, mat)
end
