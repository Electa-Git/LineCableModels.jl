# ==========================================
# 1. THE VAULT
# ==========================================
# Just strictly holds the universal boundaries.
@relax struct SolidCore{T <: Real, P <: AbstractShapeParams{T}} <: AbstractShape{T}
	r_in::T
	r_ex::T
	params::P
end

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
	cmp::Symbol,
	current_r::T,
	payload::Tuple{M, C},
) where {Target, T <: Real, M <: Material, C <: Circular}

	mat, params = payload

	current_r <= eps(T) ||
		throw(DomainError(current_r, "Topological violation: SolidCore must be at r=0."))

	shape = SolidCore(current_r, params.r, params)

	return Target(cmp, shape, mat)
end

