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
	grp::Symbol,
	prev_bound::Circular{T},  # <-- Dispatches on the primitive
	payload::Tuple{M, C},
) where {Target, T <: Real, M <: Material, C <: Circular}

	mat, params = payload

	prev_bound.r <= eps(T) ||
		throw(
			DomainError(
				prev_bound.r,
				"Topological violation: SolidCore must start at r=0.",
			),
		)

	shape = SolidCore(prev_bound.r, params.r, params)
	return Target(grp, shape, mat)
end

