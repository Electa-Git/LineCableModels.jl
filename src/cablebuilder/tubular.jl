# ==========================================
# THE VAULT
# ==========================================
@relax struct TubularLayer{T <: Real, P <: AbstractShapeParams{T}} <: AbstractShape{T}
	r_in::T
	r_ex::T
	params::P
end

# We use a Union to safely catch both Conductors and Insulators 
# without introducing dynamic `isa` checks or type-pirating the base AbstractCablePart.
@inline function validate(
	part::Union{ConductorPart{T, <:TubularLayer}, InsulatorPart{T, <:TubularLayer}},
) where {T}
	shape = part.shape

	# Cascade to intrinsic validation (checks t > 0)
	validate(shape.params)

	# Topological bounds check
	shape.r_ex > shape.r_in || throw(
		DomainError(
			shape.r_ex,
			"Physics violation: TubularLayer outer radius ($(shape.r_ex)) must be strictly greater than inner radius ($(shape.r_in)).",
		),
	)

	return part
end


# ==========================================
# THE FUNCTOR SPECIALIZATION (Spatial Collapse)
# ==========================================
@inline function build_part(
	::Type{Target},
	::Type{TubularLayer},
	cmp::Symbol,
	prev_bound::Circular{T},
	payload::Tuple{M, A},
) where {Target, T <: Real, M <: Material, A <: Annular}

	mat, params = payload

	# Conformal anchor: The tube strictly wraps the inner circular boundary.
	r_in = prev_bound.r

	# Extrusion: Expand by the intrinsic payload thickness.
	r_ex = r_in + params.t

	# Collapse the shape geometry
	shape = TubularLayer(r_in, r_ex, params)

	# Emit the atomic physics part (Zero 2D awareness here)
	return Target(cmp, shape, mat)
end