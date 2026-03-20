# ==========================================
# THE VALIDATION BOUNDARY
# ==========================================
# If a specific part doesn't define topological rules, it passes.
@inline validate(part::AbstractCablePart) = part

# ==========================================
# THE BUILDER
# ==========================================
struct PartBuilder{Target, Shape, P <: Tuple}
	grp::Symbol
	payload::P
end

# 1. THE CONSTRUCTOR (The Zero-Alloc Val Interceptor)
@inline function PartBuilder(
	::Val{Target}, ::Val{Shape}, grp::Symbol, payload...,
) where {Target, Shape}

	return PartBuilder{Target, Shape, typeof(payload)}(grp, payload)
end

# 2. THE FUNCTOR (The Spatial Collapse - Restored)
@inline function (b::PartBuilder{Target, Shape})(current_r) where {Target, Shape}
	origin = b.payload[1]
	shape_args = Base.tail(b.payload)

	# Pass origin down to the concrete shape builder
	part = build_part(Target, Shape, b.grp, origin, current_r, shape_args)
	return validate(part)
end

# ==========================================
# THE DSL HOOK
# ==========================================
@inline function Builder(
	::Type{Target}, ::Type{Shape}, grp::Symbol, args...,
) where {Target, Shape}

	# Slot 4 is always the origin. Default to concentric.
	grids = (
		Grid(Val{Target}()),
		Grid(Val{Shape}()),
		Grid(grp),
		Grid(((0.0, 0.0),)),
		map(Grid, args)...,
	)

	return Gridspace{PartBuilder}(grids)
end
