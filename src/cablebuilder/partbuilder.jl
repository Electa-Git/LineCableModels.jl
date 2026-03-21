# ==========================================
# THE VALIDATION BOUNDARY
# ==========================================
# If a specific part doesn't define topological rules, it passes.
@inline validate(part::AbstractCablePart) = part

# ==========================================
# THE ATOMIC BUILDER (Physics & Intrinsic Geometry)
# ==========================================
struct PartBuilder{Target, Shape, P <: Tuple}
	cmp::Symbol
	payload::P
end

# THE CONSTRUCTOR (The Zero-Alloc Val Interceptor)
@inline function PartBuilder(
	::Val{Target}, ::Val{Shape}, cmp::Symbol, payload...,
) where {Target, Shape}
	return PartBuilder{Target, Shape, typeof(payload)}(cmp, payload)
end

# THE FUNCTOR (The Spatial Collapse - Restored)
@inline function (b::PartBuilder{Target, Shape})(current_r) where {Target, Shape}
	# Pure coaxial materialization. Zero 2D awareness.
	part = build_part(Target, Shape, b.cmp, current_r, b.payload)
	return validate(part)
end

# ==========================================
# THE DSL HOOK
# ==========================================
@inline function Builder(
	::Type{Target}, ::Type{Shape}, cmp::Symbol, args...,
) where {Target, Shape}
	grids = (
		Grid(Val{Target}()),
		Grid(Val{Shape}()),
		Grid(cmp),
		map(Grid, args)...,
	)
	return Gridspace{PartBuilder}(grids)
end

# ==========================================
# THE GROUP BUILDER (Topology & Layout)
# ==========================================
struct GroupBuilder{P <: Tuple}
	payload::P
end

@inline function GroupBuilder(::Val{PartGroup}, payload...)
	return GroupBuilder{typeof(payload)}(payload)
end

@inline function (b::GroupBuilder)(current_r)
	origin = b.payload[1]
	n      = b.payload[2]
	m      = b.payload[3]

	inner_builders = Base.tail(Base.tail(Base.tail(b.payload)))

	# 1. Local coaxial stacking
	local_parts = build_design(0.0, inner_builders)
	local_r_ex = r_ex(local_parts[end])

	# 2. Global topological translation
	ox, oy = origin
	bound_r_in = current_r

	if n == 1 && m == 1
		bound_r_ex = max(current_r, sqrt(ox^2 + oy^2) + local_r_ex)
	else
		bound_r_ex = current_r + (m * 2 * local_r_ex)
	end

	T = promote_type(typeof(bound_r_in), typeof(bound_r_ex), typeof(ox), typeof(oy))

	# 3. Emit the anonymous topological folder
	part = PartGroup(
		convert(T, bound_r_in),
		convert(T, bound_r_ex),
		(convert(T, ox), convert(T, oy)),
		n, m,
		local_parts,
	)

	return validate(part)
end

@inline function Builder(
	::Type{PartGroup}, origin, n, m, layers::Tuple,
)
	grids = (
		Grid(Val{PartGroup}()),
		Grid((origin,)),
		Grid(n),
		Grid(m),
		layers...,
	)

	return Gridspace{GroupBuilder}(grids)
end
