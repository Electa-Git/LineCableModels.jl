
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
	# This routes to the specific shape logic you write in solidcore.jl, etc.
	part = build_part(Target, Shape, b.grp, current_r, b.payload)
	return validate(part)
end

# ==========================================
# THE DSL HOOK
# ==========================================
@inline function Builder(
	::Type{Target}, ::Type{Shape}, grp::Symbol, args...,
) where {Target, Shape}

	# Wrap the types in Val{}() to make them 100% concrete for the tuple iteration
	grids = (Grid(Val{Target}()), Grid(Val{Shape}()), Grid(grp), map(Grid, args)...)

	return Gridspace{PartBuilder}(grids)
end
