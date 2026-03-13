
# ==========================================
# THE VALIDATION BOUNDARY
# ==========================================
# Fallback: If a shape doesn't define rules, it passes. 
# (Or make it throw an error to force yourself to write rules for everything).
@inline validate(::Type{Shape}, payload::Tuple) where {Shape} = true

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
	return build_part(Target, Shape, b.grp, current_r, b.payload)
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
