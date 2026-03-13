
# ==========================================
# THE VALIDATION BOUNDARY
# ==========================================
# Fallback: If a shape doesn't define rules, it passes. 
# (Or make it throw an error to force yourself to write rules for everything).
@inline validate(::Type{Shape}, payload::Tuple) where {Shape} = true

# ==========================================
# THE BUILDER
# ==========================================
struct PartBuilder{Part, Shape, P <: Tuple}
	part::Type{Part}
	shape::Type{Shape}
	grp::Symbol
	payload::P
end

# The Strict Constructor (The Sanity Checkpoint)
@inline function PartBuilder(
	part::Type{Part}, shape::Type{Shape}, grp::Symbol, payload...,
) where {Part, Shape}

	# THE CHOKEPOINT: Validates the payload before the object is built.
	validate(shape, payload)

	return PartBuilder{Part, Shape, typeof(payload)}(part, shape, grp, payload)
end

# The DSL Hook (Remains unchanged)
@inline function Builder(
	part::Type{Part}, shape::Type{Shape}, grp::Symbol, args...,
) where {Part, Shape}
	grids = (Grid(part), Grid(shape), Grid(grp), map(Grid, args)...)
	return Gridspace{PartBuilder}(grids)
end
