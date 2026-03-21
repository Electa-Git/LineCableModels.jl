# ---------------------------------------------------------
# The Concrete Target 
# ---------------------------------------------------------
struct CableDesign{T <: Tuple}
	payload::T
end

# ---------------------------------------------------------
# The Allocation-Free Stacking Engine
# ---------------------------------------------------------
@inline build_design(r, ::Tuple{}) = ()

@inline function build_design(r, builders::Tuple)
	b = first(builders)
	target = b(r)

	# The part knows its own boundary envelope.
	# Atomic parts return their 1D coaxial radius.
	# PartGroups return their translated 2D boundary.
	part_boundary = r_ex(target)

	# The running radius must engulf the boundary
	next_r = max(r, part_boundary)

	return (target, build_design(next_r, Base.tail(builders))...)
end

# ---------------------------------------------------------
# The Constructor (Hit by the Gridspace Generator)
# ---------------------------------------------------------
# The generator splats the unrolled PartBuilders here.
@inline function CableDesign(builders...)
	parts = build_design(0.0, builders)
	return CableDesign{typeof(parts)}(parts)
end


# ---------------------------------------------------------
# The DSL Hook (Intent Capture & Auto-Grouping)
# ---------------------------------------------------------
# The user provided explicit topological groups. All good.
@inline CableDesign(layers::Tuple{Vararg{<:Gridspace{GroupBuilder}}}) =
	Gridspace{CableDesign}(layers)

# The user provided naked 1D physics parts. Default-group them.
@inline function CableDesign(layers::Tuple{Vararg{<:Gridspace{PartBuilder}}})
	grids = (
		Grid(Val{PartGroup}()),
		Grid(((0.0, 0.0),)), # Default origin at center
		Grid(1),             # Default n
		Grid(1),             # Default m
		layers...,            # Splat the naked part spaces
	)

	default_group = Gridspace{GroupBuilder}(grids)
	return Gridspace{CableDesign}((default_group,))
end

# Mixed garbage.
@inline function CableDesign(layers::Tuple)
	throw(
		ArgumentError(
			"Topological violation: You cannot mix raw parts and explicit Groups " *
			"at the top level of CableDesign. Either wrap everything in Group() " *
			"or pass raw parts exclusively.",
		),
	)
end
