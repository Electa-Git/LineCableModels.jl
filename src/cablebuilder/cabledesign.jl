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

	# Extract boundary logic respecting the new absolute origin
	ox, oy = target.origin
	part_boundary = sqrt(ox^2 + oy^2) + r_ex(target)

	# The running radius must engulf the largest eccentric boundary
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

# The Tuple hook (no splatting required at call site)
# Usage: parts = (Conductor.Solid(...), ...); CableDesign(parts)
@inline CableDesign(layers::Tuple) = Gridspace{CableDesign}(layers)