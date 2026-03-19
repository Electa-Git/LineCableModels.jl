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
	return (target, build_design(r_ex(target.shape), Base.tail(builders))...)
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
