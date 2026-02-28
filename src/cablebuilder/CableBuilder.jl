module CableBuilder

include("spec.jl")
include("grid.jl")

include("materials.jl")


include("../../src/commons/Commons.jl")
include("../../src/uncertainbessels/UncertainBessels.jl"
)
include("../../src/utils/Utils.jl")
include("../../src/validation/Validation.jl")

include("shapes.jl")
include("types.jl")
include("cabledesign.jl")


export Material,
	MaterialSpec, build, Conductor, AbstractCablePart, CableDesign, CableDesignSpec
export SolidCore, TubularShape, Concentric, ConductorPart
export Grid, DeterministicGrid, RelativeGrid, AbsoluteGrid, AbsoluteError, AbstractSpec

# ==========================================
# THE FRONTEND API (Compilation boundary)
# ==========================================
export SolidCore, TubularShape, Concentric, ConductorPart

module Conductor
	import ..CableBuilder: SolidCoreSpec, TubularPartSpec, ConductorPart, Grid
	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec

	function Solid(tag::Symbol, mat; r)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return SolidCoreSpec(Grid(ConductorPart), Grid(tag), Grid(r), mat_spec)
	end

	function Tubular(tag::Symbol, mat; t)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return TubularPartSpec(Grid(ConductorPart), Grid(tag), Grid(t), mat_spec)
	end
end

module Insulator
	import ..CableBuilder: SolidCoreSpec, TubularPartSpec, InsulatorPart, Grid
	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec

	# Insulators don't usually have solid cores, but the logic holds!
	function Tubular(tag::Symbol, mat; t)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return TubularPartSpec(Grid(InsulatorPart), Grid(tag), Grid(t), mat_spec)
	end
end


# ---------------------------------------------------------
# The Cheap-Allocation Stacking Engine
# ---------------------------------------------------------
# 1. Base Case: We ran out of builders. Return an empty tuple.
@inline build_layer(current_r) = ()

# 2. Recursive Step: Build the current layer, update radius, and recursively process the rest.
@inline function build_layer(current_r, builder, rest...)
	part = builder(current_r)
	# Construct the strictly-typed tuple entirely on the stack
	return (part, build_layer(part.shape.r_ex, rest...)...)
end

# 3. The Master Constructor
@inline function CableDesign(builders...)
	# Kick off the compile-time recursion starting at r = 0.0
	# No arrays. No for-loops. No dynamic dispatch.
	parts_tuple = build_layer(0.0, builders...)

	return CableDesign(parts_tuple)
end

end # module
