module CableBuilder

include("spec.jl")
include("grid.jl")

include("materials.jl")


include("../commons/Commons.jl")
include("../uncertainbessels/UncertainBessels.jl"
)
include("../utils/Utils.jl")
include("../validation/Validation.jl")

include("shapes.jl")
include("types.jl")
include("cabledesign.jl")


export Material,
	MaterialSpec, build, Conductor, AbstractCablePart, CableDesign, CableDesignSpec
export SolidCore, TubularShape, Concentric, ConductorPart, Enclosure
export Grid, DeterministicGrid, RelativeGrid, AbsoluteGrid, AbsoluteError, AbstractSpec,
	EnclosureSpec

# ==========================================
# THE FRONTEND API (Compilation boundary)
# ==========================================
export SolidCore, TubularShape, Concentric, ConductorPart

module Conductor
	import ..CableBuilder: SolidCoreSpec, TubularPartSpec, ConductorPart, Grid
	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec

	@inline function Solid(tag::Symbol, mat; r)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return SolidCoreSpec(ConductorPart, Grid(tag), Grid(r), mat_spec)
	end

	@inline function Tubular(tag::Symbol, mat; t)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return TubularPartSpec(ConductorPart, Grid(tag), Grid(t), mat_spec)
	end

	@inline function Pipe(tag::Symbol, mat; t, filler, offset = 0.0)
		inner = Tubular(tag, mat; t = t)          # tubular wall
		return EnclosureSpec(ConductorPart, inner, filler; offset = offset)
	end
end

module Insulator
	import ..CableBuilder: SolidCoreSpec, TubularPartSpec, InsulatorPart, Grid
	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec

	# Insulators don't usually have solid cores, but the logic holds!
	@inline function Tubular(tag::Symbol, mat; t)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return TubularPartSpec(InsulatorPart, Grid(tag), Grid(t), mat_spec)
	end
end



end # module
