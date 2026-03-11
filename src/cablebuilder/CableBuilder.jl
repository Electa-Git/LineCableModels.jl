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


export Material, CableDesign
export Conductor, Insulator
export Grid, AbsoluteError

# ==========================================
# THE FRONTEND API (Compilation boundary)
# ==========================================
export Conductor, Insulator
module Conductor
	import ..CableBuilder: SolidCoreSpec, TubularLayerSpec, ConductorPart, Grid
	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec, CircStrandedSpec

	@inline function Solid(cmp::Symbol, mat; r)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return SolidCoreSpec(ConductorPart, Grid(cmp), Grid(r), mat_spec)
	end

	@inline function Tubular(cmp::Symbol, mat; t)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return TubularLayerSpec(ConductorPart, Grid(cmp), Grid(t), mat_spec)
	end

	@inline function Pipe(cmp::Symbol, mat; t, filler, offset = 0.0)
		inner = Tubular(cmp, mat; t = t)          # tubular wall
		return EnclosureSpec(ConductorPart, inner, filler; offset = offset)
	end

	@inline function Stranded(
		cmp::Symbol,
		mat;
		pattern::Symbol = :layer,
		r_w,
		n_w,
		lay_r,
		lay_d = 1,
	)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return CircStrandedSpec(
			ConductorPart,
			Grid(cmp),
			Grid(r_w),
			Grid(n_w),
			Grid(lay_r),
			Grid(lay_d),
			mat_spec,
		)
	end

end

module Insulator
	import ..CableBuilder: SolidCoreSpec, TubularLayerSpec, InsulatorPart, Grid
	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec

	# Insulators don't usually have solid cores, but the logic holds!
	@inline function Tubular(cmp::Symbol, mat; t)
		mat_spec = convert(AbstractSpec{Material}, mat)
		return TubularLayerSpec(InsulatorPart, Grid(cmp), Grid(t), mat_spec)
	end
end



end # module
