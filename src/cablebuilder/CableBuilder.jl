module CableBuilder

include("grid.jl")
include("gridspace.jl")
include("macros.jl")

include("materials.jl")


include("../commons/Commons.jl")
include("../uncertainbessels/UncertainBessels.jl"
)
include("../utils/Utils.jl")

include("types.jl")
include("partbuilder.jl")
include("shapes.jl")
include("cabledesign.jl")


export Material, CableDesign
export Conductor, Insulator
export Grid, AbsoluteError

# ==========================================
# THE FRONTEND API (Compilation boundary)
# ==========================================
export Conductor, Insulator
module Conductor
	import ..CableBuilder: ConductorPart, Builder
	import ..CableBuilder: Circular, Rectangular
	import ..CableBuilder: SolidCore
	import ..CableBuilder: Grid
	import ..CableBuilder: Material

	@inline function Solid(grp::Symbol, mat; r)
		# The semicolon triggers the @gridspace kwarg interceptor.
		# If `r` is a Grid, this returns a Gridspace{Circular}.
		# If `r` is a Real, it returns a concrete Circular.
		params = Circular(; r = r)

		return Builder(ConductorPart, SolidCore, grp, mat, params)
	end

# @inline function Tubular(grp::Symbol, mat; t)
# 	mat_spec = convert(AbstractSpec{Material}, mat)
# 	return TubularLayerSpec(ConductorPart, Grid(grp), Grid(t), mat_spec)
# end

# @inline function Pipe(grp::Symbol, mat; t, filler, offset = 0.0)
# 	inner = Tubular(grp, mat; t = t)          # tubular wall
# 	return EnclosureSpec(ConductorPart, inner, filler; offset = offset)
# end

# @inline function Stranded(
# 	grp::Symbol,
# 	mat;
# 	pattern::Symbol = :layer,
# 	r_w,
# 	n_w,
# 	lay_r,
# 	lay_d = 1,
# )
# 	mat_spec = convert(AbstractSpec{Material}, mat)
# 	return CircStrandedSpec(
# 		ConductorPart,
# 		Grid(grp),
# 		Grid(r_w),
# 		Grid(n_w),
# 		Grid(lay_r),
# 		Grid(lay_d),
# 		mat_spec,
# 	)
# end

end

# module Insulator
# 	import ..CableBuilder: SolidCoreSpec, TubularLayerSpec, InsulatorPart, Grid
# 	import ..CableBuilder: AbstractSpec, Material, EnclosureSpec

# 	# Insulators don't usually have solid cores, but the logic holds!
# 	@inline function Tubular(grp::Symbol, mat; t)
# 		mat_spec = convert(AbstractSpec{Material}, mat)
# 		return TubularLayerSpec(InsulatorPart, Grid(grp), Grid(t), mat_spec)
# 	end
# end



end # module
