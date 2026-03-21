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


export Material, CableDesign, PartGroup
export Grid, AbsoluteError

# ==========================================
# THE FRONTEND API (Compilation boundary)
# ==========================================
export Conductor, Insulator, Group

@inline function Group(layers::Tuple; origin = (0.0, 0.0), n = 1, m = 1)
	return Builder(PartGroup, origin, n, m, layers)
end

module Conductor
	import ..CableBuilder: ConductorPart, Builder
	import ..CableBuilder: Circular, Rectangular
	import ..CableBuilder: SolidCore
	import ..CableBuilder: Grid
	import ..CableBuilder: Material

	@inline function Solid(cmp::Symbol, mat; r)
		# The semicolon triggers the @gridspace kwarg interceptor.
		# If `r` is a Grid, this returns a Gridspace{Circular}.
		# If `r` is a Real, it returns a concrete Circular.
		params = Circular(; r = r)

		return Builder(ConductorPart, SolidCore, cmp, mat, params)
	end

# @inline function Tubular(cmp::Symbol, mat; t)
# 	mat_spec = convert(AbstractSpec{Material}, mat)
# 	return TubularLayerSpec(ConductorPart, Grid(cmp), Grid(t), mat_spec)
# end

# @inline function Pipe(cmp::Symbol, mat; t, filler, offset = 0.0)
# 	inner = Tubular(cmp, mat; t = t)          # tubular wall
# 	return EnclosureSpec(ConductorPart, inner, filler; offset = offset)
# end

# @inline function Stranded(
# 	cmp::Symbol,
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
# 		Grid(cmp),
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
# 	@inline function Tubular(cmp::Symbol, mat; t)
# 		mat_spec = convert(AbstractSpec{Material}, mat)
# 		return TubularLayerSpec(InsulatorPart, Grid(cmp), Grid(t), mat_spec)
# 	end
# end



end # module
