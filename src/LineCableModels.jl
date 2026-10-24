module LineCableModels

include("commons/Commons.jl")
include("materials/Materials.jl")
include("earth/Earth.jl")
include("datamodel/DataModel.jl")

using .Commons: resolve, validate
using .Materials: Material, material
using .Earth: EarthLayer, EarthModel, homogeneous
using .DataModel: Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle,
                  Shell, BentStrip, CircularPattern, Helix, LayRatio, area, centroid,
                  support,
                  boundary, r_in, r_ex, placements, overlength, pitch

export Material, material, EarthLayer, EarthModel, homogeneous
export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Shell
export CircularPattern, Helix, LayRatio
export area, centroid, support, boundary, r_in, r_ex, placements, overlength, pitch
public Commons, Materials, Earth, DataModel, validate, resolve, BentStrip

end
