module LineCableModels

include("commons/Commons.jl")
include("materials/Materials.jl")
include("earth/Earth.jl")
include("datamodel/DataModel.jl")

using .Commons: build, resolve, validate
using .Materials: Material, material
using .Earth: EarthLayer, EarthModel, homogeneous
using .DataModel: Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle,
                  Shell, BentStrip, CircularPattern, Helix, LayRatio, area, centroid,
                  support, distance, boundary, r_in, r_ex, placements, overlength, pitch,
                  DifferenceShape, AbstractCablePart, Region, Stack, Terminal, Repeat, Lay,
                  Enclosure, PlacedRegion, CableDesign, CableSystem, CableBlueprint,
                  BlueprintConductor, BlueprintDielectric, flatten

export Material, material, EarthLayer, EarthModel, homogeneous
export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Shell
export CircularPattern, Helix, LayRatio
export area, centroid, support, distance, boundary, r_in, r_ex, placements, overlength,
       pitch
export AbstractCablePart, Region, Stack, Terminal, Repeat, Lay, Enclosure, CableDesign,
       CableSystem, build
public Commons, Materials, Earth, DataModel, validate, resolve, BentStrip, DifferenceShape
public PlacedRegion, CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten

end
