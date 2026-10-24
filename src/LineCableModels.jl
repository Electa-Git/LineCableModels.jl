module LineCableModels

include("commons/Commons.jl")
include("materials/Materials.jl")
include("earth/Earth.jl")

using .Commons: validate
using .Materials: Material, material
using .Earth: EarthLayer, EarthModel, homogeneous

export Material, material, EarthLayer, EarthModel, homogeneous
public Commons, Materials, Earth, validate

end
