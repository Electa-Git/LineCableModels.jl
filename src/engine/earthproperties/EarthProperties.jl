"""
    LineCableModels.Engine.EarthProperties

Frequency-dependent evaluation of static `EarthModel` layers.
"""
module EarthProperties

export CPEarth, evaluate

using DocStringExtensions: TYPEDEF, TYPEDSIGNATURES
import ..Engine: AbstractEarthPropertiesFormulation, description
using ...EarthProps: EarthModel

include("constant.jl")

end
