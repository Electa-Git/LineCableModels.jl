"""
    LineCableModels.DataModel

Describe the cross-section of a cable: the geometry of its parts, and the rings and helical
lays of its repeated wires.
"""
module DataModel

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: resolve, validate

export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Shell
export BentStrip
export CircularPattern, Helix, LayRatio
export area, centroid, support, boundary, r_in, r_ex, placements, overlength, pitch

include("geometry.jl")
include("placement.jl")

end
