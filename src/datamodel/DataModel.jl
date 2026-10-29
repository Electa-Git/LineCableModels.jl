"""
    LineCableModels.DataModel

Describe the cross-section of a cable. Its declared parts build a design, which holds the
placed regions and the tables that map them to terminals, repeats, lays and fills. Designs
placed in one cross-section of a line, with the conductor index of each terminal, build a
system. The macros [`@cable`](@ref) and [`@terminal`](@ref) and the builders of the layers,
the stranded conductor and the armor declare the parts of a cable.
"""
module DataModel

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: build, resolve, validate
using ..Materials: Material

export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Offset
export AnnularSector
export AbstractLayoutPattern, CircularPattern, Helix, LayRatio
export area, centroid, support, distance, boundary, r_in, r_ex, placements, overlength,
       pitch
export DifferenceShape
export AbstractCablePart, Region, Stack, Terminal, Repeat, Lay, Enclosure
export PlacedRegion, CableDesign, CableSystem
export CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten
export @cable, @terminal, core, insulation, screen, sheath, bedding, jacket, stranded, armor

include("geometry.jl")
include("placement.jl")
include("parts.jl")
include("design.jl")
include("system.jl")
include("blueprint.jl")
include("builders.jl")

end
