"""
    LineCableModels.DataModel

Describe the cross-section of a cable. Its declared parts build a design, which holds the
placed regions and the tables that map them to terminals, repeats, lays and fills. Designs
placed in one cross-section of a line, with the conductor index of each terminal, build a
system. The macros [`@design`](@ref) and [`@terminal`](@ref) and the builders of the layers,
the stranded conductor, the courses of wires and the tapes declare the parts of a cable. The
formations return the placements of the cables of a system.
"""
module DataModel

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: build, resolve, validate
using ..Materials: Material
using ..Commons: nominal
using SpecialFunctions: ellipe
import ..BaseGeometry: Contour, area, perimeter, centroid, support, distance, dimension,
                       draw
using ..BaseGeometry: Placement, Segment, Arc, EllipticArc, Curve, parallel, direction

export AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Offset
export Sector, Ellipse, Polygon
export AbstractLayoutPattern, CircularPattern, LinearPattern, Helix, LayRatio
export boundary, r_in, r_ex, placements, istangent,
       overlength,
       pitch
export PerforatedShape, ArcPolygon
export AbstractCablePart, Region, Stack, Terminal, Repeat, Lay, Enclosure, Assembly
export PlacedRegion, CableDesign, CableSystem
export CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten
export @design, @terminal, core, insulation, screen, sheath, bedding, jacket, stranded,
       wires, tape, trefoil, flat

include("geometry.jl")
include("contour.jl")
include("placement.jl")
include("parts.jl")
include("design.jl")
include("system.jl")
include("blueprint.jl")
include("builders.jl")

end
