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
using QuadGK: quadgk

export Placement, AbstractShape, AbstractPrimitive, Circle, Annulus, Rectangle, Offset
export Sector, Ellipse, Polygon
export AbstractLayoutPattern, CircularPattern, LinearPattern, Helix, LayRatio
export area, perimeter, centroid, support, distance, boundary, r_in, r_ex, placements,
       overlength,
       pitch
export PerforatedShape, ArcPolygon
public Segment, Arc, EllipticArc, Curve, Contour, turn, parallel, clearance, dimension,
       direction
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
