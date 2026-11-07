"""
    LineCableModels.BaseGeometry

Describe the plane geometry of a cross-section, free of any cable part. It defines the
placement, the segments, circular arcs and elliptic arcs, and the contour of closed
paths, with the methods on them from the areas and supports to the parallels, the
clearances and the distance across. The verb [`draw`](@ref) walks a contour through its
paths to their curves, and a target defines how to draw each curve in its own primitives.
"""
module BaseGeometry

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: resolve, validate
using ..Commons: nominal
using QuadGK: quadgk
using IntervalSets: ClosedInterval, leftendpoint, rightendpoint

export Placement, area, perimeter, centroid, support, distance
public Segment, Arc, EllipticArc, Curve, Contour, turn, parallel, clearance, dimension,
       direction, lowest, draw, overlap, vertex, apsides, winding

include("placement.jl")
include("contour.jl")

"""
$(TYPEDSIGNATURES)

Draw on `target` in the target's own primitives. The walk reads a contour's paths, then
their curves. Coordinates are in \\[m\\] and angles in \\[rad\\].

A target defines methods for [`Segment`](@ref) and [`Arc`](@ref), and can override any
level above them. Each method returns what the target made.

A contour returns a named tuple of `outer` and `holes` results.
For a path, the result is a vector of results from its curves.
"""
function draw end

function draw(target, contour::Contour)
    return (outer = draw(target, contour.outer),
        holes = [draw(target, hole) for hole in contour.holes])
end
draw(target, path::AbstractVector{<:Curve}) = [draw(target, curve) for curve in path]

end
