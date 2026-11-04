"""
    LineCableModels.CAD

Draw the outlines of the cable regions on a target.
Walk a system or design through its regions, contours and paths to their curves.
Targets define how to draw each curve in their own primitives.
"""
module CAD

using DocStringExtensions: TYPEDSIGNATURES
using ..DataModel: AbstractShape, PlacedRegion, CableDesign, CableSystem,
                   Curve, Contour

public draw

"""
$(TYPEDSIGNATURES)

Draw on `target` in the target's own primitives. The walk reads each system's or design's
already placed regions, then `placed.shape`, its contour, paths and curves. Coordinates
are in \\[m\\] and angles in \\[rad\\].

A target defines methods for [`Segment`](@ref LineCableModels.DataModel.Segment) and
[`Arc`](@ref LineCableModels.DataModel.Arc), and can override any level above them. Each
method returns what the target made.

A system or design returns a vector of results from its regions.
A placed region returns the result for `placed.shape`.
An `AbstractShape` returns the result for its contour.
A contour returns a named tuple of `outer` and `holes` results.
For a path, the result is a vector of results from its curves.
"""
function draw end

draw(target, system::CableSystem) = [draw(target, placed) for placed in system.regions]
draw(target, design::CableDesign) = [draw(target, placed) for placed in design.regions]
draw(target, placed::PlacedRegion) = draw(target, placed.shape)
draw(target, shape::AbstractShape) = draw(target, Contour(shape))
function draw(target, contour::Contour)
    return (outer = draw(target, contour.outer),
        holes = [draw(target, hole) for hole in contour.holes])
end
draw(target, path::AbstractVector{<:Curve}) = [draw(target, curve) for curve in path]

end
