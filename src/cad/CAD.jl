"""
    LineCableModels.CAD

Describe the planar curves that bound cable regions and draw their contours on a target.
Walk a system or design through its regions, contours and paths to their curves.
Targets define how to draw each curve in their own primitives.
"""
module CAD

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using ..DataModel: Placement, AbstractShape, Circle, Annulus, AnnularSector,
                   DifferenceShape, PlacedRegion, CableDesign, CableSystem
import ..Commons: validate

public Segment, Arc, Curve, Contour, draw

"""
$(TYPEDEF)

Represent a straight curve between two points in the cross-section plane.

$(TYPEDFIELDS)
"""
struct Segment{T <: Real}
    "Coordinates at `t = 0` in \\[m\\]."
    a::NTuple{2, T}
    "Coordinates at `t = 1` in \\[m\\]."
    b::NTuple{2, T}
    function Segment{T}(a::NTuple{2, T}, b::NTuple{2, T}) where {T <: Real}
        return validate(new{T}(a, b))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a segment from the points `a` and `b` \\[m\\], in the floating type that their
coordinates promote to.

# Errors

- Throws `DomainError` when a coordinate is not finite or the points coincide.
"""
function Segment(a::Tuple{Real, Real}, b::Tuple{Real, Real})
    coordinates = map(float, promote(a..., b...))
    T = eltype(coordinates)
    return Segment{T}((coordinates[1], coordinates[2]), (coordinates[3], coordinates[4]))
end

function validate(segment::Segment)
    all(isfinite, (segment.a..., segment.b...)) || throw(DomainError(
        (segment.a, segment.b), "Segment coordinates must be finite"))
    segment.a != segment.b || throw(DomainError(
        (segment.a, segment.b), "Segment endpoints must differ"))
    return segment
end

"""
$(TYPEDEF)

Represent a circular curve from the x-axis of its placement over a signed angular span.
A positive span runs counterclockwise and a negative span clockwise.

$(TYPEDFIELDS)
"""
struct Arc{T <: Real}
    "Center \\[m\\] and initial direction \\[rad\\]."
    at::Placement{T}
    "Radius \\[m\\]."
    r::T
    "Signed angular span \\[rad\\]."
    span::T
    function Arc{T}(at::Placement{T}, r::T, span::T) where {T <: Real}
        return validate(new{T}(at, r, span))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an arc at `at` with radius `r` \\[m\\] and signed angular `span` \\[rad\\].
The placement translates in \\[m\\] and rotates in \\[rad\\].
The arguments promote to a shared floating type.

# Errors

- Throws `DomainError` when `r` is not positive and finite, or when `span` is not finite
  or does not satisfy ``0 < |\\mathrm{span}| \\leq 2\\pi``.
"""
function Arc(at::Placement, r::Real, span::Real)
    T = float(promote_type(eltype(at), typeof(r), typeof(span)))
    return Arc{T}(convert(Placement{T}, at), convert(T, r), convert(T, span))
end

function validate(arc::Arc)
    isfinite(arc.r) && arc.r > zero(arc.r) ||
        throw(DomainError(arc.r, "Arc.r must be positive and finite"))
    isfinite(arc.span) && zero(arc.span) < abs(arc.span) <= 2π ||
        throw(DomainError(arc.span, "Arc.span must be finite with 0 < |span| ≤ 2π"))
    return arc
end

"""
A curve of a contour is a [`Segment`](@ref) or an [`Arc`](@ref) in the floating type `T`.
"""
const Curve{T} = Union{Segment{T}, Arc{T}}

"""
$(TYPEDEF)

Represent the closed paths that bound a region of the plane. Each path runs
counterclockwise about the region it bounds. Each hole follows the same orientation.

$(TYPEDFIELDS)
"""
struct Contour{T <: Real}
    "Outer path of curves with coordinates in \\[m\\] and angles in \\[rad\\]."
    outer::Vector{Curve{T}}
    "Paths of holes with coordinates in \\[m\\] and angles in \\[rad\\]."
    holes::Vector{Vector{Curve{T}}}
    @doc """
    $(TYPEDSIGNATURES)

    Construct a contour from its `outer` path and the paths of its `holes`, whose
    coordinates are in \\[m\\] and angles in \\[rad\\]. Each path is nonempty and closed.
    The tolerance is `sqrt(eps(T))` times the largest absolute coordinate at a path's ends.
    That scale is at least one meter.

    # Errors

    - Throws `DomainError` if a path is empty or its ends do not meet.
    """
    function Contour(outer::Vector{Curve{T}},
            holes::Vector{Vector{Curve{T}}}) where {T <: Real}
        return validate(new{T}(outer, holes))
    end
end

function validate(contour::Contour{T}) where {T}
    for path in (contour.outer, contour.holes...)
        isempty(path) && throw(DomainError(path, "Contour paths must be nonempty"))
        scale = max(
            one(T), maximum(abs(coordinate)
            for curve in path for t in (0, 1) for coordinate in curve(t)))
        tolerance = sqrt(eps(T)) * scale
        for k in eachindex(path)
            endpoint = path[k](1)
            nextpoint = path[mod1(k + 1, length(path))](0)
            all(abs.(endpoint .- nextpoint) .<= tolerance) ||
                throw(DomainError(path, "Contour path must be closed at curve $k"))
        end
    end
    return contour
end

"""
$(TYPEDSIGNATURES)

Evaluate the segment at `t` \\[dimensionless\\] in `[0, 1]`:

```math
\\boldsymbol{p}(t) = \\boldsymbol{a} + t(\\boldsymbol{b} - \\boldsymbol{a}),
```

where ``\\boldsymbol{a}`` and ``\\boldsymbol{b}`` are its endpoint coordinates \\[m\\].
Return coordinates in \\[m\\].
"""
function (segment::Segment)(t::Real)
    return segment.a .+ t .* (segment.b .- segment.a)
end

"""
$(TYPEDSIGNATURES)

Evaluate the arc at `t` \\[dimensionless\\] in `[0, 1]`:

```math
\\boldsymbol{p}(t) = \\mathrm{at}\\bigl(r\\cos(t\\,\\mathrm{span}),
    r\\sin(t\\,\\mathrm{span})\\bigr),
```

where ``r`` is its radius \\[m\\], ``\\mathrm{span}`` its signed angular span \\[rad\\], and
``\\mathrm{at}`` rotates by the initial direction \\[rad\\] and translates by the center
coordinates \\[m\\]. Return coordinates in \\[m\\].
"""
function (arc::Arc)(t::Real)
    return arc.at((arc.r * cos(t * arc.span), arc.r * sin(t * arc.span)))
end

"""
$(TYPEDSIGNATURES)

Construct the four counterclockwise quadrant arcs of `circle`, whose radius and center
coordinates are in \\[m\\]. Arc ``k = 0, \\ldots, 3`` starts at the angle ``k \\pi / 2`` and
spans ``\\pi / 2`` \\[rad\\], independent of the placement's rotation. The contour has no
holes.
"""
function Contour(circle::Circle{T}) where {T}
    c = circle.at
    outer = Curve{T}[Arc(Placement(c.x, c.y, T(k * π / 2)), circle.r, T(π / 2))
                     for k in 0:3]
    return Contour(outer, Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Construct the counterclockwise quadrants of `annulus`: its outer radius bounds the outer
path and its inner radius bounds one hole. Radii and center coordinates are in \\[m\\].
Both paths start at angle zero, independent of the placement's rotation \\[rad\\].
"""
function Contour(annulus::Annulus)
    outer = Contour(Circle(annulus.ro, annulus.at)).outer
    hole = Contour(Circle(annulus.ri, annulus.at)).outer
    return Contour(outer, [hole])
end

"""
$(TYPEDSIGNATURES)

Construct the closed path of `sector`, whose radii and center coordinates are in \\[m\\]
and whose angular span and placement rotation are in \\[rad\\]. A partial turn of span
``\\Delta\\phi`` about the rotation ``\\phi`` follows the outer arc from
``\\phi - \\Delta\\phi / 2`` with span ``\\Delta\\phi``, the radial edge, the inner arc
from ``\\phi + \\Delta\\phi / 2`` with span ``-\\Delta\\phi`` and the closing radial
edge. With a
zero inner radius, two segments join the outer arc to the center. A full turn gives the
contour of a circle or annulus.
"""
function Contour(sector::AnnularSector{T}) where {T}
    c = sector.at
    if sector.span == 2π
        if iszero(sector.ri)
            return Contour(Circle(sector.ro, c))
        end
        return Contour(Annulus(sector.ri, sector.ro, c))
    end
    outer = Arc(Placement(c.x, c.y, c.φ - sector.span / 2), sector.ro, sector.span)
    if iszero(sector.ri)
        center = (c.x, c.y)
        inward = Segment(outer(1), center)
        outward = Segment(center, outer(0))
        return Contour(Curve{T}[outer, inward, outward], Vector{Curve{T}}[])
    end
    inner = Arc(Placement(c.x, c.y, c.φ + sector.span / 2), sector.ri, -sector.span)
    inward = Segment(outer(1), inner(0))
    outward = Segment(inner(1), outer(0))
    return Contour(Curve{T}[outer, inward, inner, outward], Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Construct the contour of `shape`, whose dimensions and positions are in \\[m\\] and angles
in \\[rad\\]. Keep the contour of `shape.outer` and append `Contour(hole).outer` for each
member of `shape.holes` to its holes. All paths run counterclockwise about the regions they
bound.
"""
function Contour(shape::DifferenceShape)
    contour = Contour(shape.outer)
    holes = copy(contour.holes)
    append!(holes, [Contour(hole).outer for hole in shape.holes])
    return Contour(contour.outer, holes)
end

"""
$(TYPEDSIGNATURES)

Draw on `target` in the target's own primitives. The walk reads each system's or design's
already placed regions, then `placed.shape`, its contour, paths and curves. Coordinates
are in \\[m\\] and angles in \\[rad\\].

A target defines methods for [`Segment`](@ref) and [`Arc`](@ref), and can override any
level above them. Each method returns what the target made.

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
