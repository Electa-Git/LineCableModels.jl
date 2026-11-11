"""
$(TYPEDEF)

Represent a region bounded by one closed path of segments, circular arcs and elliptic arcs.

$(TYPEDFIELDS)
"""
struct ArcPolygon{T <: Real} <: AbstractShape{T}
    "Closed path of curves that bounds the region, counterclockwise, in the frame of the cable."
    outer::Vector{Curve{T}}
    @doc """
    $(TYPEDSIGNATURES)

    Construct a region from its counterclockwise `outer` path, with coordinates in \\[m\\]
    and angles in \\[rad\\], through `validate`.
    """
    function ArcPolygon(outer::Vector{Curve{T}}) where {T <: Real}
        return validate(new{T}(outer))
    end
end

"""
$(TYPEDSIGNATURES)

Check that the outline of `shape` is one closed path, by the checks of `Contour`, and
return `shape`.
Coordinates are in \\[m\\] and angles in \\[rad\\].
"""
function validate(shape::ArcPolygon{T}) where {T}
    Contour(shape.outer, Vector{Curve{T}}[])
    return shape
end

"""
$(TYPEDSIGNATURES)

Return `shape.outer` as the outer path with no holes.
Coordinates are in \\[m\\] and angles in \\[rad\\].
"""
function Contour(shape::ArcPolygon{T}) where {T}
    return Contour(shape.outer, Vector{Curve{T}}[])
end

area(shape::ArcPolygon) = area(Contour(shape))
centroid(shape::ArcPolygon) = centroid(Contour(shape))
perimeter(shape::ArcPolygon) = perimeter(Contour(shape))
boundary(shape::ArcPolygon) = shape

function Base.issubset(B::AbstractShape, A::Union{Rectangle, Polygon, Ellipse, ArcPolygon})
    container = Contour(A)
    content = Contour(B)
    normals = sort!(
        reduce(vcat,
            vec([map(curve -> mod(direction(curve, t) - π / 2, 2π), path)
                 for t in (0, 1), path in (container.outer, content.outer)]));
        by = nominal)
    push!(normals, first(normals) + 2π)
    for search in (false, true), k in 1:(length(normals) - 1)

        low, high = normals[k], normals[k + 1]
        φ = low
        if search
            nominal(high) > nominal(low) || continue
            while true
                left = low + (high - low) / 3
                right = high - (high - low) / 3
                (nominal(left) == nominal(low) || nominal(right) == nominal(high) ||
                 nominal(left) == nominal(right)) && break
                if nominal(support(container, left) - support(content, left)) <
                   nominal(support(container, right) - support(content, right))
                    high = right
                else
                    low = left
                end
            end
            φ = argmin(
                angle -> nominal(support(container, angle) - support(content, angle)),
                (low, (low + high) / 2, high))
        end
        radius = support(container, φ)
        nominal(support(content, φ)) <= nominal(radius) + tolerance(radius) || return false
    end
    return true
end

"""
$(TYPEDSIGNATURES)

Return `true` when the interiors of `B` and of the circle `hole` are disjoint. The
distance of `B` from the center of the circle is then at least its radius ``r`` less the
tolerance ``√eps(float(r))max(1,r)``. Lengths are in \\[m\\]. A content that touches the
hole lies apart from it, as a content that touches its container lies within it.
"""
function Base.isdisjoint(B::AbstractShape, hole::Circle)
    radius = hole.r
    return nominal(distance(resolve(Placement(-hole.at.x, -hole.at.y), B))) +
           tolerance(radius) >= nominal(radius)
end

"""
$(TYPEDSIGNATURES)

Return `true` when the interiors of `B` and `hole` are disjoint, by a direction ``ψ``
that separates them: ``support(hole,ψ)≤-support(B,ψ+π)+tol``. The directions are the
outward normals of the segments of the hole and the reversed normals of the segments of
`B`. The test is exact for polygonal outlines. Only segments give directions: beside an
arc of a hole other than a circle, the method returns `false`.
Lengths are in \\[m\\], angles in \\[rad\\], and the tolerance is that of a circular hole.
"""
function Base.isdisjoint(B::AbstractShape, hole::AbstractShape)
    content = Contour(B)
    directions = [atan(-(edge.q[1] - edge.p[1]), edge.q[2] - edge.p[2])
                  for edge in Contour(hole).outer if edge isa Segment]
    append!(directions,
        [atan(-(edge.q[1] - edge.p[1]), edge.q[2] - edge.p[2]) + π
         for edge in content.outer if edge isa Segment])
    apart = false
    for ψ in directions
        radius = support(content, ψ + π)
        apart |= nominal(support(hole, ψ)) <= -nominal(radius) + tolerance(radius)
    end
    return apart
end

"""
$(TYPEDSIGNATURES)

Return whether ``B⊆A.outer`` and `B` lies apart from every hole of `A`, by `isdisjoint`.
Lengths are in \\[m\\].
"""
function Base.issubset(B::AbstractShape, A::PerforatedShape)
    return B ⊆ A.outer && all(hole -> isdisjoint(B, hole), A.holes)
end

"""
$(TYPEDSIGNATURES)

Place the outline by ``P\\mapsto\\mathrm{at}(P)`` for segment ends and
``c\\mapsto\\mathrm{at}\\,c`` for arc placements. Coordinates are in \\[m\\]
and rotations in \\[rad\\]. Keep each elliptic arc's normal offset \\[m\\].
"""
function resolve(at::Placement, shape::ArcPolygon)
    T = promote_type(eltype(at), eltype(shape))
    return ArcPolygon(Curve{T}[resolve(at, curve) for curve in shape.outer])
end

"""
$(TYPEDSIGNATURES)

Map the curves' coordinates \\[m\\], radii \\[m\\], angles \\[rad\\] and normal offsets \\[m\\] of `shape` by ``x\\mapsto f(x)``.
"""
function Base.map(f, shape::ArcPolygon)
    outer = similar(shape.outer)
    for (k, curve) in enumerate(shape.outer)
        outer[k] = map(f, curve)
    end
    return ArcPolygon(outer)
end

"""
$(TYPEDSIGNATURES)

Construct four counterclockwise segments from
``\\mathrm{at}((-w/2,-h/2))`` around `rectangle`.
Width, height and coordinates are in \\[m\\].
"""
function Contour(rectangle::Rectangle{T}) where {T}
    w, h, at = rectangle.w, rectangle.h, rectangle.at
    points = [at((-w / 2, -h / 2)), at((w / 2, -h / 2)),
        at((w / 2, h / 2)), at((-w / 2, h / 2))]
    return Contour(Curve{T}[Segment(points[k], points[mod1(k + 1, 4)]) for k in 1:4],
        Vector{Curve{T}}[])
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

With displaced sides or fillets, let ``α=span/2``, ``g=apex\\sin α`` and ``ρ=fillet``.
The upper fillet centers ``c`` satisfy ``(\\sin α,-\\cos α)\\cdot c=g+ρ``, with
``|c|=ro-ρ`` at the back and ``ri+ρ`` at the base. The contacts are their projections on
the side and on its circle. The tip center is ``(apex+ρ/\\sin α,0)``. Reflection gives
the lower contacts, and the placement moves every curve to the frame of the cable.
"""
function Contour(sector::Sector{T}) where {T}
    c = sector.at
    if !(iszero(sector.apex) && iszero(nominal(sector.fillet)))
        α = sector.span / 2
        g = sector.apex * sin(α)
        ρ = sector.fillet
        β = α - asin((g + ρ) / (sector.ro - ρ))
        outer = Arc(Placement(c.x, c.y, c.φ - β), sector.ro, 2β)
        if iszero(nominal(ρ))
            if iszero(nominal(sector.ri))
                tip = c((sector.apex, zero(sector.apex)))
                return Contour(
                    Curve{T}[outer, Segment(outer(1), tip),
                        Segment(tip, outer(0))], Vector{Curve{T}}[])
            end
            γ = α - asin(g / sector.ri)
            inner = Arc(Placement(c.x, c.y, c.φ + γ), sector.ri, -2γ)
            return Contour(
                Curve{T}[outer, Segment(outer(1), inner(0)), inner,
                    Segment(inner(1), outer(0))],
                Vector{Curve{T}}[])
        end
        # The upper side has inward normal (sin α, -cos α).
        # Its inner parallel meets the back circle at the back fillet's center.
        d = sqrt((sector.ro - ρ)^2 - (g + ρ)^2)
        x = d * cos(α) + (g + ρ) * sin(α)
        y = d * sin(α) - (g + ρ) * cos(α)
        upper = c((x, y))
        lower = c((x, -y))
        upper_back = Arc(Placement(upper..., c.φ + β), ρ, α + π / 2 - β)
        lower_back = Arc(Placement(lower..., c.φ - α - π / 2), ρ, α + π / 2 - β)
        if nominal(sector.ri) > zero(nominal(sector.ri))
            d = sqrt((sector.ri + ρ)^2 - (g + ρ)^2)
            x = d * cos(α) + (g + ρ) * sin(α)
            y = d * sin(α) - (g + ρ) * cos(α)
            upper = c((x, y))
            lower = c((x, -y))
            γ = α - asin((g + ρ) / (sector.ri + ρ))
            upper_base = Arc(Placement(upper..., c.φ + α + π / 2),
                ρ, γ + π / 2 - α)
            inner = Arc(Placement(c.x, c.y, c.φ + γ), sector.ri, -2γ)
            lower_base = Arc(Placement(lower..., c.φ + π - γ),
                ρ, γ + π / 2 - α)
            return Contour(
                Curve{T}[outer, upper_back,
                    Segment(upper_back(1), upper_base(0)), upper_base, inner, lower_base,
                    Segment(lower_base(1), lower_back(0)), lower_back],
                Vector{Curve{T}}[])
        end
        tip = c((sector.apex + ρ / sin(α), zero(ρ)))
        tip_arc = Arc(Placement(tip..., c.φ + α + π / 2), ρ, π - sector.span)
        return Contour(
            Curve{T}[outer, upper_back,
                Segment(upper_back(1), tip_arc(0)), tip_arc,
                Segment(tip_arc(1), lower_back(0)), lower_back],
            Vector{Curve{T}}[])
    end
    if sector.span == 2π
        if iszero(nominal(sector.ri))
            return Contour(Circle(sector.ro, c))
        end
        return Contour(Annulus(sector.ri, sector.ro, c))
    end
    outer = Arc(Placement(c.x, c.y, c.φ - sector.span / 2), sector.ro, sector.span)
    if iszero(nominal(sector.ri))
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
function Contour(shape::PerforatedShape)
    contour = Contour(shape.outer)
    holes = copy(contour.holes)
    append!(holes, [Contour(hole).outer for hole in shape.holes])
    return Contour(contour.outer, holes)
end

"""
$(TYPEDSIGNATURES)

Construct the four counterclockwise elliptic quadrants of `ellipse`,
``φ0=kπ/2``, ``span=π/2`` \\[rad\\], for ``k=0,…,3``.
Semi-axes and center coordinates are in \\[m\\].
"""
function Contour(ellipse::Ellipse{T}) where {T}
    return Contour(
        Curve{T}[EllipticArc(ellipse.at, ellipse.a, ellipse.b, T(k * π / 2),
                     T(π / 2)) for k in 0:3],
        Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Construct the placed edges ``Segment(at(v_i),at(v_{i+1}))`` \\[m\\] of `polygon`,
Close the last edge to the first vertex. The contour has no holes.
"""
function Contour(polygon::Polygon{T}) where {T}
    points = [polygon.at(vertex) for vertex in polygon.vertices]
    return Contour(
        Curve{T}[Segment(points[i], points[mod1(i + 1, length(points))])
                 for i in eachindex(points)],
        Vector{Curve{T}}[])
end

"""
$(TYPEDSIGNATURES)

Return the distance across the metal of `shape` \\[m\\]: the distance across the region
that its contour bounds, twice the depth at which the inward parallel of its outline first
degenerates.
"""
dimension(shape::AbstractShape) = dimension(Contour(shape))
