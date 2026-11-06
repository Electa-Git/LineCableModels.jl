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
support(shape::ArcPolygon) = support(Contour(shape))
support(shape::ArcPolygon, φ::Real) = support(Contour(shape), φ)
distance(shape::ArcPolygon) = distance(Contour(shape))
boundary(shape::ArcPolygon) = shape

function Base.issubset(B::AbstractShape, A::Union{Rectangle, Polygon, Ellipse, ArcPolygon})
    T = promote_type(eltype(A), eltype(B))
    paths = (Contour(A).outer, Contour(B).outer)
    for edge in [curve for path in paths for curve in path if curve isa Segment]
        dx, dy = edge.q .- edge.p
        φ = atan(-dx, dy)
        radius = support(A, φ)
        tolerance = sqrt(eps(float(nominal(radius)))) *
                    max(one(nominal(radius)), nominal(radius))
        nominal(support(B, φ)) <= nominal(radius) + tolerance || return false
    end
    arcs = [[(at = curve.at, a = curve.r, b = curve.r,
                 φ0 = zero(T), span = curve.span, offset = zero(T))
             for curve in path if curve isa Arc] for path in paths]
    for (i, path) in enumerate(paths)
        append!(arcs[i],
            [(at = curve.at, a = curve.a, b = curve.b,
                 φ0 = curve.φ0, span = curve.span, offset = curve.offset)
             for curve in path if curve isa EllipticArc])
    end
    for arc in Iterators.flatten(arcs)
        start = arc.at.φ + atan(arc.a * sin(arc.φ0), arc.b * cos(arc.φ0))
        finish = arc.at.φ + atan(arc.a * sin(arc.φ0 + arc.span),
            arc.b * cos(arc.φ0 + arc.span))
        span = sign(arc.span) * mod(sign(arc.span) * (finish - start), T(2π))
        abs(nominal(arc.span)) >= 2π && (span = sign(arc.span) * T(2π))
        nominal(arc.span) < 0 && (start += T(π))
        low = start + (nominal(span) < 0 ? span : zero(T))
        high = start + (nominal(span) > 0 ? span : zero(T))
        for angle in (low, high)
            radius = support(A, angle)
            tolerance = sqrt(eps(float(nominal(radius)))) *
                        max(one(nominal(radius)), nominal(radius))
            nominal(support(B, angle)) <= nominal(radius) + tolerance || return false
        end
        while true
            left = low + (high - low) / 3
            right = high - (high - low) / 3
            (nominal(left) == nominal(low) || nominal(right) == nominal(high) ||
             nominal(left) == nominal(right)) && break
            if nominal(support(A, left) - support(B, left)) <
               nominal(support(A, right) - support(B, right))
                high = right
            else
                low = left
            end
        end
        φ = argmin(angle -> nominal(support(A, angle) - support(B, angle)),
            (low, (low + high) / 2, high))
        radius = support(A, φ)
        tolerance = sqrt(eps(float(nominal(radius)))) *
                    max(one(nominal(radius)), nominal(radius))
        nominal(support(B, φ)) <= nominal(radius) + tolerance || return false
    end
    return true
end

"""
$(TYPEDSIGNATURES)

Return whether ``B⊆A.outer`` and `B` lies apart from every hole.
A circular hole of radius ``r`` needs ``distance(B-c)+tol≥r``.
Another hole needs ``support(hole,ψ)≤-support(B,ψ+π)+tol`` in the outward normal of a
segment of the hole, or in the reversed normal of a segment of `B`.
Lengths and tolerances are in \\[m\\], angles in \\[rad\\].
The test is exact for polygonal outlines. Only segments give directions, so the method
refuses a content beside an arc of a hole other than a circle.
"""
function Base.issubset(B::AbstractShape, A::PerforatedShape)
    B ⊆ A.outer || return false
    for hole in [hole for hole in A.holes if hole isa Circle]
        radius = hole.r
        tolerance = sqrt(eps(float(nominal(radius)))) *
                    max(one(nominal(radius)), nominal(radius))
        nominal(distance(resolve(Placement(-hole.at.x, -hole.at.y), B))) + tolerance >=
        nominal(radius) || return false
    end
    for hole in [hole for hole in A.holes if !(hole isa Circle)]
        directions = [atan(-(edge.q[1] - edge.p[1]), edge.q[2] - edge.p[2])
                      for edge in Contour(hole).outer if edge isa Segment]
        append!(directions,
            [atan(-(edge.q[1] - edge.p[1]), edge.q[2] - edge.p[2]) + π
             for edge in Contour(B).outer if edge isa Segment])
        apart = false
        for ψ in directions
            radius = support(B, ψ + π)
            tolerance = sqrt(eps(float(nominal(radius)))) *
                        max(one(nominal(radius)), nominal(radius))
            apart |= nominal(support(hole, ψ)) <= -nominal(radius) + tolerance
        end
        apart || return false
    end
    return true
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
