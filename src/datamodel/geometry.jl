"""
$(TYPEDEF)

Supertype of the exact cross-sections of the parts of a cable, in the floating type `T`.
"""
abstract type AbstractShape{T <: Real} end

"""
$(TYPEDEF)

Supertype of the cross-sections that a part declares. A primitive states its dimensions and
its placement.
"""
abstract type AbstractPrimitive{T <: Real} <: AbstractShape{T} end

"""
$(TYPEDEF)

Represent a circular cross-section, such as a solid conductor or a round wire.

$(TYPEDFIELDS)
"""
struct Circle{T <: Real} <: AbstractPrimitive{T}
    "Radius \\[m\\]."
    r::T
    "Placement of the center in the cable cross-section."
    at::Placement{T}
    function Circle{T}(r::T, at::Placement{T}) where {T <: Real}
        return validate(new{T}(r, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a circle of radius `r` \\[m\\] centered at the origin.
"""
Circle(r::Real) = Circle(r, Placement(zero(r), zero(r)))

"""
$(TYPEDSIGNATURES)

Construct a circle of radius `r` \\[m\\] at the placement `at`, in the type that both
promote to with `Float64`.

# Errors

- Throws `DomainError` when `r` is not positive and finite.
"""
function Circle(r::Real, at::Placement)
    T = promote_type(Float64, typeof(r), eltype(at))
    return Circle{T}(convert(T, r), convert(Placement{T}, at))
end

function validate(circle::Circle)
    isfinite(circle.r) && circle.r > zero(circle.r) ||
        throw(DomainError(circle.r, "Circle.r must be positive and finite"))
    return circle
end

"""
$(TYPEDEF)

Represent an annular cross-section, such as a tubular conductor or an insulation layer.

$(TYPEDFIELDS)
"""
struct Annulus{T <: Real} <: AbstractPrimitive{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Placement of the center in the cable cross-section."
    at::Placement{T}
    function Annulus{T}(ri::T, ro::T, at::Placement{T}) where {T <: Real}
        return validate(new{T}(ri, ro, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an annulus between the radii `ri` and `ro` \\[m\\], centered at the origin.
"""
Annulus(ri::Real, ro::Real) = Annulus(ri, ro, Placement(zero(ri), zero(ro)))

"""
$(TYPEDSIGNATURES)

Construct an annulus between the radii `ri` and `ro` \\[m\\] at the placement `at`, in the
type that they promote to with `Float64`.

# Errors

- Throws `DomainError` when `ri` is negative or not finite, or when `ro` is not finite or does
  not exceed `ri`.
"""
function Annulus(ri::Real, ro::Real, at::Placement)
    T = promote_type(Float64, typeof(ri), typeof(ro), eltype(at))
    return Annulus{T}(convert(T, ri), convert(T, ro), convert(Placement{T}, at))
end

function validate(annulus::Annulus)
    isfinite(annulus.ri) && annulus.ri >= zero(annulus.ri) ||
        throw(DomainError(annulus.ri, "Annulus.ri must be nonnegative and finite"))
    isfinite(annulus.ro) && annulus.ro > annulus.ri ||
        throw(DomainError(annulus.ro, "Annulus.ro must be finite and exceed Annulus.ri"))
    return annulus
end

# An annulus that follows a part begins at or outside the part's outer radius.
function validate(annulus::Annulus, inner::Union{Circle, Annulus})
    radius = r_ex(inner)
    tolerance = sqrt(eps(float(radius))) * max(one(radius), radius)
    annulus.ri + tolerance >= radius || throw(DomainError(annulus.ri,
        "Annulus.ri must not be less than $radius m, the outer radius of the part inside it"))
    return annulus
end

"""
$(TYPEDEF)

Represent a rectangular cross-section: a container, or a strand before a circular pattern
bends it onto its course. The width of a strand lies along the course and its height
across it.

$(TYPEDFIELDS)
"""
struct Rectangle{T <: Real} <: AbstractPrimitive{T}
    "Width along the local x-axis \\[m\\]."
    w::T
    "Height along the local y-axis \\[m\\]."
    h::T
    "Placement of the center in the cable cross-section."
    at::Placement{T}
    function Rectangle{T}(w::T, h::T, at::Placement{T}) where {T <: Real}
        return validate(new{T}(w, h, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a rectangle of width `w` and height `h` \\[m\\] centered at the origin, in the
type that they promote to with `Float64`.

# Errors

- Throws `DomainError` when `w` or `h` is not positive and finite.
"""
Rectangle(w::Real, h::Real) = Rectangle(w, h, Placement(zero(w), zero(h)))

"""
$(TYPEDSIGNATURES)

Construct a rectangle of width `w` and height `h` \\[m\\] at the placement `at`, in the
type that they promote to with `Float64`.

# Errors

- Throws `DomainError` when `w` or `h` is not positive and finite.
"""
function Rectangle(w::Real, h::Real, at::Placement)
    T = promote_type(Float64, typeof(w), typeof(h), eltype(at))
    return Rectangle{T}(convert(T, w), convert(T, h), convert(Placement{T}, at))
end

function validate(rectangle::Rectangle)
    isfinite(rectangle.w) && rectangle.w > zero(rectangle.w) ||
        throw(DomainError(rectangle.w, "Rectangle.w must be positive and finite"))
    isfinite(rectangle.h) && rectangle.h > zero(rectangle.h) ||
        throw(DomainError(rectangle.h, "Rectangle.h must be positive and finite"))
    return rectangle
end

"""
$(TYPEDEF)

Declare a distance outward from the outline of the part inside, which the offset resolves
against. As the cross-section of a region, it is a layer of that thickness around the
outline, such as an insulation, a semiconducting screen or a jacket. As the radius of a
[`CircularPattern`](@ref), it is the clearance between the outline and the members.

$(TYPEDFIELDS)
"""
struct Offset{T <: Real}
    "Distance outward from the outline \\[m\\]."
    t::T
    function Offset{T}(t::T) where {T <: Real}
        return validate(new{T}(t))
    end
end

"""
$(TYPEDSIGNATURES)

Declare the distance `t` \\[m\\] outward from the outline of the part inside.

# Errors

- Throws `DomainError` when `t` is negative or not finite.
"""
function Offset(t::Real)
    T = promote_type(Float64, typeof(t))
    return Offset{T}(convert(T, t))
end

function validate(offset::Offset)
    isfinite(offset.t) && offset.t >= zero(offset.t) ||
        throw(DomainError(offset.t, "Offset.t must be nonnegative and finite"))
    return offset
end

"""
$(TYPEDEF)

Represent a sector: the region between the inner radius ``r_i`` and the outer
radius ``r_o`` over the angular span ``\\Delta\\phi``, symmetric about the x-axis of its
placement. The sides meet at `apex` and each convex corner has radius `fillet`
\\[m\\]. A circular pattern bends a rectangular strand of width ``w`` and height ``h``
into an annular sector of the same area:

```math
\\frac{\\Delta\\phi}{2}\\left(r_o^2 - r_i^2\\right) = w h.
```

$(TYPEDFIELDS)
"""
struct Sector{T <: Real} <: AbstractPrimitive{T}
    "Inner radius \\[m\\], zero when the sides meet."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Angular span \\[rad\\]."
    span::T
    "Distance \\[m\\] from the center to the meeting point of the sides, zero if radial."
    apex::T
    "Radius \\[m\\] of the arc that rounds each convex corner, zero for sharp corners."
    fillet::T
    "Placement of the center of curvature. Its x-axis is the symmetry axis of the sector."
    at::Placement{T}
    function Sector{T}(ri::T, ro::T, span::T, apex::T, fillet::T,
            at::Placement{T}) where {T <: Real}
        return validate(new{T}(ri, ro, span, apex, fillet, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a sector from its radii, apex distance and fillet radius \\[m\\], angular
span \\[rad\\] and placement, in the type that they promote to with `Float64`.

# Errors

- Throws `DomainError` for the first field outside the bounds of a sector, or for a
  fillet that does not fit.
"""
function Sector(ri::Real, ro::Real, span::Real, apex::Real, fillet::Real, at::Placement)
    T = promote_type(Float64, typeof(ri), typeof(ro), typeof(span), typeof(apex),
        typeof(fillet), eltype(at))
    return Sector{T}(convert(T, ri), convert(T, ro), convert(T, span), convert(T, apex),
        convert(T, fillet), convert(Placement{T}, at))
end

function validate(sector::Sector)
    isfinite(nominal(sector.ri)) && nominal(sector.ri) >= zero(nominal(sector.ri)) ||
        throw(DomainError(sector.ri, "Sector.ri must be nonnegative and finite"))
    isfinite(nominal(sector.ro)) && nominal(sector.ro) > nominal(sector.ri) ||
        throw(DomainError(sector.ro, "Sector.ro must be finite and exceed Sector.ri"))
    zero(nominal(sector.span)) < nominal(sector.span) <= 2π ||
        throw(DomainError(sector.span, "Sector.span must lie in (0, 2π]"))
    isfinite(nominal(sector.apex)) && nominal(sector.apex) >= zero(nominal(sector.apex)) &&
    (nominal(sector.span) >= π ||
     nominal(sector.apex) < nominal(sector.ro) * cos(nominal(sector.span) / 2)) ||
        throw(DomainError(sector.apex,
            "Sector.apex must be nonnegative and finite, and below ro cos(span / 2) when span < π"))
    isfinite(nominal(sector.fillet)) &&
    nominal(sector.fillet) >= zero(nominal(sector.fillet)) ||
        throw(DomainError(sector.fillet, "Sector.fillet must be nonnegative and finite"))
    if nominal(sector.fillet) > zero(nominal(sector.fillet))
        nominal(sector.apex) * sin(nominal(sector.span) / 2) + 2nominal(sector.fillet) <
        nominal(sector.ro) ||
            throw(DomainError(sector.fillet,
                "Sector.fillet must satisfy apex sin(span / 2) + 2 fillet < ro"))
        if iszero(nominal(sector.ri))
            nominal(sector.fillet) * (csc(nominal(sector.span) / 2) - 1) +
            nominal(sector.apex) <
            nominal(sector.ro) - 2nominal(sector.fillet) ||
                throw(DomainError(sector.fillet,
                    "Sector.fillet must satisfy fillet (csc(span / 2) - 1) + apex < ro - 2 fillet"))
        end
    end
    return sector
end

"""
$(TYPEDEF)

Represent an ellipse with semi-axes `a` and `b` along the axes of its placement.

$(TYPEDFIELDS)
"""
struct Ellipse{T <: Real} <: AbstractPrimitive{T}
    "Semi-axis along the local x-axis \\[m\\]."
    a::T
    "Semi-axis along the local y-axis \\[m\\]."
    b::T
    "Placement of the center and axes."
    at::Placement{T}
    @doc """
    $(TYPEDSIGNATURES)

    Construct an ellipse through `validate`, with semi-axes in \\[m\\].
    """
    function Ellipse{T}(a::T, b::T, at::Placement{T}) where {T <: Real}
        return validate(new{T}(a, b, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct an ellipse at the origin with semi-axes `a` and `b` \\[m\\], promoted with `Float64`.
"""
function Ellipse(a::Real, b::Real)
    T = promote_type(Float64, typeof(a), typeof(b))
    return Ellipse{T}(convert(T, a), convert(T, b), Placement(zero(T), zero(T)))
end

"""
$(TYPEDSIGNATURES)

Construct an ellipse at `at` with semi-axes `a` and `b` \\[m\\], promoted with `Float64`.
"""
function Ellipse(a::Real, b::Real, at::Placement)
    T = promote_type(Float64, typeof(a), typeof(b), eltype(at))
    return Ellipse{T}(convert(T, a), convert(T, b), convert(Placement{T}, at))
end

"""
$(TYPEDSIGNATURES)

Check ``a>0`` and ``b>0`` finite semi-axes \\[m\\], and return `ellipse`.
Throw `DomainError` for a semi-axis outside these bounds.
"""
function validate(ellipse::Ellipse)
    isfinite(nominal(ellipse.a)) && nominal(ellipse.a) > zero(nominal(ellipse.a)) ||
        throw(DomainError(ellipse.a, "Ellipse.a must be positive and finite"))
    isfinite(nominal(ellipse.b)) && nominal(ellipse.b) > zero(nominal(ellipse.b)) ||
        throw(DomainError(ellipse.b, "Ellipse.b must be positive and finite"))
    return ellipse
end

"""
$(TYPEDEF)

Represent a simple polygon with counterclockwise vertices in its placement's frame.

$(TYPEDFIELDS)
"""
struct Polygon{T <: Real} <: AbstractPrimitive{T}
    "Counterclockwise vertices in the local frame \\[m\\]."
    vertices::Vector{Tuple{T, T}}
    "Placement of the local frame."
    at::Placement{T}
    @doc """
    $(TYPEDSIGNATURES)

    Construct a polygon through `validate`, with vertex coordinates in \\[m\\].
    """
    function Polygon{T}(vertices::Vector{Tuple{T, T}}, at::Placement{T}) where {T <: Real}
        return validate(new{T}(vertices, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a polygon at the origin from `vertices` \\[m\\], promoted with `Float64`.
"""
function Polygon(vertices::AbstractVector{<:Tuple{Real, Real}})
    T = promote_type(Float64, (typeof(x) for vertex in vertices for x in vertex)...)
    return Polygon{T}(Tuple{T, T}[convert.(T, vertex) for vertex in vertices],
        Placement(zero(T), zero(T)))
end

"""
$(TYPEDSIGNATURES)

Construct a polygon from local `vertices` \\[m\\] at `at`, promoted with `Float64`.
"""
function Polygon(vertices::AbstractVector{<:Tuple{Real, Real}}, at::Placement)
    T = promote_type(Float64, eltype(at), (typeof(x) for vertex in vertices
    for x in vertex)...)
    return Polygon{T}(Tuple{T, T}[convert.(T, vertex) for vertex in vertices],
        convert(Placement{T}, at))
end

"""
$(TYPEDSIGNATURES)

Check for at least three finite vertices, a positive `area` \\[m²\\], and edges that
meet only at adjacent shared ends. Segment intersections use the signs of
``(Q-P)\\times(R-P)``. Return `polygon`, or throw `DomainError`.
"""
function validate(polygon::Polygon)
    vertices = polygon.vertices
    n = length(vertices)
    n >= 3 || throw(DomainError(vertices, "Polygon needs at least three vertices"))
    all(isfinite(nominal(x)) for vertex in vertices for x in vertex) ||
        throw(DomainError(vertices, "Polygon coordinates must be finite"))
    nominal(area(polygon)) > 0 ||
        throw(DomainError(vertices, "Polygon signed area must be positive"))
    for i in 1:n
        P, Q = vertices[i], vertices[mod1(i + 1, n)]
        any(!iszero(nominal(x)) for x in Q .- P) ||
            throw(DomainError(vertices, "Polygon edges must have distinct ends"))
        for j in (i + 1):n
            R, S = vertices[j], vertices[mod1(j + 1, n)]
            if j == i + 1 || (i == 1 && j == n)
                a, b, c = j == i + 1 ? (P, Q, S) : (R, S, Q)
                cross = (b[1] - a[1]) * (c[2] - b[2]) -
                        (b[2] - a[2]) * (c[1] - b[1])
                dot = sum((b .- a) .* (c .- b))
                iszero(nominal(cross)) && nominal(dot) <= zero(nominal(dot)) &&
                    throw(DomainError(vertices, "Polygon adjacent edges must not overlap"))
                continue
            end
            crosses = ((Q[1] - P[1]) * (R[2] - P[2]) - (Q[2] - P[2]) * (R[1] - P[1]),
                (Q[1] - P[1]) * (S[2] - P[2]) - (Q[2] - P[2]) * (S[1] - P[1]),
                (S[1] - R[1]) * (P[2] - R[2]) - (S[2] - R[2]) * (P[1] - R[1]),
                (S[1] - R[1]) * (Q[2] - R[2]) - (S[2] - R[2]) * (Q[1] - R[1]))
            if min(nominal(P[1]), nominal(Q[1])) <= max(nominal(R[1]), nominal(S[1])) &&
               min(nominal(R[1]), nominal(S[1])) <= max(nominal(P[1]), nominal(Q[1])) &&
               min(nominal(P[2]), nominal(Q[2])) <= max(nominal(R[2]), nominal(S[2])) &&
               min(nominal(R[2]), nominal(S[2])) <= max(nominal(P[2]), nominal(Q[2])) &&
               nominal(crosses[1]) * nominal(crosses[2]) <= 0 &&
               nominal(crosses[3]) * nominal(crosses[4]) <= 0
                throw(DomainError(vertices, "Polygon edges must not cross or touch"))
            end
        end
    end
    return polygon
end

"""
$(TYPEDEF)

Represent a cross-section with holes: the outer cross-section without the cross-sections
inside it, such as the fill between the wires of a course.

$(TYPEDFIELDS)
"""
struct PerforatedShape{
    T <: Real, O <: AbstractShape{T}, H <: AbstractVector{<:AbstractShape}} <:
       AbstractShape{T}
    "Outer cross-section."
    outer::O
    "Outlines removed from the outer one. They lie inside it, and each lies outside the others."
    holes::H
    function PerforatedShape(outer::O,
            holes::H) where {
            T <: Real, O <: AbstractShape{T}, H <: AbstractVector{<:AbstractShape}}
        return new{T, O, H}(outer, holes)
    end
end

"""
$(TYPEDEF)

Represent the outline at the center of a cable, before its first part. The first part keeps
its declared cross-section.
"""
struct EmptyBoundary
    EmptyBoundary() = new()
end

Base.eltype(::AbstractShape{T}) where {T} = T

"""
    issubset(B::AbstractShape, A::AbstractShape)

Return whether the outline `B` lies within the region of `A`, written `B ⊆ A`.
For a circle or annulus `A`, translate `B` to the frame centered on `A` as
`C = resolve(Placement(-A.at.x, -A.at.y), B)`. Compare ``support(C)≤r_ex(A)+tol``
and ``distance(C)+tol≥r_in(A)``, with the tolerance computed from `r_ex(A)`.

For a rectangle, polygon, ellipse or arc polygon `A`, compare
``support(B,φ)≤support(A,φ)+tol`` at the outward normal at each end of every curve of
both contours. A sector takes this test through the arc polygon of its contour. Between
two consecutive normals both supports are smooth. A ternary search finds the least
``g(φ)=support(A,φ)-support(B,φ)`` between them, and the test compares it too. This is
exact for convex `A`.

The arithmetic uses the promoted coordinate type. Decisions and the stop criterion use
`nominal` values. The support tolerance is ``√eps(float(x))max(1,x)`` for the
support ``x`` compared \\[m\\].

A sector also requires ``distance(C)+tol≥A.ri``, with
`C = resolve(Placement(-A.at.x, -A.at.y), B)` and the tolerance computed from `A.ri`.
For `A.ri = 0`, this condition holds at once. An annular sector is its convex hull
without the disk of radius `A.ri` about its center, so the test is exact for it.

An arc polygon uses the support of its convex hull.
A content in a concave pocket can pass. The design's cover check tests it.
A perforated `B` uses the directional support of its outer boundary.
"""
function Base.issubset(B::AbstractShape, A::Union{Circle, Annulus})
    C = resolve(Placement(-A.at.x, -A.at.y), B)
    radius = r_ex(A)
    tolerance = sqrt(eps(float(nominal(radius)))) *
                max(one(nominal(radius)), nominal(radius))
    return nominal(support(C)) <= nominal(radius) + tolerance &&
           nominal(distance(C)) + tolerance >= nominal(r_in(A))
end

function Base.issubset(B::AbstractShape, A::Sector)
    B ⊆ ArcPolygon(Contour(A).outer) || return false
    iszero(nominal(A.ri)) && return true
    C = resolve(Placement(-A.at.x, -A.at.y), B)
    radius = A.ri
    tolerance = sqrt(eps(float(nominal(radius)))) *
                max(one(nominal(radius)), nominal(radius))
    return nominal(distance(C)) + tolerance >= nominal(radius)
end

"""
    area(shape)

Return the cross-sectional area of `shape` \\[m²\\]:

```math
A_\\mathrm{circle} = \\pi r^2, \\qquad
A_\\mathrm{annulus} = \\pi \\left(r_o^2 - r_i^2\\right), \\qquad
A_\\mathrm{rectangle} = w h, \\qquad
A_\\mathrm{sector} = \\frac{\\Delta\\phi}{2} \\left(r_o^2 - r_i^2\\right),
```

where ``r`` is the radius of a circle, ``r_i`` and ``r_o`` are inner and outer radii, ``w``
and ``h`` are the width and height of a rectangle \\[m\\], and ``\\Delta\\phi`` \\[rad\\] is
the angular span of a sector. A [`PerforatedShape`](@ref) has the area of its outer
cross-section less the areas of its holes. A sector with displaced sides or fillets, and
an arc polygon, take the area of their contour by Green's theorem.

An ellipse has ``A=πab``. A polygon uses
``A=\\sum_i(x_i y_{i+1}-y_i x_{i+1})/2``.
A layer subtracts the wrapped outline's area and gives ``Pt+πt²`` for a convex outline.
"""
area(circle::Circle) = π * circle.r^2
area(annulus::Annulus) = π * (annulus.ro^2 - annulus.ri^2)
area(rectangle::Rectangle) = rectangle.w * rectangle.h
function area(sector::Sector)
    if iszero(sector.apex) && iszero(nominal(sector.fillet))
        return sector.span * (sector.ro^2 - sector.ri^2) / 2
    end
    return area(Contour(sector))
end
function area(shape::PerforatedShape{T}) where {T}
    return area(shape.outer) - (sum(area, shape.holes; init = zero(T))::T)
end

"""
    centroid(shape)

Return the centroid `(x, y)` of `shape` in the cable cross-section \\[m\\]. For an annular
sector, the centroid lies on its symmetry axis at the distance

```math
d = \\frac{4 \\sin(\\Delta\\phi / 2) \\left(r_o^3 - r_i^3\\right)}
    {3 \\Delta\\phi \\left(r_o^2 - r_i^2\\right)}
```

from its center of curvature, where ``\\Delta\\phi`` is its angular span and ``r_i`` and ``r_o``
are its inner and outer radii. A circle and an annulus have their centroid at their center. A
[`PerforatedShape`](@ref) of area ``A`` has the centroid

```math
\\mathbf c = \\frac{A_o \\mathbf c_o - \\sum_h A_h \\mathbf c_h}{A},
```

where ``A_o`` and ``\\mathbf c_o`` are the area and centroid of its outer cross-section, and
``A_h`` and ``\\mathbf c_h`` those of each hole. A rectangle's centroid is its center.
A rounded sector or an arc polygon takes the centroid of its contour.

An ellipse's centroid is its center. A polygon uses
``c=\\sum_i(v_i+v_{i+1})(v_i\\times v_{i+1})/(6A)``, then applies its placement.
The layer's centroid subtracts the wrapped outline's moment and divides by its area.
"""
centroid(shape::Union{Circle, Annulus}) = (shape.at.x, shape.at.y)

function centroid(sector::Sector)
    if !(iszero(sector.apex) && iszero(nominal(sector.fillet)))
        return centroid(Contour(sector))
    end
    d = 4sin(sector.span / 2) * (sector.ro^3 - sector.ri^3) /
        (3sector.span * (sector.ro^2 - sector.ri^2))
    return sector.at((d, zero(d)))
end

function centroid(shape::PerforatedShape{T}) where {T}
    x, y = area(shape.outer) .* centroid(shape.outer)
    for hole in shape.holes
        x, y = (x, y) .- (area(hole)::T) .* (centroid(hole)::Tuple{T, T})
    end
    return (x, y) ./ area(shape)
end

"""
    support(shape)

Return the radius of the smallest circle about the origin that contains `shape` \\[m\\],
the distance of its farthest point from the origin. For a circle or an annulus, it is

```math
S = \\sqrt{x_c^2 + y_c^2} + r_o,
```

where ``(x_c, y_c)`` is the center of `shape` and ``r_o`` its outer radius \\[m\\]. The
support of a [`PerforatedShape`](@ref) is that of its outer cross-section. Of a polygon,
it is the farthest vertex. A rectangle, a sector or an arc polygon takes the support of its
contour. An ellipse takes the point of its curve farthest from the origin.

`support(shape, φ)` returns the support function ``\\max(x\\cos φ+y\\sin φ)`` \\[m\\]
of `shape` in the direction `φ` \\[rad\\]. For a circle or an annulus, it is the
projection of the center plus the outer radius. A rectangle takes its farthest corner. A
sector or an arc polygon takes the support of its contour. The ellipse's directional
support is ``c⋅u+√(a²cos²ψ+b²sin²ψ)``, ``ψ=φ-at.φ``. A layer has the support of its
parallel boundary.

Holes leave the support unchanged, so for an annulus or a sector it is the support of the
convex hull.
"""
support(circle::Circle) = distance(circle.at) + circle.r
support(annulus::Annulus) = distance(annulus.at) + annulus.ro
support(sector::Sector) = support(Contour(sector))
support(shape::PerforatedShape) = support(shape.outer)

"""
    distance(shape)

Return the radius of a circle about the origin that `shape` does not enter \\[m\\]:

```math
\\max(0, d - r_o, r_i - d),
```

where ``d`` is the distance from the origin to the center of `shape`, and ``r_i`` and
``r_o`` are its inner and outer radii. The inner radius of a circle is zero, and the center of
a sector is its center of curvature. For a sector, ``r_i`` is the nearest point of its
inner side: `ri` for an inner arc, `apex` for a sharp tip and
``apex+fillet/\\sin(span/2)-fillet`` for a rounded tip.

A rectangle, a polygon, an arc polygon or an ellipse takes the distance of its contour. A
layer keeps its outer boundary's distance bound.
"""
distance(circle::Circle) = max(zero(circle.r), distance(circle.at) - circle.r)

function distance(shape::Annulus)
    offset = distance(shape.at)
    return max(zero(offset), offset - shape.ro, shape.ri - offset)
end

distance(shape::PerforatedShape) = distance(shape.outer)

"""
    boundary(shape)

Return the outline of `shape`. A circle, rectangle, sector, ellipse, polygon and arc polygon are their own outlines. The
outline of an annulus is the circle of its outer radius. The outline of a
[`PerforatedShape`](@ref) is the outline of its outer cross-section. The next part resolves
against the outline, and an enclosure removes it from its fill.
"""
function boundary end

boundary(circle::Circle) = circle
boundary(annulus::Annulus) = Circle(annulus.ro, annulus.at)
boundary(sector::Sector) = sector
boundary(shape::PerforatedShape) = boundary(shape.outer)

"""
    r_in(shape)

Return the inner radius of a circular shape \\[m\\]: zero for a circle.
"""
function r_in end

r_in(circle::Circle) = zero(circle.r)
r_in(annulus::Annulus) = annulus.ri
r_in(sector::Sector) = sector.ri

"""
    r_ex(shape)

Return the outer radius of a circular shape \\[m\\].
"""
function r_ex end

r_ex(circle::Circle) = circle.r
r_ex(annulus::Annulus) = annulus.ro
r_ex(sector::Sector) = sector.ro

# The first part of a cable keeps its declared cross-section.
resolve(::EmptyBoundary, primitive::AbstractPrimitive) = primitive

# An offset around a circle or an annulus is the annulus of its thickness at the same
# placement.
function resolve(inner::Union{Circle, Annulus}, offset::Offset)
    return Annulus(r_ex(inner), r_ex(inner) + offset.t, inner.at)
end

# A declared annulus keeps its radii and takes the placement of the part inside it.
function resolve(inner::Union{Circle, Annulus}, annulus::Annulus)
    return Annulus(validate(annulus, inner).ri, annulus.ro, inner.at)
end

# Each method composes `at` with the placement of the part that it places.
resolve(at::Placement, circle::Circle) = Circle(circle.r, at * circle.at)
function resolve(at::Placement, annulus::Annulus)
    return Annulus(annulus.ri, annulus.ro, at * annulus.at)
end
function resolve(at::Placement, sector::Sector)
    return Sector(
        sector.ri, sector.ro, sector.span, sector.apex, sector.fillet, at * sector.at)
end
function resolve(at::Placement, shape::PerforatedShape)
    outer = resolve(at, shape.outer)
    return PerforatedShape(outer,
        AbstractShape{eltype(outer)}[resolve(at, hole) for hole in shape.holes])
end

Base.map(f, circle::Circle) = Circle(f(circle.r), map(f, circle.at))
function Base.map(f, annulus::Annulus)
    return Annulus(f(annulus.ri), f(annulus.ro), map(f, annulus.at))
end
function Base.map(f, rectangle::Rectangle)
    w, h = f(rectangle.w), f(rectangle.h)
    at = map(f, rectangle.at)
    T = promote_type(Float64, typeof(w), typeof(h), eltype(at))
    return Rectangle{T}(convert(T, w), convert(T, h), convert(Placement{T}, at))
end
Base.map(f, offset::Offset) = Offset(f(offset.t))
function Base.map(f, sector::Sector)
    return Sector(f(sector.ri), f(sector.ro), f(sector.span), f(sector.apex),
        f(sector.fillet), map(f, sector.at))
end
function Base.map(f, shape::PerforatedShape)
    outer = map(f, shape.outer)
    return PerforatedShape(outer,
        AbstractShape{eltype(outer)}[map(f, hole) for hole in shape.holes])
end

centroid(rectangle::Rectangle) = (rectangle.at.x, rectangle.at.y)
support(rectangle::Rectangle) = support(Contour(rectangle))
boundary(rectangle::Rectangle) = rectangle
distance(rectangle::Rectangle) = distance(Contour(rectangle))

function distance(sector::Sector)
    offset = hypot(sector.at.x, sector.at.y)
    if iszero(sector.apex) && iszero(nominal(sector.fillet))
        return max(zero(offset), offset - sector.ro, sector.ri - offset)
    end
    inner = if nominal(sector.ri) > zero(nominal(sector.ri))
        sector.ri
    else
        sector.apex + sector.fillet / sin(sector.span / 2) - sector.fillet
    end
    return max(zero(offset), offset - sector.ro, inner - offset)
end

support(circle::Circle, φ::Real) = circle.at.x * cos(φ) + circle.at.y * sin(φ) + circle.r
function support(annulus::Annulus, φ::Real)
    annulus.at.x * cos(φ) + annulus.at.y * sin(φ) + annulus.ro
end
support(sector::Sector, φ::Real) = support(Contour(sector), φ)
support(rectangle::Rectangle, φ::Real) = support(Contour(rectangle), φ)
support(shape::PerforatedShape, φ::Real) = support(shape.outer, φ)

"""
    perimeter(shape)

Return the length of every outer and hole curve of `shape` \\[m\\]:

``P_\\mathrm{circle}=2πr``, ``P_\\mathrm{annulus}=2π(r_i+r_o)`` and
``P_\\mathrm{rectangle}=2(w+h)``, where the radii, width and height are in \\[m\\].
A sector, an arc polygon and a [`PerforatedShape`](@ref) sum the curves of their contour.
An ellipse has ``P=4aE(1-b²/a²)``, with ``E`` the complete elliptic integral of the
second kind.

The formula also holds for ``b>a``. Polygon perimeters sum their edges' lengths.
A layer sums the perimeters of its parallel boundary and its hole.
"""
perimeter(circle::Circle) = 2π * circle.r
perimeter(annulus::Annulus) = 2π * (annulus.ri + annulus.ro)
perimeter(rectangle::Rectangle) = 2 * (rectangle.w + rectangle.h)
perimeter(sector::Sector) = perimeter(Contour(sector))
perimeter(shape::PerforatedShape) = perimeter(Contour(shape))

"""
$(TYPEDSIGNATURES)

Return the rectangle placed by `at`, its placement
``\\mathrm{at}\\,\\mathrm{at}_\\mathrm{rectangle}``, in the type that its numbers promote
to with `Float64`. Dimensions and translations are in \\[m\\], rotations in \\[rad\\].
"""
function resolve(at::Placement, rectangle::Rectangle)
    placed = at * rectangle.at
    T = promote_type(Float64, eltype(rectangle), eltype(placed))
    return Rectangle{T}(convert(T, rectangle.w), convert(T, rectangle.h),
        convert(Placement{T}, placed))
end

"""
$(TYPEDSIGNATURES)

Resolve the outward offset `t` \\[m\\] of a sector by
``(r_i,r_o,a,f)\\mapsto(\\max(r_i-t,0),r_o+t,a-t/\\sin(α),f+t)``,
where ``α=span/2`` and all four lengths are in \\[m\\].
Keep its angular span and placement. Return the layer of thickness `t`,
``PerforatedShape(Sector′,[inner])``, whose `boundary` is the parallel.
"""
function resolve(inner::Sector, offset::Offset)
    validate(offset, inner)
    T = promote_type(eltype(inner), typeof(offset.t))
    return PerforatedShape(
        Sector(max(inner.ri - offset.t, zero(inner.ri)), inner.ro + offset.t,
            inner.span, inner.apex - offset.t / sin(inner.span / 2),
            inner.fillet + offset.t, inner.at),
        AbstractShape{T}[inner])
end

"""
$(TYPEDSIGNATURES)

Check that ``t\\leq apex\\sin(span/2)``, the largest offset thickness that fits
\\[m\\], and return `offset`. Throw `DomainError` above that limit.
"""
function validate(offset::Offset, sector::Sector)
    limit = sector.apex * sin(sector.span / 2)
    offset.t <= limit || throw(DomainError(offset.t,
        "Offset.t must not exceed $limit m, apex sin(span / 2), the thickness that fits"))
    return offset
end

area(ellipse::Ellipse) = π * ellipse.a * ellipse.b
function area(polygon::Polygon)
    vertices = polygon.vertices
    return sum(vertices[i][1] * vertices[mod1(i + 1, length(vertices))][2] -
               vertices[i][2] * vertices[mod1(i + 1, length(vertices))][1]
    for i in eachindex(vertices)) / 2
end
centroid(ellipse::Ellipse) = (ellipse.at.x, ellipse.at.y)
function centroid(polygon::Polygon{T}) where {T}
    x, y = zero(T), zero(T)
    vertices = polygon.vertices
    for i in eachindex(vertices)
        P, Q = vertices[i], vertices[mod1(i + 1, length(vertices))]
        cross = P[1] * Q[2] - P[2] * Q[1]
        x += (P[1] + Q[1]) * cross
        y += (P[2] + Q[2]) * cross
    end
    return polygon.at((x, y) ./ (6area(polygon)))
end
support(ellipse::Ellipse) = support(EllipticArc(ellipse.at, ellipse.a, ellipse.b, 0, 2π))
function support(ellipse::Ellipse, φ::Real)
    ψ = φ - ellipse.at.φ
    return ellipse.at.x * cos(φ) + ellipse.at.y * sin(φ) +
           sqrt(ellipse.a^2 * cos(ψ)^2 + ellipse.b^2 * sin(ψ)^2)
end
function support(polygon::Polygon)
    return maximum(hypot(polygon.at(vertex)...) for vertex in polygon.vertices)
end
function support(polygon::Polygon, φ::Real)
    return maximum(point[1] * cos(φ) + point[2] * sin(φ)
    for point in (polygon.at(vertex) for vertex in polygon.vertices))
end
distance(ellipse::Ellipse) = distance(Contour(ellipse))
distance(polygon::Polygon) = distance(Contour(polygon))
perimeter(ellipse::Ellipse) = 4ellipse.a * ellipe(1 - ellipse.b^2 / ellipse.a^2)
function perimeter(polygon::Polygon)
    return sum(hypot((polygon.vertices[mod1(i + 1, length(polygon.vertices))] .-
                      polygon.vertices[i])...) for i in eachindex(polygon.vertices))
end
boundary(ellipse::Ellipse) = ellipse
boundary(polygon::Polygon) = polygon

"""
$(TYPEDSIGNATURES)

Place the ellipse by ``at_{new}=at\\,at_{ellipse}`` and keep its semi-axes \\[m\\].
"""
function resolve(at::Placement, ellipse::Ellipse)
    placed = at * ellipse.at
    T = promote_type(Float64, eltype(ellipse), eltype(placed))
    return Ellipse{T}(convert(T, ellipse.a), convert(T, ellipse.b), convert(Placement{T}, placed))
end

"""
$(TYPEDSIGNATURES)

Place the polygon by ``at_{new}=at\\,at_{polygon}`` and keep its local vertices \\[m\\].
"""
function resolve(at::Placement, polygon::Polygon)
    return Polygon(polygon.vertices, at * polygon.at)
end

"""
$(TYPEDSIGNATURES)

Map the semi-axes \\[m\\] and placement of `ellipse` by ``x\\mapsto f(x)``.
"""
function Base.map(f, ellipse::Ellipse)
    return Ellipse(f(ellipse.a), f(ellipse.b), map(f, ellipse.at))
end

"""
$(TYPEDSIGNATURES)

Map the local vertex coordinates \\[m\\] and placement of `polygon` by ``x\\mapsto f(x)``.
"""
function Base.map(f, polygon::Polygon)
    return Polygon([map(f, vertex) for vertex in polygon.vertices], map(f, polygon.at))
end

"""
$(TYPEDSIGNATURES)

Resolve a distance ``t`` \\[m\\] outward over the outer contour of `shape`.
Return the layer of thickness `t`, ``PerforatedShape(parallel,[boundary(shape)])``.
The layer rests on the parallel of the contour.
Its `boundary` is the parallel.
Throw `ArgumentError` for a clockwise arc or corner:
"an offset follows a circle, an annulus, a sector, a rectangle, an ellipse or a convex polygon".
"""
function resolve(shape::AbstractShape, offset::Offset)
    path = Contour(shape).outer
    T = promote_type(eltype(shape), typeof(offset.t))
    normals = Tuple{T, T}[]
    indices = Int[]
    for (k, curve) in [(k, curve) for (k, curve) in enumerate(path) if curve isa Segment]
        push!(indices, k)
        dx, dy = curve.q .- curve.p
        φ = atan(-dx, dy)
        push!(normals, (φ, φ))
    end
    for (k, curve) in [(k, curve) for (k, curve) in enumerate(path) if curve isa Arc]
        nominal(curve.span) > 0 ||
            throw(ArgumentError("an offset follows a circle, an annulus, a sector, a rectangle, an ellipse or a convex polygon"))
        push!(indices, k)
        push!(normals, (curve.at.φ, curve.at.φ + curve.span))
    end
    for (k, curve) in
        [(k, curve) for (k, curve) in enumerate(path) if curve isa EllipticArc]
        nominal(curve.span) > 0 ||
            throw(ArgumentError("an offset follows a circle, an annulus, a sector, a rectangle, an ellipse or a convex polygon"))
        θ₁, θ₂ = curve.φ0, curve.φ0 + curve.span
        push!(indices, k)
        push!(normals,
            (atan(curve.a * sin(θ₁), curve.b * cos(θ₁)) + curve.at.φ,
                atan(curve.a * sin(θ₂), curve.b * cos(θ₂)) + curve.at.φ))
    end
    order = sortperm(indices)
    normals = normals[order]
    for k in eachindex(path)
        φ = normals[k][2]
        Δ = mod(normals[mod1(k + 1, length(path))][1] - φ + π, 2π) - π
        scale = max(one(nominal(φ)), abs(nominal(φ)),
            abs(nominal(normals[mod1(k + 1, length(path))][1])))
        tolerance = sqrt(eps(float(scale))) * max(one(scale), scale)
        nominal(Δ) >= -tolerance ||
            throw(ArgumentError("an offset follows a circle, an annulus, a sector, a rectangle, an ellipse or a convex polygon"))
    end
    return PerforatedShape(ArcPolygon(parallel(Contour(shape), offset.t).outer),
        AbstractShape{T}[boundary(shape)])
end
