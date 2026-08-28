"""
$(TYPEDEF)

Supertype for resolved cross-sectional primitives.

Each primitive contains its absolute pose in the completed cable coordinate
system. Relative placement remains part of the authoritative `Group`,
`Assembly`, or `Enclosure` declaration and is composed during [`build`](@ref).
"""
abstract type AbstractPrimitive{T <: Real} end

Base.eltype(::AbstractPrimitive{T}) where {T} = T
Base.eltype(::Type{<:AbstractPrimitive{T}}) where {T} = T

"Explicit initial state for resolving the first member of an outward stack."
struct EmptyBoundary end

"Resolved solid circle."
struct Disk{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    r::T
    at::P
end

"Resolved rectangle."
struct Rectangle{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    w::T
    h::T
    at::P
end

"Resolved ellipse."
struct Ellipse{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    a::T
    b::T
    at::P
end

"Resolved circular or annular sector."
struct Sector{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    ri::T
    ro::T
    φ0::T
    span::T
    at::P
end

"Resolved circular annulus."
struct Annulus{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    ri::T
    ro::T
    at::P
end

"Resolved polygon."
struct Polygon{T <: Real, V <: Tuple, P <: Pose2{T}} <: AbstractPrimitive{T}
    points::V
    at::P
end

function Disk(r::Real, at::Pose2)
    T = promote_type(typeof(r), eltype(at))
    return Disk{T, Pose2{T}}(convert(T, r), convert(Pose2{T}, at))
end
function Rectangle(w::Real, h::Real, at::Pose2)
    T = promote_type(typeof(w), typeof(h), eltype(at))
    return Rectangle{T, Pose2{T}}(
        convert(T, w), convert(T, h), convert(Pose2{T}, at))
end
function Ellipse(a::Real, b::Real, at::Pose2)
    T = promote_type(typeof(a), typeof(b), eltype(at))
    return Ellipse{T, Pose2{T}}(
        convert(T, a), convert(T, b), convert(Pose2{T}, at))
end
function Sector(ri::Real, ro::Real, φ0::Real, span::Real, at::Pose2)
    T = promote_type(
        typeof(ri), typeof(ro), typeof(φ0), typeof(span), eltype(at))
    return Sector{T, Pose2{T}}(
        convert(T, ri), convert(T, ro), convert(T, φ0), convert(T, span),
        convert(Pose2{T}, at))
end
function Annulus(ri::Real, ro::Real, at::Pose2)
    T = promote_type(typeof(ri), typeof(ro), eltype(at))
    return Annulus{T, Pose2{T}}(
        convert(T, ri), convert(T, ro), convert(Pose2{T}, at))
end
function _polygon(points::Tuple, at::Pose2)
    coordinate_types = (
        (typeof(coordinate) for point in points for coordinate in point)...,
        eltype(at)
    )
    T = promote_type(coordinate_types...)
    vertices = map(point -> (convert(T, point[1]), convert(T, point[2])), points)
    return Polygon{T, typeof(vertices), Pose2{T}}(
        vertices, convert(Pose2{T}, at))
end

_origin(::Type{T}) where {T <: Real} = Pose2(zero(T), zero(T), zero(T))

Disk(r::T) where {T <: Real} = Disk(r, _origin(T))
Rectangle(w::T, h::T) where {T <: Real} = Rectangle(w, h, _origin(T))
Ellipse(a::T, b::T) where {T <: Real} = Ellipse(a, b, _origin(T))
Sector(ri::T, ro::T, φ0::T, span::T) where {T <: Real} =
    Sector(ri, ro, φ0, span, _origin(T))
Annulus(ri::T, ro::T) where {T <: Real} = Annulus(ri, ro, _origin(T))
function Polygon(points::V) where {V <: Tuple}
    T = typeof(first(points)[1])
    return Polygon{T, V, Pose2{T}}(points, _origin(T))
end

"Return the outer resolved boundary of `primitive`."
function boundary end

"Return cross-sectional area \\[m²\\]."
function area end

"Return the absolute cross-sectional centroid as `(x, y)` \\[m\\]."
function centroid end

"Return the absolute directional support coordinate at angle `φ` \\[m\\]."
function support end

"Return the inner radius of radial geometry \\[m\\]."
function r_in end

"Return the outer radius of radial geometry \\[m\\]."
function r_ex end

"Return the radial thickness of radial geometry \\[m\\]."
function thickness end

@inline _local_centroid(::Union{Disk, Rectangle, Ellipse, Annulus}) = (0, 0)
function _local_centroid(primitive::Sector)
    full = isapprox(primitive.span, oftype(primitive.span, 2π))
    full && return (zero(primitive.span), zero(primitive.span))
    radial = 4 * sin(primitive.span / 2) *
             (primitive.ro^3 - primitive.ri^3) /
             (3 * primitive.span * (primitive.ro^2 - primitive.ri^2))
    angle = primitive.φ0 + primitive.span / 2
    return (radial * cos(angle), radial * sin(angle))
end
function _local_centroid(primitive::Polygon)
    signed_twice_area = zero(eltype(primitive))
    xmoment = zero(eltype(primitive))
    ymoment = zero(eltype(primitive))
    for index in eachindex(primitive.points)
        next = mod1(index + 1, length(primitive.points))
        left = primitive.points[index]
        right = primitive.points[next]
        cross = left[1] * right[2] - right[1] * left[2]
        signed_twice_area += cross
        xmoment += (left[1] + right[1]) * cross
        ymoment += (left[2] + right[2]) * cross
    end
    return (
        xmoment / (3 * signed_twice_area),
        ymoment / (3 * signed_twice_area)
    )
end

function centroid(primitive::AbstractPrimitive)
    x, y = _local_centroid(primitive)
    at = primitive.at
    return (
        at.x + cos(at.φ) * x - sin(at.φ) * y,
        at.y + sin(at.φ) * x + cos(at.φ) * y
    )
end

_local_support(primitive::Disk, φ::Real) = primitive.r
_local_support(primitive::Annulus, φ::Real) = primitive.ro
_local_support(primitive::Rectangle, φ::Real) =
    abs(cos(φ)) * primitive.w / 2 + abs(sin(φ)) * primitive.h / 2
_local_support(primitive::Ellipse, φ::Real) =
    hypot(primitive.a * cos(φ), primitive.b * sin(φ))
function _local_support(primitive::Sector, φ::Real)
    candidates = (primitive.φ0, primitive.φ0 + primitive.span, φ, φ + π)
    values = map(candidates) do angle
        relative = mod(angle - primitive.φ0, oftype(primitive.span, 2π))
        relative <= primitive.span || return -oftype(primitive.span, Inf)
        direction = cos(angle - φ)
        radius = direction >= zero(direction) ? primitive.ro : primitive.ri
        return radius * direction
    end
    return maximum(values)
end
_local_support(primitive::Polygon, φ::Real) = maximum(
    point[1] * cos(φ) + point[2] * sin(φ) for point in primitive.points
)

function support(primitive::AbstractPrimitive, φ::Real)
    at = primitive.at
    return at.x * cos(φ) + at.y * sin(φ) +
           _local_support(primitive, φ - at.φ)
end

_local_extent(primitive::Disk) = primitive.r
_local_extent(primitive::Annulus) = primitive.ro
_local_extent(primitive::Rectangle) = hypot(primitive.w, primitive.h) / 2
_local_extent(primitive::Ellipse) = max(primitive.a, primitive.b)
_local_extent(primitive::Sector) = primitive.ro
_local_extent(primitive::Polygon) = maximum(hypot(point...) for point in primitive.points)
support(primitive::AbstractPrimitive) = hypot(primitive.at.x, primitive.at.y) +
                                        _local_extent(primitive)

boundary(primitive::Union{Disk, Rectangle, Ellipse, Sector, Polygon}) = primitive
boundary(primitive::Annulus) = Disk(primitive.ro, primitive.at)

area(primitive::Disk) = π * primitive.r^2
area(primitive::Rectangle) = primitive.w * primitive.h
area(primitive::Ellipse) = π * primitive.a * primitive.b
area(primitive::Sector) =
    primitive.span * (primitive.ro^2 - primitive.ri^2) / 2
area(primitive::Annulus) = π * (primitive.ro^2 - primitive.ri^2)
function area(primitive::Polygon)
    twice = sum(eachindex(primitive.points)) do index
        next = mod1(index + 1, length(primitive.points))
        primitive.points[index][1] * primitive.points[next][2] -
            primitive.points[next][1] * primitive.points[index][2]
    end
    return abs(twice) / 2
end

r_in(primitive::Disk) = zero(primitive.r)
r_in(primitive::Annulus) = primitive.ri
r_in(primitive::Sector) = primitive.ri
r_ex(primitive::Disk) = primitive.r
r_ex(primitive::Annulus) = primitive.ro
r_ex(primitive::Sector) = primitive.ro
thickness(primitive::Disk) = primitive.r
thickness(primitive::Annulus) = primitive.ro - primitive.ri
thickness(primitive::Sector) = primitive.ro - primitive.ri

_with_pose(primitive::Disk, at::Pose2) = Disk(primitive.r, at)
_with_pose(primitive::Rectangle, at::Pose2) = Rectangle(primitive.w, primitive.h, at)
_with_pose(primitive::Ellipse, at::Pose2) = Ellipse(primitive.a, primitive.b, at)
_with_pose(primitive::Sector, at::Pose2) =
    Sector(primitive.ri, primitive.ro, primitive.φ0, primitive.span, at)
_with_pose(primitive::Annulus, at::Pose2) = Annulus(primitive.ri, primitive.ro, at)
_with_pose(primitive::Polygon, at::Pose2) = _polygon(primitive.points, at)
