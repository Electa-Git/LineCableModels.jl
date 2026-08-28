"""
$(TYPEDEF)

Supertype for resolved absolute cross-sectional geometry.
"""
abstract type AbstractShape{T <: Real} end

Base.eltype(::AbstractShape{T}) where {T} = T
Base.eltype(::Type{<:AbstractShape{T}}) where {T} = T

"Explicit initial state for resolving the first member of an outward stack."
struct EmptyBoundary end

"Resolved solid circle."
struct DiskShape{T <: Real} <: AbstractShape{T}
    r::T
end

"Resolved centered rectangle."
struct RectangleShape{T <: Real} <: AbstractShape{T}
    w::T
    h::T
end

"Resolved centered ellipse."
struct EllipseShape{T <: Real} <: AbstractShape{T}
    a::T
    b::T
end

"Resolved circular or annular sector."
struct SectorShape{T <: Real} <: AbstractShape{T}
    ri::T
    ro::T
    φ0::T
    span::T
end

"Resolved circular annulus."
struct AnnulusShape{T <: Real} <: AbstractShape{T}
    ri::T
    ro::T
end

"Resolved polygon."
struct PolygonShape{T <: Real, P <: Tuple} <: AbstractShape{T}
    points::P
end

function PolygonShape(points::P) where {P <: Tuple}
    T = typeof(first(points)[1])
    return PolygonShape{T, P}(points)
end

"Resolved shape placed in a parent frame."
struct PlacedShape{T <: Real, S <: AbstractShape, P} <: AbstractShape{T}
    shape::S
    at::P
end

"Return the outer resolved boundary of `shape`."
function boundary end

"Return cross-sectional area \\[m²\\]."
function area end

"Return the cross-sectional centroid as `(x, y)` \\[m\\]."
function centroid end

"Return the directional support coordinate at angle `φ` \\[m\\]."
function support end

"Return the inner radius of radial geometry \\[m\\]."
function r_in end

"Return the outer radius of radial geometry \\[m\\]."
function r_ex end

"Return the radial thickness of radial geometry \\[m\\]."
function thickness end

boundary(shape::DiskShape) = shape
boundary(shape::RectangleShape) = shape
boundary(shape::EllipseShape) = shape
boundary(shape::SectorShape) = shape
boundary(shape::PolygonShape) = shape
boundary(shape::AnnulusShape) = DiskShape(shape.ro)

area(shape::DiskShape) = π * shape.r^2
area(shape::RectangleShape) = shape.w * shape.h
area(shape::EllipseShape) = π * shape.a * shape.b
area(shape::SectorShape) = shape.span * (shape.ro^2 - shape.ri^2) / 2
area(shape::AnnulusShape) = π * (shape.ro^2 - shape.ri^2)
function area(shape::PolygonShape)
    twice = sum(eachindex(shape.points)) do index
        next = mod1(index + 1, length(shape.points))
        shape.points[index][1] * shape.points[next][2] -
            shape.points[next][1] * shape.points[index][2]
    end
    return abs(twice) / 2
end

centroid(shape::Union{DiskShape, RectangleShape, EllipseShape, AnnulusShape}) =
    (zero(eltype(shape)), zero(eltype(shape)))
function centroid(shape::SectorShape)
    full = isapprox(shape.span, oftype(shape.span, 2π))
    full && return (zero(shape.span), zero(shape.span))
    radial = 4 * sin(shape.span / 2) * (shape.ro^3 - shape.ri^3) /
             (3 * shape.span * (shape.ro^2 - shape.ri^2))
    angle = shape.φ0 + shape.span / 2
    return (radial * cos(angle), radial * sin(angle))
end
function centroid(shape::PolygonShape)
    signed_twice_area = zero(eltype(shape))
    xmoment = zero(eltype(shape))
    ymoment = zero(eltype(shape))
    for index in eachindex(shape.points)
        next = mod1(index + 1, length(shape.points))
        left = shape.points[index]
        right = shape.points[next]
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

support(shape::DiskShape, φ::Real) = shape.r
support(shape::AnnulusShape, φ::Real) = shape.ro
support(shape::RectangleShape, φ::Real) =
    abs(cos(φ)) * shape.w / 2 + abs(sin(φ)) * shape.h / 2
support(shape::EllipseShape, φ::Real) =
    hypot(shape.a * cos(φ), shape.b * sin(φ))
function support(shape::SectorShape, φ::Real)
    candidates = (shape.φ0, shape.φ0 + shape.span, φ, φ + π)
    values = map(candidates) do angle
        relative = mod(angle - shape.φ0, oftype(shape.span, 2π))
        relative <= shape.span || return -oftype(shape.span, Inf)
        direction = cos(angle - φ)
        radius = direction >= zero(direction) ? shape.ro : shape.ri
        return radius * direction
    end
    return maximum(values)
end
support(shape::PolygonShape, φ::Real) = maximum(
    point[1] * cos(φ) + point[2] * sin(φ) for point in shape.points
)

support(shape::DiskShape) = shape.r
support(shape::AnnulusShape) = shape.ro
support(shape::RectangleShape) = hypot(shape.w, shape.h) / 2
support(shape::EllipseShape) = max(shape.a, shape.b)
support(shape::SectorShape) = shape.ro
support(shape::PolygonShape) = maximum(hypot(point...) for point in shape.points)

r_in(shape::DiskShape) = zero(shape.r)
r_in(shape::AnnulusShape) = shape.ri
r_in(shape::SectorShape) = shape.ri
r_ex(shape::DiskShape) = shape.r
r_ex(shape::AnnulusShape) = shape.ro
r_ex(shape::SectorShape) = shape.ro
thickness(shape::DiskShape) = shape.r
thickness(shape::AnnulusShape) = shape.ro - shape.ri
thickness(shape::SectorShape) = shape.ro - shape.ri
