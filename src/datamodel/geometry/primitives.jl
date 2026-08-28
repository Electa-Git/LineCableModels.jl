"""
$(TYPEDEF)

Supertype for exact cross-sectional geometry.

Intrinsic primitives and derived domains such as conformal shells share this
interface. A shape is independent of any plotting tessellation.
"""
abstract type AbstractShape{T <: Real} end

Base.eltype(::AbstractShape{T}) where {T} = T
Base.eltype(::Type{<:AbstractShape{T}}) where {T} = T

"""
$(TYPEDEF)

Supertype for intrinsic cross-sectional primitives.

A primitive states material-neutral dimensions. During construction, simple
primitives also carry the `Pose2` produced by coordinate composition. Shapes
that require additional exact contact geometry use a separate resolved
`AbstractShape` implementation.
"""
abstract type AbstractPrimitive{T <: Real} <: AbstractShape{T} end

Base.eltype(::AbstractPrimitive{T}) where {T} = T
Base.eltype(::Type{<:AbstractPrimitive{T}}) where {T} = T

"Explicit initial state for resolving the first member of an outward stack."
struct EmptyBoundary end

"""
$(TYPEDEF)

Represent a circular cross-section and its composed pose.

$(TYPEDFIELDS)
"""
struct Disk{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    "Radius \\[m\\]."
    r::T
    "Pose in the completed cable coordinate system."
    at::P

    function Disk{T, P}(r::T, at::P) where {T <: Real, P <: Pose2{T}}
        isfinite(r) && r > zero(r) ||
            throw(DomainError(r, "disk radius must be positive and finite"))
        return new{T, P}(r, at)
    end
end

"""
$(TYPEDEF)

Represent a rectangular cross-section and its composed pose.

$(TYPEDFIELDS)
"""
struct Rectangle{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    "Width along the local x-axis \\[m\\]."
    w::T
    "Height along the local y-axis \\[m\\]."
    h::T
    "Pose in the completed cable coordinate system."
    at::P

    function Rectangle{T, P}(w::T, h::T, at::P) where {T <: Real, P <: Pose2{T}}
        isfinite(w) && w > zero(w) ||
            throw(DomainError(w, "rectangle width must be positive and finite"))
        isfinite(h) && h > zero(h) ||
            throw(DomainError(h, "rectangle height must be positive and finite"))
        return new{T, P}(w, h, at)
    end
end

"""
$(TYPEDEF)

Represent an elliptical cross-section and its composed pose.

$(TYPEDFIELDS)
"""
struct Ellipse{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    "Semi-axis along the local x-axis \\[m\\]."
    a::T
    "Semi-axis along the local y-axis \\[m\\]."
    b::T
    "Pose in the completed cable coordinate system."
    at::P

    function Ellipse{T, P}(a::T, b::T, at::P) where {T <: Real, P <: Pose2{T}}
        isfinite(a) && a > zero(a) ||
            throw(DomainError(a, "ellipse semi-axis a must be positive and finite"))
        isfinite(b) && b > zero(b) ||
            throw(DomainError(b, "ellipse semi-axis b must be positive and finite"))
        return new{T, P}(a, b, at)
    end
end

"""
$(TYPEDEF)

Represent a circular or annular sector and its composed pose.

$(TYPEDFIELDS)
"""
struct Sector{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Starting angle in the local frame \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T
    "Pose in the completed cable coordinate system."
    at::P

    function Sector{T, P}(
            ri::T, ro::T, φ0::T, span::T, at::P
    ) where {T <: Real, P <: Pose2{T}}
        isfinite(ri) && ri >= zero(ri) || throw(DomainError(
            ri, "sector inner radius must be nonnegative and finite"
        ))
        isfinite(ro) && ro > ri || throw(DomainError(
            ro, "sector outer radius must exceed its inner radius"
        ))
        isfinite(φ0) || throw(DomainError(φ0, "sector start angle must be finite"))
        isfinite(span) && zero(span) < span <= oftype(span, 2π) || throw(DomainError(
            span, "sector span must lie in (0, 2π]"
        ))
        return new{T, P}(ri, ro, φ0, span, at)
    end
end

"""
$(TYPEDEF)

Represent an annular cross-section and its composed pose.

$(TYPEDFIELDS)
"""
struct Annulus{T <: Real, P <: Pose2{T}} <: AbstractPrimitive{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Pose in the completed cable coordinate system."
    at::P

    function Annulus{T, P}(ri::T, ro::T, at::P) where {T <: Real, P <: Pose2{T}}
        isfinite(ri) && ri >= zero(ri) || throw(DomainError(
            ri, "annulus inner radius must be nonnegative and finite"
        ))
        isfinite(ro) && ro > ri || throw(DomainError(
            ro, "annulus outer radius must exceed its inner radius"
        ))
        return new{T, P}(ri, ro, at)
    end
end

"""
$(TYPEDEF)

Represent a polygonal cross-section and its composed pose.

$(TYPEDFIELDS)
"""
struct Polygon{T <: Real, V <: Tuple, P <: Pose2{T}} <: AbstractPrimitive{T}
    "Ordered local `(x, y)` vertices \\[m\\]."
    points::V
    "Pose in the completed cable coordinate system."
    at::P

    function Polygon{T, V, P}(points::V, at::P) where {
            T <: Real, V <: Tuple, P <: Pose2{T}
    }
        length(points) >= 3 ||
            throw(ArgumentError("a polygon requires at least three vertices"))
        all(point -> length(point) == 2 && all(isfinite, point), points) ||
            throw(ArgumentError("polygon vertices must be finite coordinate pairs"))
        twice_area = sum(eachindex(points)) do index
            next = mod1(index + 1, length(points))
            points[index][1] * points[next][2] -
                points[next][1] * points[index][2]
        end
        iszero(twice_area) && throw(DomainError(
            twice_area, "polygon area must be nonzero"
        ))
        return new{T, V, P}(points, at)
    end
end

struct _DifferencePrimitive{
        T <: Real,
        O <: AbstractShape,
        H <: Tuple,
        P <: Pose2{T}
} <: AbstractShape{T}
    outer::O
    holes::H
    at::P

    function _DifferencePrimitive(
            outer::O,
            holes::H
    ) where {O <: AbstractShape, H <: Tuple}
        all(hole -> hole isa AbstractShape, holes) || throw(ArgumentError(
            "difference holes must be exact resolved shapes"
        ))
        T = promote_type(eltype(outer), map(eltype, holes)...)
        at = convert(Pose2{T}, outer.at)
        return new{T, O, H, typeof(at)}(outer, holes, at)
    end
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
function Rectangle(w::Real, h::Real)
    T = promote_type(typeof(w), typeof(h))
    return Rectangle(convert(T, w), convert(T, h), _origin(T))
end
function Ellipse(a::Real, b::Real)
    T = promote_type(typeof(a), typeof(b))
    return Ellipse(convert(T, a), convert(T, b), _origin(T))
end
function Sector(ri::Real, ro::Real, φ0::Real, span::Real)
    T = promote_type(typeof(ri), typeof(ro), typeof(φ0), typeof(span))
    return Sector(
        convert(T, ri), convert(T, ro), convert(T, φ0), convert(T, span), _origin(T)
    )
end
function Annulus(ri::Real, ro::Real)
    T = promote_type(typeof(ri), typeof(ro))
    return Annulus(convert(T, ri), convert(T, ro), _origin(T))
end
function Polygon(points::V) where {V <: Tuple}
    return _polygon(points, _origin(typeof(float(first(points)[1]))))
end

function Base.convert(
        ::Type{<:AbstractPrimitive{T}}, value::Disk
) where {T <: Real}
    return Disk(convert(T, value.r), convert(Pose2{T}, value.at))
end
function Base.convert(
        ::Type{<:AbstractPrimitive{T}}, value::Rectangle
) where {T <: Real}
    return Rectangle(
        convert(T, value.w), convert(T, value.h), convert(Pose2{T}, value.at)
    )
end
function Base.convert(
        ::Type{<:AbstractPrimitive{T}}, value::Ellipse
) where {T <: Real}
    return Ellipse(
        convert(T, value.a), convert(T, value.b), convert(Pose2{T}, value.at)
    )
end
function Base.convert(
        ::Type{<:AbstractPrimitive{T}}, value::Sector
) where {T <: Real}
    return Sector(
        convert(T, value.ri), convert(T, value.ro), convert(T, value.φ0),
        convert(T, value.span), convert(Pose2{T}, value.at)
    )
end
function Base.convert(
        ::Type{<:AbstractPrimitive{T}}, value::Annulus
) where {T <: Real}
    return Annulus(
        convert(T, value.ri), convert(T, value.ro), convert(Pose2{T}, value.at)
    )
end
function Base.convert(
        ::Type{<:AbstractPrimitive{T}}, value::Polygon
) where {T <: Real}
    points = map(point -> (convert(T, point[1]), convert(T, point[2])), value.points)
    return _polygon(points, convert(Pose2{T}, value.at))
end

"Return the outer resolved boundary of `primitive`."
function boundary end

"Return cross-sectional area \\[m²\\]."
function area end

"Return exact boundary perimeter in metres."
function perimeter end

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

perimeter(primitive::Disk) = 2π * primitive.r
perimeter(primitive::Rectangle) = 2 * (primitive.w + primitive.h)
perimeter(primitive::Sector) =
    primitive.span * (primitive.ro + primitive.ri) +
    2 * (primitive.ro - primitive.ri)
perimeter(primitive::Annulus) = 2π * (primitive.ro + primitive.ri)
function perimeter(primitive::Polygon)
    return sum(eachindex(primitive.points)) do index
        next = mod1(index + 1, length(primitive.points))
        hypot(
            primitive.points[next][1] - primitive.points[index][1],
            primitive.points[next][2] - primitive.points[index][2]
        )
    end
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

boundary(primitive::_DifferencePrimitive) = boundary(primitive.outer)
support(primitive::_DifferencePrimitive, φ::Real) = support(primitive.outer, φ)
support(primitive::_DifferencePrimitive) = support(primitive.outer)
area(primitive::_DifferencePrimitive) =
    area(primitive.outer) - sum(area, primitive.holes; init = zero(eltype(primitive)))
function centroid(primitive::_DifferencePrimitive)
    total_area = area(primitive)
    total_area > zero(total_area) || throw(DomainError(
        total_area, "difference primitive must have positive area"
    ))
    outer_area = area(primitive.outer)
    outer_centroid = centroid(primitive.outer)
    xmoment = outer_area * outer_centroid[1]
    ymoment = outer_area * outer_centroid[2]
    for hole in primitive.holes
        hole_area = area(hole)
        hole_centroid = centroid(hole)
        xmoment -= hole_area * hole_centroid[1]
        ymoment -= hole_area * hole_centroid[2]
    end
    return (xmoment / total_area, ymoment / total_area)
end
