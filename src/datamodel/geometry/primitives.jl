"""
$(TYPEDEF)

Supertype for intrinsic, unresolved cross-sectional geometry.
"""
abstract type AbstractPrimitive{T <: Real} end

Base.eltype(::AbstractPrimitive{T}) where {T} = T
Base.eltype(::Type{<:AbstractPrimitive{T}}) where {T} = T

"""
$(TYPEDEF)

Represent a solid circular primitive.

$(TYPEDFIELDS)
"""
struct Disk{T <: Real} <: AbstractPrimitive{T}
    "Radius \\[m\\]."
    r::T

    function Disk{T}(r::T) where {T <: Real}
        isfinite(r) && r > zero(r) ||
            throw(DomainError(r, "disk radius must be positive and finite"))
        return new{T}(r)
    end
end

Disk(r::Real) = Disk{typeof(float(r))}(float(r))

"""
$(TYPEDEF)

Represent a centered rectangular primitive.

$(TYPEDFIELDS)
"""
struct Rectangle{T <: Real} <: AbstractPrimitive{T}
    "Width along the local x-axis \\[m\\]."
    w::T
    "Height along the local y-axis \\[m\\]."
    h::T

    function Rectangle{T}(w::T, h::T) where {T <: Real}
        isfinite(w) && w > zero(w) ||
            throw(DomainError(w, "rectangle width must be positive and finite"))
        isfinite(h) && h > zero(h) ||
            throw(DomainError(h, "rectangle height must be positive and finite"))
        return new{T}(w, h)
    end
end

function Rectangle(w::Real, h::Real)
    values = map(float, promote(w, h))
    return Rectangle{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Represent a centered elliptical primitive.

$(TYPEDFIELDS)
"""
struct Ellipse{T <: Real} <: AbstractPrimitive{T}
    "Semi-axis along the local x-axis \\[m\\]."
    a::T
    "Semi-axis along the local y-axis \\[m\\]."
    b::T

    function Ellipse{T}(a::T, b::T) where {T <: Real}
        isfinite(a) && a > zero(a) ||
            throw(DomainError(a, "ellipse semi-axis a must be positive and finite"))
        isfinite(b) && b > zero(b) ||
            throw(DomainError(b, "ellipse semi-axis b must be positive and finite"))
        return new{T}(a, b)
    end
end

function Ellipse(a::Real, b::Real)
    values = map(float, promote(a, b))
    return Ellipse{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Represent a circular sector or annular sector in local polar coordinates.

$(TYPEDFIELDS)
"""
struct Sector{T <: Real} <: AbstractPrimitive{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Starting angle \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T

    function Sector{T}(ri::T, ro::T, φ0::T, span::T) where {T <: Real}
        isfinite(ri) && ri >= zero(ri) ||
            throw(DomainError(ri, "sector inner radius must be nonnegative and finite"))
        isfinite(ro) && ro > ri ||
            throw(DomainError(ro, "sector outer radius must exceed its inner radius"))
        isfinite(φ0) || throw(DomainError(φ0, "sector start angle must be finite"))
        isfinite(span) && zero(span) < span <= oftype(span, 2π) ||
            throw(DomainError(span, "sector span must lie in (0, 2π]"))
        return new{T}(ri, ro, φ0, span)
    end
end

function Sector(ri::Real, ro::Real, φ0::Real, span::Real)
    values = map(float, promote(ri, ro, φ0, span))
    return Sector{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Represent an absolute circular annulus.

$(TYPEDFIELDS)
"""
struct Annulus{T <: Real} <: AbstractPrimitive{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T

    function Annulus{T}(ri::T, ro::T) where {T <: Real}
        isfinite(ri) && ri >= zero(ri) ||
            throw(DomainError(ri, "annulus inner radius must be nonnegative and finite"))
        isfinite(ro) && ro > ri ||
            throw(DomainError(ro, "annulus outer radius must exceed its inner radius"))
        return new{T}(ri, ro)
    end
end


function Annulus(ri::Real, ro::Real)
    values = map(float, promote(ri, ro))
    return Annulus{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Represent an outward conformal layer relative to a resolved boundary.

$(TYPEDFIELDS)
"""
struct Shell{T <: Real} <: AbstractPrimitive{T}
    "Normal layer thickness \\[m\\]."
    t::T

    function Shell{T}(t::T) where {T <: Real}
        isfinite(t) && t > zero(t) ||
            throw(DomainError(t, "shell thickness must be positive and finite"))
        return new{T}(t)
    end
end

Shell(t::Real) = Shell{typeof(float(t))}(float(t))

"""
$(TYPEDEF)

Represent an explicit simple polygon in local coordinates.

$(TYPEDFIELDS)
"""
struct Polygon{T <: Real, P <: Tuple} <: AbstractPrimitive{T}
    "Ordered vertices as `(x, y)` pairs \\[m\\]."
    points::P

    function Polygon{T, P}(points::P) where {T <: Real, P <: Tuple}
        length(points) >= 3 ||
            throw(ArgumentError("a polygon requires at least three vertices"))
        all(point -> length(point) == 2, points) ||
            throw(ArgumentError("polygon vertices must contain two coordinates"))
        all(point -> all(isfinite, point), points) ||
            throw(ArgumentError("polygon coordinates must be finite"))
        signed_twice_area = sum(eachindex(points)) do index
            next = mod1(index + 1, length(points))
            points[index][1] * points[next][2] -
                points[next][1] * points[index][2]
        end
        iszero(signed_twice_area) &&
            throw(DomainError(signed_twice_area, "polygon area must be nonzero"))
        return new{T, P}(points)
    end
end

function Polygon(points)
    vertices = Tuple(points)
    length(vertices) >= 3 ||
        throw(ArgumentError("a polygon requires at least three vertices"))
    all(point -> applicable(length, point) && length(point) == 2, vertices) ||
        throw(ArgumentError("polygon vertices must contain two coordinates"))
    all(point -> point[1] isa Real && point[2] isa Real, vertices) ||
        throw(ArgumentError("polygon coordinates must be real numbers"))
    T = promote_type(map(point -> typeof(float(point[1])), vertices)...,
        map(point -> typeof(float(point[2])), vertices)...)
    converted = map(point -> (convert(T, point[1]), convert(T, point[2])), vertices)
    return Polygon{T, typeof(converted)}(converted)
end

Base.convert(::Type{<:AbstractPrimitive{T}}, value::Disk) where {T <: Real} =
    Disk{T}(convert(T, value.r))
Base.convert(::Type{<:AbstractPrimitive{T}}, value::Rectangle) where {T <: Real} =
    Rectangle{T}(convert(T, value.w), convert(T, value.h))
Base.convert(::Type{<:AbstractPrimitive{T}}, value::Ellipse) where {T <: Real} =
    Ellipse{T}(convert(T, value.a), convert(T, value.b))
Base.convert(::Type{<:AbstractPrimitive{T}}, value::Sector) where {T <: Real} =
    Sector{T}(
        convert(T, value.ri), convert(T, value.ro),
        convert(T, value.φ0), convert(T, value.span)
    )
Base.convert(::Type{<:AbstractPrimitive{T}}, value::Annulus) where {T <: Real} =
    Annulus{T}(convert(T, value.ri), convert(T, value.ro))
Base.convert(::Type{<:AbstractPrimitive{T}}, value::Shell) where {T <: Real} =
    Shell{T}(convert(T, value.t))

function Base.convert(::Type{<:AbstractPrimitive{T}}, value::Polygon) where {T <: Real}
    points = map(point -> (convert(T, point[1]), convert(T, point[2])), value.points)
    return Polygon{T, typeof(points)}(points)
end
