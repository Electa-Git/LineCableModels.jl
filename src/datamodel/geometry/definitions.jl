"""
$(TYPEDEF)

Supertype for intrinsic, unresolved cross-sectional geometry.
"""
abstract type AbstractPrimitiveDefinition{T <: Real} end

Base.eltype(::AbstractPrimitiveDefinition{T}) where {T} = T
Base.eltype(::Type{<:AbstractPrimitiveDefinition{T}}) where {T} = T

"""
$(TYPEDEF)

Define an intrinsic solid circle by its radius.

$(TYPEDFIELDS)
"""
struct DiskDefinition{T <: Real} <: AbstractPrimitiveDefinition{T}
    "Radius \\[m\\]."
    r::T

    function DiskDefinition{T}(r::T) where {T <: Real}
        isfinite(r) && r > zero(r) ||
            throw(DomainError(r, "disk radius must be positive and finite"))
        return new{T}(r)
    end
end

DiskDefinition(r::Real) = DiskDefinition{typeof(float(r))}(float(r))

"""
$(TYPEDEF)

Define an intrinsic rectangle centered on its local origin.

$(TYPEDFIELDS)
"""
struct RectangleDefinition{T <: Real} <: AbstractPrimitiveDefinition{T}
    "Width along the local x-axis \\[m\\]."
    w::T
    "Height along the local y-axis \\[m\\]."
    h::T

    function RectangleDefinition{T}(w::T, h::T) where {T <: Real}
        isfinite(w) && w > zero(w) ||
            throw(DomainError(w, "rectangle width must be positive and finite"))
        isfinite(h) && h > zero(h) ||
            throw(DomainError(h, "rectangle height must be positive and finite"))
        return new{T}(w, h)
    end
end

function RectangleDefinition(w::Real, h::Real)
    values = map(float, promote(w, h))
    return RectangleDefinition{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Define an intrinsic ellipse centered on its local origin.

$(TYPEDFIELDS)
"""
struct EllipseDefinition{T <: Real} <: AbstractPrimitiveDefinition{T}
    "Semi-axis along the local x-axis \\[m\\]."
    a::T
    "Semi-axis along the local y-axis \\[m\\]."
    b::T

    function EllipseDefinition{T}(a::T, b::T) where {T <: Real}
        isfinite(a) && a > zero(a) ||
            throw(DomainError(a, "ellipse semi-axis a must be positive and finite"))
        isfinite(b) && b > zero(b) ||
            throw(DomainError(b, "ellipse semi-axis b must be positive and finite"))
        return new{T}(a, b)
    end
end

function EllipseDefinition(a::Real, b::Real)
    values = map(float, promote(a, b))
    return EllipseDefinition{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Define an intrinsic circular or annular sector in local polar coordinates.

$(TYPEDFIELDS)
"""
struct SectorDefinition{T <: Real} <: AbstractPrimitiveDefinition{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Starting angle \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T

    function SectorDefinition{T}(ri::T, ro::T, φ0::T, span::T) where {T <: Real}
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

function SectorDefinition(ri::Real, ro::Real, φ0::Real, span::Real)
    values = map(float, promote(ri, ro, φ0, span))
    return SectorDefinition{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Define an intrinsic circular annulus.

$(TYPEDFIELDS)
"""
struct AnnulusDefinition{T <: Real} <: AbstractPrimitiveDefinition{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T

    function AnnulusDefinition{T}(ri::T, ro::T) where {T <: Real}
        isfinite(ri) && ri >= zero(ri) ||
            throw(DomainError(ri, "annulus inner radius must be nonnegative and finite"))
        isfinite(ro) && ro > ri ||
            throw(DomainError(ro, "annulus outer radius must exceed its inner radius"))
        return new{T}(ri, ro)
    end
end


function AnnulusDefinition(ri::Real, ro::Real)
    values = map(float, promote(ri, ro))
    return AnnulusDefinition{typeof(first(values))}(values...)
end

"""
$(TYPEDEF)

Define an outward conformal layer relative to a resolved boundary.

$(TYPEDFIELDS)
"""
struct ShellDefinition{T <: Real} <: AbstractPrimitiveDefinition{T}
    "Normal layer thickness \\[m\\]."
    t::T

    function ShellDefinition{T}(t::T) where {T <: Real}
        isfinite(t) && t > zero(t) ||
            throw(DomainError(t, "shell thickness must be positive and finite"))
        return new{T}(t)
    end
end

ShellDefinition(t::Real) = ShellDefinition{typeof(float(t))}(float(t))

"""
$(TYPEDEF)

Define an explicit simple polygon in local coordinates.

$(TYPEDFIELDS)
"""
struct PolygonDefinition{T <: Real, P <: Tuple} <: AbstractPrimitiveDefinition{T}
    "Ordered vertices as `(x, y)` pairs \\[m\\]."
    points::P

    function PolygonDefinition{T, P}(points::P) where {T <: Real, P <: Tuple}
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

function PolygonDefinition(points)
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
    return PolygonDefinition{T, typeof(converted)}(converted)
end

Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::DiskDefinition) where {T <: Real} =
    DiskDefinition{T}(convert(T, value.r))
Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::RectangleDefinition) where {T <: Real} =
    RectangleDefinition{T}(convert(T, value.w), convert(T, value.h))
Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::EllipseDefinition) where {T <: Real} =
    EllipseDefinition{T}(convert(T, value.a), convert(T, value.b))
Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::SectorDefinition) where {T <: Real} =
    SectorDefinition{T}(
        convert(T, value.ri), convert(T, value.ro),
        convert(T, value.φ0), convert(T, value.span)
    )
Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::AnnulusDefinition) where {T <: Real} =
    AnnulusDefinition{T}(convert(T, value.ri), convert(T, value.ro))
Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::ShellDefinition) where {T <: Real} =
    ShellDefinition{T}(convert(T, value.t))

function Base.convert(::Type{<:AbstractPrimitiveDefinition{T}}, value::PolygonDefinition) where {T <: Real}
    points = map(point -> (convert(T, point[1]), convert(T, point[2])), value.points)
    return PolygonDefinition{T, typeof(points)}(points)
end
