"""
$(TYPEDEF)

Represent a planar translation and a rotation relative to a parent frame.

$(TYPEDFIELDS)
"""
struct Placement{T <: Real}
    "Horizontal translation \\[m\\]."
    x::T
    "Vertical translation \\[m\\]."
    y::T
    "Counterclockwise rotation \\[rad\\]."
    φ::T
    function Placement{T}(x::T, y::T, φ::T) where {T <: Real}
        return validate(new{T}(x, y, φ))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a placement after floating and promoting its coordinates. The arguments are the
fields of [`Placement`](@ref), and the rotation is zero unless given.

# Errors

- Throws `DomainError` when a coordinate or the rotation is not finite.
"""
function Placement(x::Real, y::Real, φ::Real = 0)
    values = map(float, promote(x, y, φ))
    return Placement{eltype(values)}(values...)
end

function validate(placement::Placement)
    all(isfinite, (placement.x, placement.y, placement.φ)) || throw(DomainError(
        (placement.x, placement.y, placement.φ),
        "Placement coordinates and rotation must be finite"))
    return placement
end

Base.eltype(::Placement{T}) where {T} = T

function Base.convert(::Type{Placement{T}}, placement::Placement) where {T <: Real}
    return Placement{T}(
        convert(T, placement.x), convert(T, placement.y), convert(T, placement.φ))
end

"""
$(TYPEDSIGNATURES)

Compose `child`, a placement in the frame of `parent`, into the frame that holds `parent`.
"""
function Base.:*(parent::Placement, child::Placement)
    return Placement(parent.x + cos(parent.φ) * child.x - sin(parent.φ) * child.y,
        parent.y + sin(parent.φ) * child.x + cos(parent.φ) * child.y, parent.φ + child.φ)
end

"""
$(TYPEDSIGNATURES)

Rotate `point` by the rotation of `placement`, then translate it by the position of
`placement`.

# Returns

- The point `(x, y)` in the parent frame \\[m\\].
"""
function (placement::Placement)(point::Tuple{Real, Real})
    return (placement.x + cos(placement.φ) * point[1] - sin(placement.φ) * point[2],
        placement.y + sin(placement.φ) * point[1] + cos(placement.φ) * point[2])
end

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

Construct a circle of radius `r` \\[m\\] at the placement `at`, in the floating type that both
promote to.

# Errors

- Throws `DomainError` when `r` is not positive and finite.
"""
function Circle(r::Real, at::Placement)
    T = float(promote_type(typeof(r), eltype(at)))
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
floating type that they promote to.

# Errors

- Throws `DomainError` when `ri` is negative or not finite, or when `ro` is not finite or does
  not exceed `ri`.
"""
function Annulus(ri::Real, ro::Real, at::Placement)
    T = float(promote_type(typeof(ri), typeof(ro), eltype(at)))
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

Represent a rectangular strand before a circular pattern bends it onto its course: its width
lies along the course and its height across it.

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

Construct a rectangle of width `w` and height `h` \\[m\\] centered at the origin, in the floating
type that they promote to.

# Errors

- Throws `DomainError` when `w` or `h` is not positive and finite.
"""
function Rectangle(w::Real, h::Real)
    T = float(promote_type(typeof(w), typeof(h)))
    return Rectangle{T}(convert(T, w), convert(T, h), Placement(zero(T), zero(T)))
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

Declare a layer of uniform thickness around the part that it follows, such as an insulation,
a semiconducting screen or a jacket.

$(TYPEDFIELDS)
"""
struct Shell{T <: Real}
    "Thickness \\[m\\]."
    t::T
    function Shell{T}(t::T) where {T <: Real}
        return validate(new{T}(t))
    end
end

"""
$(TYPEDSIGNATURES)

Declare a layer of thickness `t` \\[m\\].

# Errors

- Throws `DomainError` when `t` is not positive and finite.
"""
Shell(t::Real) = Shell{typeof(float(t))}(float(t))

function validate(shell::Shell)
    isfinite(shell.t) && shell.t > zero(shell.t) ||
        throw(DomainError(shell.t, "Shell.t must be positive and finite"))
    return shell
end

"""
$(TYPEDEF)

Represent a rectangular strand bent onto its course. The strip is an annular sector,
symmetric about the x-axis of its placement. Its area equals the area of the rectangle of width
``w`` and height ``h``:

```math
\\frac{\\Delta\\phi}{2}\\left(r_o^2 - r_i^2\\right) = w h,
```

where ``\\Delta\\phi`` is the angular span and ``r_i`` and ``r_o`` are the inner and outer
radii.

$(TYPEDFIELDS)
"""
struct BentStrip{T <: Real} <: AbstractShape{T}
    "Inner radius \\[m\\]."
    ri::T
    "Outer radius \\[m\\]."
    ro::T
    "Angular span \\[rad\\]."
    span::T
    "Placement of the center of curvature. Its x-axis is the symmetry axis of the strip."
    at::Placement{T}
    function BentStrip{T}(ri::T, ro::T, span::T, at::Placement{T}) where {T <: Real}
        return validate(new{T}(ri, ro, span, at))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a bent strip from its fields, in the floating type that they promote to.

# Errors

- Throws `DomainError` when a radius is negative or not finite, when `ro` does not exceed
  `ri`, or when `span` is outside (0, 2π\\].
"""
function BentStrip(ri::Real, ro::Real, span::Real, at::Placement)
    T = float(promote_type(typeof(ri), typeof(ro), typeof(span), eltype(at)))
    return BentStrip{T}(convert(T, ri), convert(T, ro), convert(T, span),
        convert(Placement{T}, at))
end

function validate(strip::BentStrip)
    isfinite(strip.ri) && strip.ri >= zero(strip.ri) ||
        throw(DomainError(strip.ri, "BentStrip.ri must be nonnegative and finite"))
    isfinite(strip.ro) && strip.ro > strip.ri ||
        throw(DomainError(strip.ro, "BentStrip.ro must be finite and exceed BentStrip.ri"))
    zero(strip.span) < strip.span <= 2π ||
        throw(DomainError(strip.span, "BentStrip.span must lie in (0, 2π]"))
    return strip
end

"""
    area(shape)

Return the cross-sectional area of `shape` \\[m²\\].
"""
function area end

area(circle::Circle) = π * circle.r^2
area(annulus::Annulus) = π * (annulus.ro^2 - annulus.ri^2)
area(rectangle::Rectangle) = rectangle.w * rectangle.h
area(strip::BentStrip) = strip.span * (strip.ro^2 - strip.ri^2) / 2

"""
    centroid(shape)

Return the centroid `(x, y)` of `shape` in the cable cross-section \\[m\\]. For a bent strip, the
centroid lies on its symmetry axis at the distance

```math
d = \\frac{4 \\sin(\\Delta\\phi / 2) \\left(r_o^3 - r_i^3\\right)}
    {3 \\Delta\\phi \\left(r_o^2 - r_i^2\\right)}
```

from its center of curvature, where ``\\Delta\\phi`` is its angular span and ``r_i`` and ``r_o``
are its inner and outer radii.
"""
function centroid end

centroid(shape::Union{Circle, Annulus}) = (shape.at.x, shape.at.y)

function centroid(strip::BentStrip)
    distance = 4sin(strip.span / 2) * (strip.ro^3 - strip.ri^3) /
               (3strip.span * (strip.ro^2 - strip.ri^2))
    return strip.at((distance, zero(distance)))
end

"""
    support(shape)

Return the distance from the origin to the center of `shape` plus its outer radius \\[m\\]: the
radius of a circle about the origin that contains `shape`.
"""
function support end

support(circle::Circle) = hypot(circle.at.x, circle.at.y) + circle.r
support(annulus::Annulus) = hypot(annulus.at.x, annulus.at.y) + annulus.ro
support(strip::BentStrip) = hypot(strip.at.x, strip.at.y) + strip.ro

"""
    boundary(shape)

Return the outline of `shape` for the next part. The outline of a circle is the circle itself.
The outline of an annulus is the circle of its outer radius.
"""
function boundary end

boundary(circle::Circle) = circle
boundary(annulus::Annulus) = Circle(annulus.ro, annulus.at)

"""
    r_in(shape)

Return the inner radius of a circular shape \\[m\\]: zero for a circle.
"""
function r_in end

r_in(circle::Circle) = zero(circle.r)
r_in(annulus::Annulus) = annulus.ri
r_in(strip::BentStrip) = strip.ri

"""
    r_ex(shape)

Return the outer radius of a circular shape \\[m\\].
"""
function r_ex end

r_ex(circle::Circle) = circle.r
r_ex(annulus::Annulus) = annulus.ro
r_ex(strip::BentStrip) = strip.ro

# A shell around a circle or an annulus is the annulus of its thickness at the same placement.
function resolve(inner::Union{Circle, Annulus}, shell::Shell)
    return Annulus(r_ex(inner), r_ex(inner) + shell.t, inner.at)
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
function resolve(at::Placement, strip::BentStrip)
    return BentStrip(strip.ri, strip.ro, strip.span, at * strip.at)
end
