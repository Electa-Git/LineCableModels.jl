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

Construct a placement in the type that its coordinates promote to with `Float64`. The
arguments are the fields of [`Placement`](@ref), and the rotation is zero unless given.

# Errors

- Throws `DomainError` when a coordinate or the rotation is not finite.
"""
function Placement(x::Real, y::Real, φ::Real = 0)
    values = promote(x, y, φ)
    T = promote_type(Float64, eltype(values))
    return Placement{T}(convert.(T, values)...)
end

function validate(placement::Placement)
    all(isfinite, (placement.x, placement.y, placement.φ)) || throw(DomainError(
        (placement.x, placement.y, placement.φ),
        "Placement coordinates and rotation must be finite"))
    return placement
end

"""
$(TYPEDSIGNATURES)

Return, in \\[m\\], the distance of `point` from the origin:

```math
d = \\sqrt{x^2 + y^2},
```

where ``x`` and ``y`` are the coordinates of `point`, in \\[m\\]. When `iszero` is true for
both coordinates, the result is the exact zero `zero(x)`. A number type with derivatives,
such as a dual number, then has zero derivatives there, where those of `hypot` are
undefined.
"""
function distance(point::Tuple{Real, Real})
    iszero(point[1]) && iszero(point[2]) ? zero(point[1]) : hypot(point...)
end

"""
$(TYPEDSIGNATURES)

Return, in \\[m\\], the distance of the position of `at` from the origin of its parent
frame. It is the distance of the point ``(x, y)``, the horizontal and vertical translations
of `at`.
"""
distance(at::Placement) = distance((at.x, at.y))

Base.eltype(::Placement{T}) where {T} = T

function Base.map(f, placement::Placement)
    return Placement(f(placement.x), f(placement.y), f(placement.φ))
end

function Base.convert(::Type{Placement{T}}, placement::Placement) where {T <: Real}
    return Placement{T}(
        convert(T, placement.x), convert(T, placement.y), convert(T, placement.φ))
end

"""
$(TYPEDSIGNATURES)

Compose `child`, a placement in the frame of `parent`, into the frame that holds `parent`:

```math
x = x_p + x_c \\cos\\phi_p - y_c \\sin\\phi_p, \\qquad
y = y_p + x_c \\sin\\phi_p + y_c \\cos\\phi_p, \\qquad
\\phi = \\phi_p + \\phi_c,
```

where ``(x_p, y_p)`` and ``(x_c, y_c)`` are the positions of `parent` and `child` \\[m\\], and
``\\phi_p`` and ``\\phi_c`` are their rotations \\[rad\\].
"""
function Base.:*(parent::Placement, child::Placement)
    return Placement(parent.x + cos(parent.φ) * child.x - sin(parent.φ) * child.y,
        parent.y + sin(parent.φ) * child.x + cos(parent.φ) * child.y, parent.φ + child.φ)
end

"""
$(TYPEDSIGNATURES)

Return the placement that undoes `p`, so that `inv(p) * p` and `p * inv(p)` are the
identity:

```math
x = -(x_p \\cos\\phi_p + y_p \\sin\\phi_p), \\qquad y = x_p \\sin\\phi_p - y_p \\cos\\phi_p,
\\qquad \\phi = -\\phi_p,
```

where ``(x_p, y_p)`` is the position of `p` \\[m\\] and ``\\phi_p`` its rotation \\[rad\\].
"""
function Base.inv(p::Placement)
    return Placement(-(p.x * cos(p.φ) + p.y * sin(p.φ)),
        p.x * sin(p.φ) - p.y * cos(p.φ), -p.φ)
end

"""
$(TYPEDSIGNATURES)

Rotate `point` by the rotation of `placement`, then translate it by the position of
`placement`:

```math
x = x_0 + u \\cos\\phi - v \\sin\\phi, \\qquad y = y_0 + u \\sin\\phi + v \\cos\\phi,
```

where ``(u, v)`` is `point` \\[m\\], ``(x_0, y_0)`` is the position of `placement` \\[m\\] and
``\\phi`` is its rotation \\[rad\\].

# Returns

- The point `(x, y)` in the parent frame \\[m\\].
"""
function (placement::Placement)(point::Tuple{Real, Real})
    return (placement.x + cos(placement.φ) * point[1] - sin(placement.φ) * point[2],
        placement.y + sin(placement.φ) * point[1] + cos(placement.φ) * point[2])
end
