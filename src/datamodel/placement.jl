"""
$(TYPEDEF)

Place `n` members evenly on a circle about the center of the enclosing part.

$(TYPEDFIELDS)
"""
struct CircularPattern{T <: Real}
    "Number of members \\[dimensionless\\]."
    n::Int
    "Radius of the circle through the member centers \\[m\\]."
    r::T
    "Angle of the first member \\[rad\\]."
    φ0::T
    "Counterclockwise angular span of the members \\[rad\\]."
    span::T
    "Clearance between adjacent members, as a fraction of their width \\[dimensionless\\]."
    gap_frac::T
    function CircularPattern{T}(n::Int, r::T, φ0::T, span::T, gap_frac::T) where {T <: Real}
        return validate(new{T}(n, r, φ0, span, gap_frac))
    end
end

"""
$(TYPEDSIGNATURES)

Construct the circular pattern of `n` members on the circle of radius `r` \\[m\\]. The keywords
are the fields of [`CircularPattern`](@ref): by default, the first member lies at angle zero
and the members fill the full circle without clearance.

# Errors

- Throws `DomainError` when `n` is not positive, when `r` is negative or not finite, when
  `span` is outside (0, 2π\\], or when `φ0` or `gap_frac` is not finite or `gap_frac` is
  negative.
"""
function CircularPattern(
        n::Integer; r::Real, φ0::Real = 0, span::Real = 2π, gap_frac::Real = 0)
    values = map(float, promote(r, φ0, span, gap_frac))
    return CircularPattern{eltype(values)}(Int(n), values...)
end

function validate(pattern::CircularPattern)
    pattern.n > 0 || throw(DomainError(pattern.n, "CircularPattern.n must be positive"))
    isfinite(pattern.r) && pattern.r >= zero(pattern.r) ||
        throw(DomainError(pattern.r, "CircularPattern.r must be nonnegative and finite"))
    isfinite(pattern.φ0) ||
        throw(DomainError(pattern.φ0, "CircularPattern.φ0 must be finite"))
    zero(pattern.span) < pattern.span <= 2π ||
        throw(DomainError(pattern.span, "CircularPattern.span must lie in (0, 2π]"))
    isfinite(pattern.gap_frac) && pattern.gap_frac >= zero(pattern.gap_frac) ||
        throw(DomainError(
            pattern.gap_frac, "CircularPattern.gap_frac must be nonnegative and finite"))
    return pattern
end

# Round wires on a circle keep the requested clearance along the chord between neighbors.
function validate(pattern::CircularPattern, circle::Circle)
    pattern.n == 1 && return pattern
    intervals = isapprox(pattern.span, 2π) ? pattern.n : pattern.n - 1
    chord = 2pattern.r * sin(pattern.span / intervals / 2)
    required = 2circle.r * (1 + pattern.gap_frac)
    tolerance = sqrt(eps(float(chord))) * max(one(chord), chord)
    chord + tolerance >= required || throw(DomainError(pattern.n,
        "$(pattern.n) wires of radius $(circle.r) m do not fit on a circle of radius " *
        "$(pattern.r) m with clearance fraction $(pattern.gap_frac)"))
    return pattern
end

# Rectangular strands on a circle fit along its arc, with the requested clearance.
function validate(pattern::CircularPattern, rectangle::Rectangle)
    occupied = pattern.n * rectangle.w * (1 + pattern.gap_frac)
    available = pattern.r * pattern.span
    tolerance = sqrt(eps(float(available))) * max(one(available), available)
    occupied <= available + tolerance || throw(DomainError(pattern.n,
        "$(pattern.n) strands of width $(rectangle.w) m do not fit along a circle of " *
        "radius $(pattern.r) m and span $(pattern.span) rad"))
    return pattern
end

# A circular pattern with an explicit radius keeps it in any context.
resolve(::Union{EmptyBoundary, AbstractShape}, pattern::CircularPattern) = pattern

"""
    placements(pattern, member)

Return the course that `pattern` lays with `member`, as `(placements, shape)`: the placement of
each member in the frame of the enclosing part, and the cross-section that every member
shares. A circular pattern places a circle as it is. It bends a rectangle into a
[`BentStrip`](@ref) of the same area that spans the height of the rectangle about the pattern
radius.

# Errors

- Throws `DomainError` when the members do not fit on the pattern with its clearance.
"""
function placements end

function placements(pattern::CircularPattern, circle::Circle)
    validate(pattern, circle)
    step = pattern.n == 1 ? zero(pattern.span) :
           isapprox(pattern.span, 2π) ? pattern.span / pattern.n :
           pattern.span / (pattern.n - 1)
    members = [Placement(pattern.r * cos(pattern.φ0 + k * step),
                   pattern.r * sin(pattern.φ0 + k * step), pattern.φ0 + k * step)
               for k in 0:(pattern.n - 1)]
    return (placements = members, shape = circle)
end

function placements(pattern::CircularPattern, rectangle::Rectangle)
    validate(pattern, rectangle)
    inner = max(zero(pattern.r), pattern.r - rectangle.h / 2)
    outer = inner + rectangle.h
    span = 2area(rectangle) / (outer^2 - inner^2)
    step = pattern.n == 1 ? zero(pattern.span) : pattern.span / pattern.n
    members = [Placement(zero(span), zero(span), pattern.φ0 + k * step)
               for k in 0:(pattern.n - 1)]
    return (placements = members,
        shape = BentStrip(inner, outer, span, Placement(zero(span), zero(span))))
end

"""
$(TYPEDEF)

Store the ratio of the lay length of a helix to the mean diameter of its course.

$(TYPEDFIELDS)
"""
struct LayRatio{T <: Real}
    "Lay length over the mean lay diameter \\[dimensionless\\]."
    q::T
    function LayRatio{T}(q::T) where {T <: Real}
        return validate(new{T}(q))
    end
end

"""
$(TYPEDSIGNATURES)

Construct the lay ratio `q` \\[dimensionless\\].

# Errors

- Throws `DomainError` when `q` is not positive and finite.
"""
LayRatio(q::Real) = LayRatio{typeof(float(q))}(float(q))

function validate(lay::LayRatio)
    isfinite(lay.q) && lay.q > zero(lay.q) ||
        throw(DomainError(lay.q, "LayRatio.q must be positive and finite"))
    return lay
end

"""
$(TYPEDEF)

Represent the helical path of the members of a course along the cable.

$(TYPEDFIELDS)
"""
struct Helix{L <: LayRatio, T <: Real}
    "Lay of the helix."
    lay::L
    "Handedness of the lay, `1` or `-1` \\[dimensionless\\]."
    dir::Int
    "Angle at which the helix starts \\[rad\\]."
    φ0::T
    function Helix{L, T}(lay::L, dir::Int, φ0::T) where {L <: LayRatio, T <: Real}
        return validate(new{L, T}(lay, dir, φ0))
    end
end

"""
$(TYPEDSIGNATURES)

Construct the helix of `lay`. The keywords are the fields of [`Helix`](@ref): by default, the
handedness is `1` and the helix starts at angle zero.

# Errors

- Throws `ArgumentError` when `dir` is not `1` or `-1`, and `DomainError` when `φ0` is not
  finite.
"""
function Helix(lay::LayRatio; dir::Integer = 1, φ0::Real = 0)
    return Helix{typeof(lay), typeof(float(φ0))}(lay, Int(dir), float(φ0))
end

function validate(helix::Helix)
    helix.dir in (-1, 1) || throw(ArgumentError("Helix.dir must be 1 or -1"))
    isfinite(helix.φ0) || throw(DomainError(helix.φ0, "Helix.φ0 must be finite"))
    return helix
end

"""
$(TYPEDSIGNATURES)

Return the ratio of the length of `helix` to the length of the cable, at the lay radius
`radius` \\[m\\]:

```math
\\frac{\\ell}{\\ell_0} = \\sqrt{1 + \\left(\\frac{\\pi}{q}\\right)^2},
```

where ``q`` is the lay ratio. A path at zero radius lies straight along the axis, with ratio 1.
"""
function overlength(helix::Helix{<:LayRatio}, radius::Real)
    factor = sqrt(1 + (π / helix.lay.q)^2)
    return iszero(radius) ? one(factor) : factor
end

"""
$(TYPEDSIGNATURES)

Return the lay length of `helix` at the lay radius `radius` \\[m\\]: ``p = 2 q r``, where ``q``
is the lay ratio and ``r`` the radius.
"""
pitch(helix::Helix{<:LayRatio}, radius::Real) = 2helix.lay.q * radius
