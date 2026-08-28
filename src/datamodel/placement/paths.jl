"Store a helical lay-length to mean-diameter ratio."
struct LayRatio{T <: Real}
    q::T
    function LayRatio{T}(q::T) where {T <: Real}
        isfinite(q) && q > zero(q) ||
            throw(DomainError(q, "lay ratio must be positive and finite"))
        return new{T}(q)
    end
end
LayRatio(q::Real) = LayRatio{typeof(float(q))}(float(q))

"Store an authoritative helical pitch length \\[m\\]."
struct Pitch{T <: Real}
    p::T
    function Pitch{T}(p::T) where {T <: Real}
        isfinite(p) && p > zero(p) ||
            throw(DomainError(p, "pitch must be positive and finite"))
        return new{T}(p)
    end
end
Pitch(p::Real) = Pitch{typeof(float(p))}(float(p))

"Store an authoritative helical lay angle relative to the cable axis \\[rad\\]."
struct LayAngle{T <: Real}
    α::T
    function LayAngle{T}(α::T) where {T <: Real}
        isfinite(α) && zero(α) < α < oftype(α, π / 2) ||
            throw(DomainError(α, "lay angle must lie in (0, π/2)"))
        return new{T}(α)
    end
end
LayAngle(α::Real) = LayAngle{typeof(float(α))}(float(α))

"""
$(TYPEDEF)

Represent one helical path from one authoritative lay definition.

$(TYPEDFIELDS)
"""
struct Helix{L, T <: Real}
    "Lay definition."
    lay::L
    "Handedness: `1` or `-1` \\[dimensionless\\]."
    dir::Int
    "Initial angular position \\[rad\\]."
    φ0::T

    function Helix(lay::L; dir::Integer = 1, φ0::Real = 0) where {L}
        lay isa Union{LayRatio, Pitch, LayAngle} ||
            throw(ArgumentError("helix lay must be LayRatio, Pitch, or LayAngle"))
        dir in (-1, 1) || throw(ArgumentError("helix direction must be -1 or 1"))
        angle = float(φ0)
        isfinite(angle) || throw(DomainError(angle, "helix initial angle must be finite"))
        return new{L, typeof(angle)}(lay, Int(dir), angle)
    end
end

"Return helical pitch at local radius `radius` \\[m\\]."
function pitch end
pitch(path::Helix{<:LayRatio}, radius::Real) = path.lay.q * (2 * radius)
pitch(path::Helix{<:Pitch}, radius::Real) = path.lay.p
pitch(path::Helix{<:LayAngle}, radius::Real) =
    2π * radius / tan(path.lay.α)

"Return helical lay angle at local radius `radius` \\[rad\\]."
function angle(path::Helix, radius::Real)
    radius >= zero(radius) || throw(DomainError(radius, "helix radius must be nonnegative"))
    return atan(2π * radius, pitch(path, radius))
end

"Return the helical path-length ratio at local radius `radius`."
function overlength(path::Helix, radius::Real)
    value = angle(path, radius)
    return inv(cos(value))
end

function overlength(path::Helix{<:LayRatio}, radius::Real)
    radius >= zero(radius) || throw(DomainError(
        radius,
        "helix radius must be nonnegative"
    ))
    ratio = path.lay.q
    return sqrt(one(ratio) + ((one(ratio) * pi) / ratio)^2)
end
