"""
Store a helical lay-length to mean-diameter ratio.
"""
struct LayRatio{T <: Real}
    q::T
    LayRatio{T}(q::T) where {T <: Real} = validate(new{T}(q))
end

function validate(lay::LayRatio)
    isfinite(lay.q) && lay.q > zero(lay.q) ||
        throw(DomainError(lay.q, "lay ratio must be positive and finite"))
    return lay
end

LayRatio(q::Real) = LayRatio{typeof(float(q))}(float(q))

"""
$(TYPEDSIGNATURES)

Declare one lay-length ratio or a homogeneous schedule of lay-length ratios.

# Arguments

- `q`: one positive lay-length to mean-diameter ratio \\[dimensionless\\].
- `values`: two or more ratios, or one tuple or vector of ratios. Specify one ratio per repeated course.

# Keywords

- `combine=:product`: gridspace composition rule.

# Returns

- A `LayRatio`, a tuple of `LayRatio` values, or the corresponding `Gridspace`
  when an explicit finite source is supplied.
"""
function LayRatio(q::Union{Real, AbstractGrid, Gridspace}; combine::Symbol = :product)
    return parameterize(LayRatio, LayRatio, (q,); combine)
end
LayRatio(lay::LayRatio; combine::Symbol = :product) = lay
function LayRatio(values::Union{Tuple, AbstractVector}; combine::Symbol = :product)
    return map(value -> LayRatio(value; combine), Tuple(values))
end
function LayRatio(first, second, remaining...; combine::Symbol = :product)
    return LayRatio((first, second, remaining...); combine)
end

"""
Store an authoritative helical pitch length \\[m\\].
"""
struct Pitch{T <: Real}
    p::T
    Pitch{T}(p::T) where {T <: Real} = validate(new{T}(p))
end

function validate(lay::Pitch)
    isfinite(lay.p) && lay.p > zero(lay.p) ||
        throw(DomainError(lay.p, "pitch must be positive and finite"))
    return lay
end

Pitch(p::Real) = Pitch{typeof(float(p))}(float(p))

"""
$(TYPEDSIGNATURES)

Declare one helical pitch or a homogeneous schedule of helical pitches.

# Arguments

- `p`: one positive helical pitch \\[m\\].
- `values`: two or more pitches, or one tuple or vector of pitches. Specify one pitch per repeated course \\[m\\].

# Keywords

- `combine=:product`: gridspace composition rule.

# Returns

- A `Pitch`, a tuple of `Pitch` values, or the corresponding `Gridspace` when
  an explicit finite source is supplied.
"""
function Pitch(p::Union{Real, AbstractGrid, Gridspace}; combine::Symbol = :product)
    return parameterize(Pitch, Pitch, (p,); combine)
end
Pitch(lay::Pitch; combine::Symbol = :product) = lay
function Pitch(values::Union{Tuple, AbstractVector}; combine::Symbol = :product)
    return map(value -> Pitch(value; combine), Tuple(values))
end
function Pitch(first, second, remaining...; combine::Symbol = :product)
    return Pitch((first, second, remaining...); combine)
end

"""
Store an authoritative helical lay angle relative to the cable axis \\[rad\\].
"""
struct LayAngle{T <: Real}
    α::T
    LayAngle{T}(α::T) where {T <: Real} = validate(new{T}(α))
end

function validate(lay::LayAngle)
    isfinite(lay.α) && zero(lay.α) < lay.α < oftype(lay.α, π / 2) ||
        throw(DomainError(lay.α, "lay angle must lie in (0, π/2)"))
    return lay
end

LayAngle(α::Real) = LayAngle{typeof(float(α))}(float(α))

"""
$(TYPEDSIGNATURES)

Declare one helical lay angle or a homogeneous schedule of lay angles.

# Arguments

- `α`: one lay angle relative to the cable axis \\[rad\\].
- `values`: two or more angles, or one tuple or vector of angles. Specify one angle per repeated course \\[rad\\].

# Keywords

- `combine=:product`: gridspace composition rule.

# Returns

- A `LayAngle`, a tuple of `LayAngle` values, or the corresponding `Gridspace`
  when an explicit finite source is supplied.
"""
function LayAngle(α::Union{Real, AbstractGrid, Gridspace}; combine::Symbol = :product)
    return parameterize(LayAngle, LayAngle, (α,); combine)
end
LayAngle(lay::LayAngle; combine::Symbol = :product) = lay
function LayAngle(values::Union{Tuple, AbstractVector}; combine::Symbol = :product)
    return map(value -> LayAngle(value; combine), Tuple(values))
end
function LayAngle(first, second, remaining...; combine::Symbol = :product)
    return LayAngle((first, second, remaining...); combine)
end

"""
$(TYPEDEF)

Represent one helical path from one authoritative lay definition.

$(TYPEDFIELDS)
"""
struct Helix{L <: Union{LayRatio, Pitch, LayAngle}, T <: Real}
    "Lay definition."
    lay::L
    "Handedness: `1` or `-1` \\[dimensionless\\]."
    dir::Int
    "Initial angular position \\[rad\\]."
    φ0::T

    function Helix{L, T}(lay::L, dir::Int, φ0::T) where {
            L <: Union{LayRatio, Pitch, LayAngle}, T <: Real}
        return validate(new{L, T}(lay, dir, φ0))
    end
end

function validate(helix::Helix)
    helix.dir in (-1, 1) || throw(ArgumentError("helix direction must be -1 or 1"))
    isfinite(helix.φ0) || throw(DomainError(helix.φ0, "helix initial angle must be finite"))
    return helix
end

function Helix(
        lay::Union{LayRatio, Pitch, LayAngle},
        dir::Union{Signed, Unsigned},
        φ0::Real
)
    angle = float(φ0)
    return Helix{typeof(lay), typeof(angle)}(lay, Int(dir), angle)
end

function Helix(
        lay::Union{LayRatio, Pitch, LayAngle, AbstractGrid, Gridspace};
        dir = 1,
        φ0 = 0,
        combine::Symbol = :product
)
    return parameterize(Helix, Helix, (lay, dir, φ0); combine)
end

"""
Return helical pitch at local radius `radius` \\[m\\].
"""
function pitch end
pitch(path::Helix{<:LayRatio}, radius::Real) = path.lay.q * (2 * radius)
pitch(path::Helix{<:Pitch}, radius::Real) = path.lay.p
pitch(path::Helix{<:LayAngle}, radius::Real) = 2π * radius / tan(path.lay.α)

"""
Return helical lay angle at local radius `radius` \\[rad\\].
"""
function angle(path::Helix, radius::Real)
    radius >= zero(radius) || throw(DomainError(radius, "helix radius must be nonnegative"))
    return atan(2π * radius, pitch(path, radius))
end

"""
Return the helical path-length ratio at local radius `radius`. A path at zero radius lies
straight along the axis, with ratio 1.
"""
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
    factor = sqrt(one(ratio) + ((one(ratio) * pi) / ratio)^2)
    return iszero(radius) ? one(factor) : factor
end
