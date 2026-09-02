"Store a helical lay-length to mean-diameter ratio."
struct LayRatio{T <: Real}
    q::T
    function LayRatio{T}(q::T) where {T <: Real}
        isfinite(q) && q > zero(q) ||
            throw(DomainError(q, "lay ratio must be positive and finite"))
        return new{T}(q)
    end
end

_lay_ratio(q) = LayRatio{typeof(float(q))}(float(q))

function _normalize_schedule(wrapper, values; combine::Symbol)
    return map(Tuple(values)) do value
        value isa wrapper ? value : wrapper(value; combine)
    end
end

"""
$(TYPEDSIGNATURES)

Declare one lay-length ratio or a homogeneous schedule of lay-length ratios.

# Arguments

- `q`: One positive lay-length to mean-diameter ratio \\[dimensionless\\].
- `values`: Two or more ratios, or one tuple or vector of ratios, defining one
  ratio per repeated course.

# Keywords

- `combine=:product`: Gridspace composition rule.

# Returns

- A `LayRatio`, a tuple of `LayRatio` values, or the corresponding `Gridspace`
  when an explicit finite source is supplied.
"""
LayRatio(q; combine::Symbol = :product) =
    _construction(LayRatio, _lay_ratio, (q,); combine)
LayRatio(values::Union{Tuple, AbstractVector}; combine::Symbol = :product) =
    _normalize_schedule(LayRatio, values; combine)
LayRatio(first, second, remaining...; combine::Symbol = :product) =
    _normalize_schedule(LayRatio, (first, second, remaining...); combine)

"Store an authoritative helical pitch length \\[m\\]."
struct Pitch{T <: Real}
    p::T
    function Pitch{T}(p::T) where {T <: Real}
        isfinite(p) && p > zero(p) ||
            throw(DomainError(p, "pitch must be positive and finite"))
        return new{T}(p)
    end
end

_pitch(p) = Pitch{typeof(float(p))}(float(p))

"""
$(TYPEDSIGNATURES)

Declare one helical pitch or a homogeneous schedule of helical pitches.

# Arguments

- `p`: One positive helical pitch \\[m\\].
- `values`: Two or more pitches, or one tuple or vector of pitches, defining
  one pitch per repeated course \\[m\\].

# Keywords

- `combine=:product`: Gridspace composition rule.

# Returns

- A `Pitch`, a tuple of `Pitch` values, or the corresponding `Gridspace` when
  an explicit finite source is supplied.
"""
Pitch(p; combine::Symbol = :product) =
    _construction(Pitch, _pitch, (p,); combine)
Pitch(values::Union{Tuple, AbstractVector}; combine::Symbol = :product) =
    _normalize_schedule(Pitch, values; combine)
Pitch(first, second, remaining...; combine::Symbol = :product) =
    _normalize_schedule(Pitch, (first, second, remaining...); combine)

"Store an authoritative helical lay angle relative to the cable axis \\[rad\\]."
struct LayAngle{T <: Real}
    α::T
    function LayAngle{T}(α::T) where {T <: Real}
        isfinite(α) && zero(α) < α < oftype(α, π / 2) ||
            throw(DomainError(α, "lay angle must lie in (0, π/2)"))
        return new{T}(α)
    end
end

_lay_angle(α) = LayAngle{typeof(float(α))}(float(α))

"""
$(TYPEDSIGNATURES)

Declare one helical lay angle or a homogeneous schedule of lay angles.

# Arguments

- `α`: One lay angle relative to the cable axis \\[rad\\].
- `values`: Two or more angles, or one tuple or vector of angles, defining one
  angle per repeated course \\[rad\\].

# Keywords

- `combine=:product`: Gridspace composition rule.

# Returns

- A `LayAngle`, a tuple of `LayAngle` values, or the corresponding `Gridspace`
  when an explicit finite source is supplied.
"""
LayAngle(α; combine::Symbol = :product) =
    _construction(LayAngle, _lay_angle, (α,); combine)
LayAngle(values::Union{Tuple, AbstractVector}; combine::Symbol = :product) =
    _normalize_schedule(LayAngle, values; combine)
LayAngle(first, second, remaining...; combine::Symbol = :product) =
    _normalize_schedule(LayAngle, (first, second, remaining...); combine)

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

    function Helix{L, T}(lay::L, dir::Int, φ0::T) where {L, T <: Real}
        lay isa Union{LayRatio, Pitch, LayAngle} ||
            throw(ArgumentError("helix lay must be LayRatio, Pitch, or LayAngle"))
        dir in (-1, 1) || throw(ArgumentError("helix direction must be -1 or 1"))
        isfinite(φ0) || throw(DomainError(φ0, "helix initial angle must be finite"))
        return new{L, T}(lay, dir, φ0)
    end
end

function _helix(lay, dir, φ0)
    dir isa Integer && !(dir isa Bool) || throw(ArgumentError(
        "helix direction must be an integer"
    ))
    angle = float(φ0)
    return Helix{typeof(lay), typeof(angle)}(lay, Int(dir), angle)
end

function Helix(
        lay;
        dir = 1,
        φ0 = 0,
        combine::Symbol = :product
)
    return _construction(Helix, _helix, (lay, dir, φ0); combine)
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
