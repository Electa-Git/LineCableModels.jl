"""
$(TYPEDEF)

Represent planar translation and orientation relative to a parent frame.

$(TYPEDFIELDS)
"""
struct Pose2{T <: Real}
    "Horizontal translation \\[m\\]."
    x::T
    "Vertical translation \\[m\\]."
    y::T
    "Counter-clockwise orientation \\[rad\\]."
    φ::T

    function Pose2{T}(x::T, y::T, φ::T) where {T <: Real}
        all(isfinite, (x, y, φ)) ||
            throw(ArgumentError("pose coordinates and orientation must be finite"))
        return new{T}(x, y, φ)
    end
end

function _pose2(x, y, φ)
    values = map(float, promote(x, y, φ))
    return Pose2{typeof(first(values))}(values...)
end

function Pose2(x, y, φ; combine::Symbol = :product)
    return parameterize(Pose2, _pose2, (x, y, φ); combine)
end

Pose2(x, y; combine::Symbol = :product) = Pose2(x, y, 0; combine)

Base.eltype(::Pose2{T}) where {T} = T
Base.eltype(::Type{<:Pose2{T}}) where {T} = T

function Base.convert(::Type{<:Pose2{T}}, pose::Pose2) where {T <: Real}
    return Pose2{T}(convert(T, pose.x), convert(T, pose.y), convert(T, pose.φ))
end

"Compose `child` relative to `parent`."
function Base.:*(parent::Pose2, child::Pose2)
    values = map(float, promote(
        parent.x, parent.y, parent.φ, child.x, child.y, child.φ
    ))
    px, py, pφ, cx, cy, cφ = values
    return Pose2(
        px + cos(pφ) * cx - sin(pφ) * cy,
        py + sin(pφ) * cx + cos(pφ) * cy,
        pφ + cφ
    )
end
