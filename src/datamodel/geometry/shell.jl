"""
$(TYPEDEF)

Declare an outward conformal layer relative to an exact resolved boundary.

`Shell` is contextual construction input rather than an intrinsic primitive.
Calling [`resolve`](@ref) against the preceding boundary produces the exact
material domain occupied by the layer.

$(TYPEDFIELDS)
"""
struct Shell{T <: Real}
    "Normal layer thickness \\[m\\]."
    t::T

    function Shell{T}(t::T) where {T <: Real}
        isfinite(t) && t > zero(t) ||
            throw(DomainError(t, "shell thickness must be positive and finite"))
        return new{T}(t)
    end
end

Shell(t::Real) = Shell{typeof(float(t))}(float(t))

Base.eltype(::Shell{T}) where {T} = T
Base.eltype(::Type{<:Shell{T}}) where {T} = T

function Base.convert(::Type{<:Shell{T}}, value::Shell) where {T <: Real}
    return Shell{T}(convert(T, value.t))
end
