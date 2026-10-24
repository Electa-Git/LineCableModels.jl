"""
$(TYPEDEF)

Supertype for materialized cable parts.
"""
abstract type AbstractCablePart end

# A constructor call that holds a Grid or a Gridspace, and that no scalar method accepts,
# lifts through the grammar.
function (::Type{T})(arguments...; combine = :product) where {T <: AbstractCablePart}
    return parameterize(T, arguments; combine)
end
