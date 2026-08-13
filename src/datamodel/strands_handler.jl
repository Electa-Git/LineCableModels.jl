struct MaxFill end # The promise proxy

import ..Validation: maxfill

"""
$(TYPEDSIGNATURES)

Fallback method for the [`maxfill`](@ref) interface. Throws an explicit error indicating that the component type `T` has not implemented its physical capacity limit.

# Arguments

- `::Type{T}`: Component type \\[dimensionless\\].
- `args...`: Geometric parameters required for the calculation \\[dimensionless\\].

# Returns

- Nothing. Always throws an `ArgumentError`.
"""
@noinline function maxfill(::Type{T}, args...) where {T}
	throw(
		ArgumentError(
			"[$(_typename(T))] `maxfill` is not implemented. Any component using `MaxFill()` " *
			"or the `PhysicalFillLimit` rule must overload `maxfill(::Type{$(_typename(T))}, args...)`.",
		),
	)
end

# The plumbing (Never needs to be touched again)
@inline _resolve_strands(x::Int, ::Type{T}, args...) where {T} = x
@inline _resolve_strands(::MaxFill, ::Type{T}, args...) where {T} = maxfill(T, args...)
