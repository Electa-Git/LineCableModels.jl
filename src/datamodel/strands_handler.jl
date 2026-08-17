"""
$(TYPEDSIGNATURES)

Report that a cable-part type has no packing-limit implementation.
"""
@noinline function maxfill(::Type{T}, args...) where {T}
    throw(ArgumentError(
        "maxfill is not implemented for $(nameof(T)) with geometry " *
        "$(map(typeof, args))",
    ))
end
