"""
$(TYPEDSIGNATURES)

Validate `element` as the element type of the result space being built.

The element type must be concrete and cannot itself be a result-space
envelope. Concrete external result types are accepted without requiring them
to subtype [`AbstractCoreResult`](@ref).

# Arguments

- `element`: proposed result-space element type.
- The result-space type being built, or `AbstractResultSpace` for a collection of
  results.

# Returns

- `element` when it satisfies the result-space invariant.

# Errors

- `ArgumentError`: `element` is abstract, is `Any` or subtypes
  [`AbstractResultSpace`](@ref).
"""
function validate(element::Type{T}, ::Type{<:AbstractResultSpace}) where {T}
    isconcretetype(element) || throw(ArgumentError(
        "result-space element type must be concrete; got $element",
    ))
    element <: AbstractResultSpace && throw(ArgumentError(
        "a result space cannot contain another result-space envelope",
    ))
    return element
end

#! explicit-imports: off
# Base's iterator trait protocol exposes these values without public bindings.
Base.IteratorSize(::Type{<:AbstractResultSpace}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:AbstractResultSpace}) = Base.HasEltype()
#! explicit-imports: on
Base.eltype(::Type{<:AbstractResultSpace{T}}) where {T} = T
