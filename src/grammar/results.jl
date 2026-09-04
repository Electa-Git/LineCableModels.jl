"""
$(TYPEDSIGNATURES)

Check whether `T` can be the element type of a result space.

The element type must be concrete and cannot itself be a result-space
envelope. Concrete external result types are accepted without requiring them
to subtype [`AbstractCoreResult`](@ref).

# Arguments

- `T`: Proposed result-space element type.

# Returns

- `nothing` when `T` satisfies the result-space invariant.

# Errors

- `ArgumentError`: `T` is abstract, is `Any`, or subtypes
  [`AbstractResultSpace`](@ref).
"""
function check_core_result(::Type{T}) where {T}
    isconcretetype(T) || throw(ArgumentError(
        "result-space element type must be concrete; got $T",
    ))
    T <: AbstractResultSpace && throw(ArgumentError(
        "a result space cannot contain another result-space envelope",
    ))
    return nothing
end

#! explicit-imports: off
# Base's iterator trait protocol exposes these values without public bindings.
Base.IteratorSize(::Type{<:AbstractResultSpace}) = Base.HasShape{1}()
Base.IteratorEltype(::Type{<:AbstractResultSpace}) = Base.HasEltype()
#! explicit-imports: on
Base.eltype(::Type{<:AbstractResultSpace{T}}) where {T} = T
