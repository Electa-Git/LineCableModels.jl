"""
$(TYPEDEF)

Abstract type representing a generic cable part.
"""
abstract type AbstractCablePart{T} end

"""
$(TYPEDEF)

Abstract type representing a conductive part of a cable.

Subtypes implement specific configurations:
- [`Tubular`](@ref)
- [`Strip`](@ref)
"""
abstract type AbstractConductorPart{T} <: AbstractCablePart{T} end

"""
$(TYPEDEF)

Abstract type representing all stranded configurations composed of grouped discrete geometric shapes.

Subtypes implement specific configurations:
- [`CircStrands`](@ref)
- [`RectStrands`](@ref)
"""
abstract type AbstractStrandsLayer{T} <: AbstractConductorPart{T} end

"""
$(TYPEDEF)

Abstract type representing an insulating part of a cable.

Subtypes implement specific configurations:
- [`Insulator`](@ref)
- [`Semicon`](@ref)
"""
abstract type AbstractInsulatorPart{T} <: AbstractCablePart{T} end

# If a correct ctor exists, Julia will pick it; this runs only when arity is wrong.
function (::Type{T})(args::Vararg{Any, N}; kwargs...) where {T <: AbstractCablePart, N}
    throw(
        ArgumentError(
        "[$(nameof(T))] constructor: invalid number of positional args ($N).",
    ),
    )
end

### Provisions for the new rectangular-strand geometry
abstract type AbstractShapeGeometry end
