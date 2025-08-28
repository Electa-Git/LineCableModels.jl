# To handle radius-related operations
abstract type AbstractRadius <: Number end

"""
$(TYPEDEF)

Represents the thickness of a cable component.

$(TYPEDFIELDS)
"""
struct Thickness{T<:Real} <: AbstractRadius
    "Numerical value of the thickness \\[m\\]."
    value::T
    function Thickness(value::T) where {T<:Real}
        value >= 0 || throw(ArgumentError("Thickness must be a non-negative number."))
        new{T}(value)
    end
end

"""
$(TYPEDEF)

Represents the diameter of a cable component.

$(TYPEDFIELDS)
"""
struct Diameter{T<:Real} <: AbstractRadius
    "Numerical value of the diameter \\[m\\]."
    value::T
    function Diameter(value::T) where {T<:Real}
        value > 0 || throw(ArgumentError("Diameter must be a positive number."))
        new{T}(value)
    end
end

"""
$(TYPEDEF)

Abstract type representing a generic cable part.
"""
abstract type AbstractCablePart{T} end

"""
$(TYPEDEF)

Abstract type representing a conductive part of a cable.

Subtypes implement specific configurations:
- [`WireArray`](@ref)
- [`Tubular`](@ref)
- [`Strip`](@ref)
"""
abstract type AbstractConductorPart{T} <: AbstractCablePart{T} end
abstract type AbstractWireArray{T} <: AbstractConductorPart{T} end

"""
$(TYPEDEF)

Abstract type representing an insulating part of a cable.

Subtypes implement specific configurations:
- [`Insulator`](@ref)
- [`Semicon`](@ref)
"""
abstract type AbstractInsulatorPart{T} <: AbstractCablePart{T} end


# If a correct ctor exists, Julia will pick it; this runs only when arity is wrong.
function (::Type{T})(args::Vararg{Any,N}; kwargs...) where {T<:AbstractCablePart,N}
    throw(ArgumentError("[$(nameof(T))] constructor: invalid number of positional args ($N)."))
end