"""
$(TYPEDEF)

Supertype for materialised cable parts.
"""
abstract type AbstractCablePart{T} end

"""
$(TYPEDEF)

Supertype for materialised conductive cable parts.

Direct implementations include:
- [`Tubular`](@ref)
- [`Strip`](@ref)
"""
abstract type AbstractConductorPart{T} <: AbstractCablePart{T} end

"""
$(TYPEDEF)

Supertype for conductor parts composed of discrete strands.

Direct implementations include:
- [`CircStrands`](@ref)
- [`RectStrands`](@ref)
"""
abstract type AbstractStrandsLayer{T} <: AbstractConductorPart{T} end

"""
$(TYPEDEF)

Supertype for materialised dielectric cable parts.

Direct implementations include:
- [`Insulator`](@ref)
- [`Semicon`](@ref)
"""
abstract type AbstractInsulatorPart{T} <: AbstractCablePart{T} end

Base.eltype(::AbstractCablePart{T}) where {T} = T
Base.eltype(::Type{<:AbstractCablePart{T}}) where {T} = T

### Provisions for the new rectangular-strand geometry
abstract type AbstractShapeGeometry end
