"""
	LineCableModels.Materials

The [`Materials`](@ref) module provides functionality for managing and utilizing material properties within the [`LineCableModels.jl`](index.md) package. This module includes definitions for material properties, a library for storing and retrieving materials, and functions for manipulating material data.

# Overview

- Defines the [`Material`](@ref) struct representing fundamental physical properties of materials.
- Provides the [`MaterialsLibrary`](@ref) mutable struct for storing a collection of materials.
- Includes functions for adding, removing, and retrieving materials from the library.
- Supports loading and saving material data from/to JSON files.
- Contains utility functions for displaying material data.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Materials

# Export public API
export Material, MaterialsLibrary

# Load common dependencies
using ..LineCableModels
include("utils/commondeps.jl")


# Module-specific dependencies
using Measurements
using DataFrames
using ..Utils

"""
$(TYPEDEF)

Defines electromagnetic and thermal properties of a material used in cable modeling:

$(TYPEDFIELDS)
"""
struct Material{T<:REALSCALAR}
    "Electrical resistivity of the material \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T
    "Reference temperature for property evaluations \\[°C\\]."
    T0::T
    "Temperature coefficient of resistivity \\[1/°C\\]."
    alpha::T

    @inline function Material{T}(rho::T, eps_r::T, mu_r::T, T0::T, alpha::T) where {T<:REALSCALAR}
        return new{T}(rho, eps_r, mu_r, T0, alpha)
    end

end

"""
$(TYPEDSIGNATURES)

Weakly-typed constructor that infers the target scalar type `T` from the arguments,
coerces values to `T`, and calls the strict numeric kernel.

# Arguments
- `rho`: Resistivity \\[Ω·m\\].
- `eps_r`: Relative permittivity \\[1\\].
- `mu_r`: Relative permeability \\[1\\].
- `T0`: Reference temperature \\[°C\\].
- `alpha`: Temperature coefficient of resistivity \\[1/°C\\].

# Returns
- `Material{T}` where `T = resolve_T(rho, eps_r, mu_r, T0, alpha)`.
"""
@inline function Material(rho, eps_r, mu_r, T0, alpha)
    T = resolve_T(rho, eps_r, mu_r, T0, alpha)
    return Material{T}(
        coerce_to_T(rho, T),
        coerce_to_T(eps_r, T),
        coerce_to_T(mu_r, T),
        coerce_to_T(T0, T),
        coerce_to_T(alpha, T),
    )
end

include("materials/materialslibrary.jl")
include("materials/dataframe.jl")
include("materials/base.jl")
include("materials/typecoercion.jl")

end # module Materials
