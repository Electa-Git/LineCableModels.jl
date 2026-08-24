# Engine-owned formulation hierarchy.
"""
$(TYPEDEF)

Supertype for Engine impedance formulations.
"""
abstract type AbstractImpedanceFormulation <: AbstractFormulation end
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

abstract type AbstractAdmittanceFormulation <: AbstractFormulation end
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

abstract type AbstractTransformFormulation <: AbstractFormulation end
abstract type AbstractEarthPropertiesFormulation <: AbstractFormulation end

"""
$(TYPEDSIGNATURES)

Construct the formulation selected by `backend`. Calling `Formulation()`
selects `:analytical`.
"""
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)
Formulation(; kwargs...) = Formulation(:analytical; kwargs...)
Formulation(::Val{:FEM}; kwargs...) = retired_fem_sector("FEM formulation")

"""
$(TYPEDEF)

Supertype for formulations that reduce a layered-earth model to equivalent
homogeneous properties.

# Implementations

- [`EnforceLayer`](@ref) selects the properties of one earth layer.
"""
abstract type AbstractEHEMFormulation <: AbstractFormulation end
