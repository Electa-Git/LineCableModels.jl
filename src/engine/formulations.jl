# Engine-owned formulation hierarchy.
"Identify the native LineCableModels numerical backend."
struct LineCableModelsEngine end

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

"Route an explicit external formulation tag to its `Val` dispatch method."
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)

"""
$(TYPEDEF)

Supertype for formulations that reduce a layered-earth model to equivalent
homogeneous properties.

# Implementations

- [`EnforceLayer`](@ref) selects the properties of one earth layer.
"""
abstract type AbstractEHEMFormulation <: AbstractFormulation end
