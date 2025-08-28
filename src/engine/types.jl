"""
$(TYPEDEF)

Abstract base type for all problem definitions in the [`LineCableModels.jl`](index.md) computation framework.
"""
abstract type ProblemDefinition end

# Formulation abstract types
abstract type AbstractFormulationSet end
abstract type AbstractFormulationOptions end

abstract type AbstractImpedanceFormulation <: AbstractFormulationSet end
abstract type InternalImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type InsulationImpedanceFormulation <: AbstractImpedanceFormulation end
abstract type EarthImpedanceFormulation <: AbstractImpedanceFormulation end

abstract type AbstractAdmittanceFormulation <: AbstractFormulationSet end
abstract type InsulationAdmittanceFormulation <: AbstractAdmittanceFormulation end
abstract type EarthAdmittanceFormulation <: AbstractAdmittanceFormulation end

"""
$(TYPEDEF)

Abstract type representing different equivalent homogeneous earth models (EHEM). Used in the multi-dispatch implementation of [`_calc_ehem_properties!`](@ref).

# Currently available formulations

- [`EnforceLayer`](@ref): Effective parameters defined according to a specific earth layer.
"""
abstract type AbstractEHEMFormulation end