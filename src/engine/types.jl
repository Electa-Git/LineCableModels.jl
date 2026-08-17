"""
$(TYPEDEF)

Abstract type for numerical problem definitions in
[`LineCableModels.jl`](index.md).
"""
abstract type ProblemDefinition end

# Formulation abstract types
"""
$(TYPEDEF)

Abstract type for formulations that select physical and numerical methods.
"""
abstract type AbstractFormulation end

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

Construct a formulation selected by `engine`, using the EMT formulation by
default.
"""
Formulation(engine::Symbol; kwargs...) = Formulation(Val(engine); kwargs...)
Formulation(; kwargs...) = Formulation(Val(:EMT); kwargs...)
Formulation(::Val{:FEM}; kwargs...) = retired_fem_sector("FEM formulation")

"""
$(TYPEDEF)

Abstract type representing equivalent homogeneous earth models used by the
multi-dispatch earth-property evaluation.

# Currently available formulations

- [`EnforceLayer`](@ref): Effective parameters defined according to a specific earth layer.
"""
abstract type AbstractEHEMFormulation <: AbstractFormulation end

abstract type AbstractFormulationOptions end
