"""
$(TYPEDEF)

Abstract type for numerical problem definitions in
[`LineCableModels.jl`](index.md).
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

Construct a formulation selected by `backend`, using the analytical formulation by
default.
"""
Formulation(backend::Symbol; kwargs...) = Formulation(Val(backend); kwargs...)
Formulation(; kwargs...) = Formulation(:analytical; kwargs...)
Formulation(::Val{:FEM}; kwargs...) = retired_fem_sector("FEM formulation")

"""
$(TYPEDEF)

Abstract type representing equivalent homogeneous earth models used by the
multi-dispatch earth-property evaluation.

# Currently available formulations

- [`EnforceLayer`](@ref): Effective parameters defined according to a specific earth layer.
"""
abstract type AbstractEHEMFormulation <: AbstractFormulation end
