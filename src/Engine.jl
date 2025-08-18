"""
    LineCableModels.Engine

The [`Engine`](@ref) module provides the main functionalities of the [`LineCableModels.jl`](index.md) package. This module implements data structures, methods and functions for calculating frequency-dependent electrical parameters (Z/Y matrices) of line and cable systems with uncertainty quantification. 

# Overview

- Calculation of frequency-dependent series impedance (Z) and shunt admittance (Y) matrices.
- Uncertainty propagation for geometric and material parameters using `Measurements.jl`.
- Internal impedance computation for solid, tubular and multi-layered coaxial conductors.
- Earth return impedances/admittances for overhead lines and underground cables (valid up to 10 MHz).
- Support for frequency-dependent soil properties.
- Handling of arbitrary polyphase systems with multiple conductors per phase.
- Phase and sequence domain calculations with uncertainty quantification.
- Novel N-layer concentric cable formulation with semiconductor modeling.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Engine

# Export public API
export LineParametersProblem, LineParameters
export CoaxialFormulation, FormulationSet

export compute!

# Load common dependencies
include("commondeps.jl")
using ..Utils
using ..Materials
using ..EarthProps
using ..DataModel
using ..LineCableModels # For physical constants (f₀, μ₀, ε₀, ρ₀, T₀, TOL, ΔTmax)
import ..LineCableModels: FormulationSet, _get_description, REALTYPES, COMPLEXTYPES, NUMERICTYPES

# Module-specific dependencies
using Measurements
using LinearAlgebra
using SpecialFunctions

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

# Problem definitions
include("Engine/problemdefs.jl")

# Submodule `InternalImpedance`
include("Engine/InternalImpedance.jl")
@force using .InternalImpedance

# Submodule `InsulationImpedance`
include("Engine/InsulationImpedance.jl")
@force using .InsulationImpedance

# Submodule `EarthImpedance`
include("Engine/EarthImpedance.jl")
@force using .EarthImpedance

# Submodule `InsulationAdmittance`
include("Engine/InsulationAdmittance.jl")
@force using .InsulationAdmittance

# Submodule `EarthAdmittance`
include("Engine/EarthAdmittance.jl")
@force using .EarthAdmittance

# Submodule `EHEM`
include("Engine/EHEM.jl")
@force using .EHEM

# Helpers
include("Engine/utils.jl")

# Workspace definition
include("Engine/workspace.jl")

# Computation methods
include("Engine/solver.jl")

# Override I/O methods
include("Engine/io.jl")

@reexport using .InternalImpedance, .InsulationImpedance, .EarthImpedance,
    .InsulationAdmittance, .EarthAdmittance, .EHEM

end # module Engine