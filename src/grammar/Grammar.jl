"""
    LineCableModels.Grammar

Define calculation supertypes and functions shared by Engine,
ParametricBuilder, UQ, and external implementations.

# Public actions

- `formulation_options` and `computation_options` normalise owner-specific options.
- `computation_details` normalises supplemental output from a registered
  computation owner, and `details` reads retained supplemental output.
- `compute` evaluates a problem through a selected formulation.
- `observe` and `@observe` read native numerical values from completed results.
- `observables` publishes explicitly requested scientific values.
"""
module Grammar

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractCoreResult, AbstractResultSpace
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables
export nominal, standard_uncertainty

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF
import ..LineCableModels: basis
using ..Units: UnitExpr, quantity, native_unit, display_unit, scale_factor

include("types.jl")
include("results.jl")
include("interfaces.jl")
include("uncertainty.jl")

public check_core_result, computation_owner, detach

end # module Grammar
