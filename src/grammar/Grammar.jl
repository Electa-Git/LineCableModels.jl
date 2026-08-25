"""
    LineCableModels.Grammar

Define calculation supertypes and functions shared by Engine,
ParametricBuilder, UQ, and external implementations.

# Public actions

- `formulation_options` and `computation_options` normalise owner-specific options.
- `compute` evaluates a problem through a selected formulation.
- `observe` reads native numerical values from completed results.
- `observables` publishes explicitly requested scientific values.
- `primitives` and `preprocess` are unimplemented extension points for future
  higher-order calculations.
"""
module Grammar

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables, primitives, preprocess
export nominal, standard_uncertainty

using DocStringExtensions: SIGNATURES, TYPEDEF
import ..LineCableModels: basis
using ..Units: UnitExpr, quantity, native_unit, display_unit, scale_factor

include("types.jl")
include("interfaces.jl")
include("uncertainty.jl")

end # module Grammar
