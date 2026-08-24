"""
    LineCableModels.Grammar

Define calculation supertypes and functions shared by Engine,
ParametricBuilder, UQ, and external implementations.

# Public actions

- `formulation_options` and `computation_options` normalise owner-specific options.
- `compute` evaluates a problem through a selected formulation.
- `observables` publishes immutable scientific result data.
- `primitives` and `preprocess` are unimplemented extension points for future
  higher-order calculations.
"""
module Grammar

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions
export formulation_options, computation_options
export compute, observables, primitives, preprocess
export nominal, standard_uncertainty

using DocStringExtensions: SIGNATURES, TYPEDEF

include("types.jl")
include("interfaces.jl")
include("uncertainty.jl")

end # module Grammar
