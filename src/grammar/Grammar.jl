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
- `validate_observables` and `unit_targets` align publication requests and
  display units for presentation consumers.
"""
module Grammar

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractCoreResult, AbstractResultSpace
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables
export nominal, uncertainty
export @orchestrator

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF
import ..LineCableModels: basis
import Tables
import ..Units
using ..Units: UnitExpr, quantity, native_unit, display_unit, scale_factor

include("types.jl")
include("orchestrators.jl")
include("results.jl")
include("interfaces.jl")
include("observables.jl")
include("uncertainty.jl")

public check_core_result, computation_owner
public validate_observables, unit_targets, detach
public observation_request, observation_indices, materialize_observation
public request_identity, request_quantity, request_indices
public ObservationPublication, publication_table
public orchestrator_root, orchestrator_method

end # module Grammar
