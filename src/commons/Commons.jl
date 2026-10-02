"""
    LineCableModels.Commons

Define calculation supertypes and functions shared by Engine,
ParametricBuilder, UQ, and external implementations.

# Physical constants

- `vacuum_permittivity` and `vacuum_permeability` return the vacuum constants in
  a requested scalar type.

# Matrix reductions

- `ReductionPlan` fixes the terminal reorder, bundle merge, Kron elimination and
  ideal transposition of primitive line matrices.
- `reduce_line_matrices!` applies a plan to one frequency and inverts the potential
  coefficients to the shunt admittance, in `ReductionBuffers`.

# Public actions

- `formulation_options` and `computation_options` normalize owner-specific options.
- `computation_details` normalizes supplemental output from a registered
  computation owner, and `details` reads retained supplemental output.
- `compute` evaluates a problem through a selected formulation.
- `observe` and `@observe` read native numerical values from completed results.
- `observables` publishes explicitly requested scientific values.
- `validate_observables` and `unit_targets` align publication requests and
  display units for presentation consumers.
"""
module Commons

export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractCoreResult, AbstractResultSpace
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables
export nominal, uncertainty

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
import ..LineCableModels: basis
import UUIDs
import Random
import ..Units
using LinearAlgebra: I, axpy!, checksquare, cond, diag, ldiv!, lu!, mul!, norm
using ..Units: UnitExpr, quantity, native_unit, display_unit, scale_factor

include("consts.jl")
include("matrixops.jl")
include("types.jl")
include("base.jl")
include("results.jl")
include("interfaces.jl")
include("formulas.jl")
include("observables.jl")
include("uncertainty.jl")
include("gridpoint.jl")
include("observed_validation.jl")
include("observedresult.jl")
include("retained_products.jl")

public vacuum_permittivity, vacuum_permeability
public ideal_transposition!, reorder_indices, kron_reduce, kron_reduce!,
       bundle_operations, merge_bundles!
public ReductionPlan, ReductionBuffers, reduce_line_matrices!
public check_core_result
public FormulaDefinition, FormulaMethod
public validate_observables, unit_targets, detach
public observation_request, observation_indices, materialize_observation
public observation_resolution
public input_fields, observation_gridpoint, observation_requests, observation_quantity
export ObservedResult
public observation_groups, observation_labels, observation_product, gridpoint_id
public observation_selection
public request_identity, request_quantity, request_indices
public normalize_observation_selector
end # module Commons
