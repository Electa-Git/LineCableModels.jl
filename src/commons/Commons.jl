"""
    LineCableModels.Commons

Define computation supertypes and functions shared by Engine,
ParametricBuilder, UQ, and external implementations.

# Physical constants

- `vacuum_permittivity` returns the vacuum permittivity in a requested scalar type.

# Public actions

- `formulation_options` and `computation_options` normalize owner-specific options.
- `formulas` lists the identifiers that a formula family registers.
- `initialize_buffers` builds the reusable arrays of a computation into its buffer
  record, one method per formula or shared component.
- `computation_details` normalizes supplemental output from a registered
  computation owner, and `details` reads retained supplemental output.
- `validate` checks a materialized input before its consumer uses it and returns it
  unchanged.
- `compute` evaluates a problem through a selected formulation, each a scalar or a space,
  through the runner `traverse`. `resolve` gives a backend what one formulation needs from
  a problem as a `Resolved` value, and `retain_gridpoint` records the identity of a point
  in its result. `ParametricProblem`, `Combinatorial` and `ParametricResult` declare and
  hold a parametric run. `on_error` decides whether a study records a failed point or stops.
- `verbosity` and `VerbosityLogger` select and filter log records by execution verbosity.
- `Grid` declares a finite source, `Gridspace` a lazy space of finite sources combined
  by `Product` or `Zip`, and `parameterize` lifts a call over them. `has_grid` tells
  whether a value contains a source. `ScalarOrGrid{T}` is a `T` or a space of `T`.
- `observe` and `@observe` read native numerical values from completed results.
- `observables` publishes explicitly requested scientific values.
- `validate_observables` and `unit_targets` align publication requests and
  display units for presentation consumers.
"""
module Commons

export AbstractProblemDefinition, AbstractFormulation, AbstractTraversal, AbstractProblemResult
export AbstractCoreResult, AbstractResultSpace
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables, validate
export nominal, uncertainty
export add!, build, homogenize, constitutive, line_length, description
export R, L, C, resistance, inductance, capacitance
export formula, formula_id
export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace, has_uncertainty, has_grid
export Combinatorial, ParametricProblem, ParametricResult

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using RequiredInterfaces: @required
import Logging
import UUIDs
import Random
import ..Units
using ..Units: UnitExpr, quantity, native_unit, display_unit, scale_factor

include("consts.jl")
include("types.jl")
include("logging.jl")
include("base.jl")
include("results.jl")
include("interfaces.jl")
include("formulas.jl")
include("selections.jl")
include("observables.jl")
include("uncertainty.jl")
include("grid.jl")
include("gridspace.jl")
include("gridpoint.jl")
include("parametric.jl")
include("runner.jl")
include("observed_validation.jl")
include("observedresult.jl")
include("retained_products.jl")

public vacuum_permittivity
public initialize_buffers
public FormulaDefinition, Expression, Functor, formulas
public validate_observables, unit_targets, detach
public observation_request, observation_indices, materialize_observation
public observation_resolution
public input_fields, observation_gridpoint, observation_requests, observation_quantity
export ObservedResult
public observation_groups, observation_labels, observation_product, gridpoint_id
public observation_selection
public request_identity, request_quantity, request_indices
public normalize_observation_selector
public Gridpoint, Composition, Product, Zip, ScalarOrGrid, parameterize, materialize, points
public Builder, FormulaBuilder, target
public realize, realize_arguments, sample_uncertainty
public resolve, Resolved, traverse, retain_gridpoint, on_error
end # module Commons
