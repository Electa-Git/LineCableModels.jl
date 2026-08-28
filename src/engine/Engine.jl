"""
    LineCableModels.Engine

Calculate cable constants and frequency-dependent line-parameter matrices from
materialised cable systems.

# Overview

- Define materialised problems, formulations, and core results.
- Calculate conductor, insulation, and earth-return impedance and admittance.
- Assemble phase-domain series-impedance and shunt-admittance matrices.
- Apply bundle reduction, Kron elimination, transposition, and modal
  transformations.
- Compare, tabulate, and describe plots of completed line-parameter results.

# Dependencies

$(IMPORTS)

"""
module Engine

# Export public API
export LineParametersProblem,
       LineParameters, SeriesImpedance, ShuntAdmittance,
       RMSError, LineParametersBenchmark, compare,
       absolute_error, relative_error,
       Z, Y, R, X, L, G, B, C,
       series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance,
       frequencies, nconductors, nfrequencies, basis,
       kronify
export AbstractFormulation, LineParametersFormulation, Formulation
export LineCableModelsEngine, LineParametersWorkspace
export constitutive
export EarthProperties, CPEarth, verbosity
export InternalImpedance, InsulationImpedance, EarthImpedance
export InsulationAdmittance, EarthAdmittance, EHEM, Transforms

export compute, flatten, plot

# Module-specific dependencies
using LinearAlgebra: I, diag, ldiv!, lu!
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, build, flatten, R, L, C,
                          resistance, inductance, capacitance
import ..LineCableModels: nominal, uncertainty
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  AbstractProblemResult, AbstractCoreResult,
                  FormulationOptions, ComputationOptions,
                  ComputationDetails,
                  formulation_options, computation_options, computation_details, details,
                  compute, observe, @observe, observables, validate_observables, unit_targets,
                  observation_request, observation_indices,
                  request_identity, request_quantity, request_indices,
                  computation_owner

using ..Units
using ..PlotBuilder
using ..Materials
using ..EarthProps: EarthModel
using ..DataModel: CableDesign, LineCableSystem, ncables, nphases
import ..DataModel
import ..LineCableModels: validate
import ..Validation
import Logging
using Logging: AbstractLogger, ConsoleLogger, with_logger
import SpecialFunctions
using QuadGK: quadgk

include("interfaces.jl")
include("formulations.jl")
include("earthkernels.jl")

include("earthproperties/EarthProperties.jl")
using .EarthProperties: CPEarth

# Problem and native formulation definitions
include("problems.jl")
include("options.jl")

# Line-parameter results and their protocols
include("lineparameters/lineparameters.jl")
include("lineparameters/quantities.jl")
include("lineparameters/benchmark.jl")

# Submodule `InternalImpedance`
include("internalimpedance/InternalImpedance.jl")
using .InternalImpedance: InternalImpedance

# Submodule `InsulationImpedance`
include("insulationimpedance/InsulationImpedance.jl")
using .InsulationImpedance: InsulationImpedance

# Submodule `EarthImpedance`
include("earthimpedance/EarthImpedance.jl")
using .EarthImpedance: EarthImpedance

# Submodule `InsulationAdmittance`
include("insulationadmittance/InsulationAdmittance.jl")
using .InsulationAdmittance: InsulationAdmittance

# Submodule `EarthAdmittance`
include("earthadmittance/EarthAdmittance.jl")
using .EarthAdmittance: EarthAdmittance

# Submodule `Transforms`
include("transforms/Transforms.jl")
using .Transforms
using .Transforms: reciprocity!, ideal_transposition!

# Submodule `EHEM`
include("ehem/EHEM.jl")
using .EHEM

# Native workspace and numerical action
include("native_adapter.jl")
include("input.jl")
include("logging.jl")
include("earthreturn.jl")
include("impedance.jl")
include("admittance.jl")
include("compute.jl")
include("reduction.jl")

# Line-parameter protocols and renderer-independent plot definitions
include("lineparameters/base.jl")
include("lineparameters/plot.jl")
include("lineparameters/plotdefinition.jl")
include("lineparameters/comparisonplot.jl")

public has_uncertainty_type

end # module Engine
