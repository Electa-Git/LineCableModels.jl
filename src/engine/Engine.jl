"""
    LineCableModels.Engine

Calculate cable constants and frequency-dependent line-parameter matrices from
completed cable declarations and Engine-owned numerical blueprints.

# Overview

- Define scalar problems, formulations, and core results.
- Calculate conductor, insulation, and earth-return impedance and admittance.
- Assemble phase-domain series-impedance and shunt-admittance matrices.
- Apply bundle reduction, Kron elimination, and ideal transposition.
- Compare, tabulate, and describe plots of completed line-parameter results.

# Dependencies

$(IMPORTS)

"""
module Engine
import ..Grammar: ObservationPublication

# Export public API
export LineParametersProblem, CableConstantsProblem,
       LineParameters, CableConstants, SeriesImpedance, ShuntAdmittance,
       RMSError, LineParametersBenchmark, compare,
       absolute_error, relative_error,
       Z, Y, R, X, L, G, B, C,
       series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance,
       frequencies, nconductors, nfrequencies, basis,
       kronify
export AbstractFormulation, LineParametersFormulation, CableConstantsFormulation,
       Formulation
export LineCableModelsCoaxial, LineCableModelsFEM,
       LineCableModelsFEMError, LineParametersWorkspace
export constitutive, formula_id, EarthPair
export verbosity
export InternalImpedance, InsulationImpedance, EarthImpedance, PipeImpedance
export InsulationAdmittance, SemiconAdmittance, EarthAdmittance
export ShuntModel, BoundarySolveError

export compute

# Module-specific dependencies
using LinearAlgebra: I, checksquare, diag, ldiv!, lu!, mul!
import LinearAlgebra: norm
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, build, R, L, C,
                          resistance, inductance, capacitance
import ..LineCableModels: nominal
import ..LineCableModels: constitutive, formula, formula_id,
                          FormulaMethod, FormulaDefinition
import ..LineCableModels: parameterize
import ..LineCableModels: performance_sample_active
#! explicit-imports: off
import ..LineCableModels: description
#! explicit-imports: on
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  AbstractProblemResult, AbstractCoreResult,
                  FormulationOptions, ComputationOptions,
                  ComputationDetails,
                  formulation_options, computation_options, computation_details, details,
                  compute, observe, observables,
                  observation_request, observation_indices, observation_resolution,
                  uncertainty,
                  request_identity, request_indices,
                  publication_table

using ..Units
import ..Grammar
using ..Materials
using ..Materials: TemperatureDependent
import ..Earth
using ..Earth: EarthMaterial, EarthModel, EquivalentHomogeneous
using ..DataModel: CableDesign, LineCableSystem, ncables, nphases
import ..DataModel
import ..TextDisplay
import ..LineCableModels: validate
import Logging
using Logging: AbstractLogger, ConsoleLogger, with_logger
import SpecialFunctions
using QuadGK: alloc_segbuf, quadgk

include("interfaces.jl")
include("formulations.jl")
include("specialfunctions.jl")

# Problem and coaxial formulation definitions
include("problems.jl")
include("options.jl")
include("integration.jl")

# Line-parameter results and their protocols
include("lineparameters/lineparameters.jl")
include("lineparameters/quantities.jl")
include("lineparameters/resolution.jl")
include("lineparameters/benchmark.jl")
include("matrixops.jl")

# Submodule `InternalImpedance`
include("internalimpedance/InternalImpedance.jl")
using .InternalImpedance: InternalImpedance

include("pipeimpedance/PipeImpedance.jl")

# Submodule `InsulationImpedance`
include("insulationimpedance/InsulationImpedance.jl")
using .InsulationImpedance: InsulationImpedance

# Submodule `EarthImpedance`
include("earthimpedance/EarthImpedance.jl")
using .EarthImpedance: EarthImpedance

# Submodule `InsulationAdmittance`
include("insulationadmittance/InsulationAdmittance.jl")
using .InsulationAdmittance: InsulationAdmittance

# Submodule `SemiconAdmittance`
include("semiconadmittance/SemiconAdmittance.jl")
using .SemiconAdmittance: SemiconAdmittance

# Submodule `EarthAdmittance`
include("earthadmittance/EarthAdmittance.jl")
using .EarthAdmittance: EarthAdmittance

# Native workspace and numerical action
include("blueprint.jl")
include("shuntmodel/ShuntModel.jl")
using .ShuntModel: BoundarySolveError
include("blueprint_shunt.jl")
include("input.jl")
include("logging.jl")
include("earthreturn.jl")
include("impedance.jl")
include("admittance.jl")
include("lineparameters.jl")
include("reduction.jl")
include("cableconstants.jl")

# Line-parameter protocols and observation publication
include("lineparameters/base.jl")
include("lineparameters/publication.jl")
include("textdisplay.jl")

public SpectralIntegral, integrate, integration_workspace
public earth_bindings, initialize_buffers, earth!, materials!, homogenize!,
       same_physical_state, layer_index, computation_type
public OBSERVABLE_RESOLUTION_REVISION
public has_uncertainty_type, numerical_magnitude
public internal_shunt_response, blueprint_dependencies
public InternalImpedanceFormulation, InsulationImpedanceFormulation,
       PipeImpedanceFormulation,
       EarthImpedanceFormulation, InsulationAdmittanceFormulation,
       SemiconAdmittanceFormulation,
       EarthAdmittanceFormulation, ShuntModelFormulation
public reduce_primitive_matrices
public layer_admittance
public ConsoleVerbosityLogger
public CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten, lineinput,
       earth_pairs

end # module Engine
