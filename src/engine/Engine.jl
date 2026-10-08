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
import ..LineCableModels

# Export public API
export LineParametersProblem, CableConstantsProblem,
       LineParameters, CableConstants, SeriesImpedance, ShuntAdmittance,
       RMSError, LineParametersBenchmark, compare,
       absolute_error, relative_error,
       Z, Y, R, X, L, G, B, C,
       series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance,
       frequencies, nconductors, nfrequencies, basis
export AbstractFormulation, LineParametersFormulation, CableConstantsFormulation,
       Formulation
export LineCableModelsCoaxial, LineCableModelsFEM,
       LineCableModelsFEMError, LineParametersWorkspace
export constitutive, formula_id, EarthPair
export verbosity
export InternalImpedance, InsulationImpedance, EarthImpedance, PipeImpedance
export InsulationAdmittance, SemiconAdmittance, EarthAdmittance
export ShuntModel, BoundarySolveError
export ModalAnalysis

export compute

# Module-specific dependencies
using LinearAlgebra: diag
import LinearAlgebra: norm
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, line_length, build, R, L, C,
                          resistance, inductance, capacitance
import ..LineCableModels: nominal
import ..LineCableModels: constitutive, formula, formula_id,
                          Expression, FormulaDefinition
import ..LineCableModels: parameterize
import ..LineCableModels: verbosity, VerbosityLogger
#! explicit-imports: off
import ..LineCableModels: description
#! explicit-imports: on
import ..Commons: AbstractProblemDefinition, AbstractFormulation,
                  AbstractProblemResult, AbstractCoreResult,
                  FormulationOptions, ComputationOptions,
                  ComputationDetails,
                  formulation_options, computation_options, computation_details, details,
                  compute, observe, observables,
                  observation_indices, observation_resolution,
                  uncertainty,
                  request_identity, request_indices

using ..Units
import ..Commons
using ..Commons: vacuum_permittivity, vacuum_permeability
using ..Commons: kron_reduce!, ReductionPlan, reduce_line_matrices!
import ..Commons: Functor
import ..Commons: initialize_buffers
using ..Materials
using ..Materials: TemperatureDependent
import ..Earth
using ..Earth: EarthMaterial, EarthModel, EquivalentHomogeneous
using ..DataModel: CableDesign, LineCableSystem, ncables, nphases
import ..DataModel
import ..TextDisplay
import ..LineCableModels: validate
import Logging
using Logging: with_logger
import SpecialFunctions
using QuadGK: alloc_segbuf, quadgk

include("interfaces.jl")
include("formulations.jl")
include("earthplan.jl")
include("specialfunctions.jl")

# Problem and coaxial formulation definitions
include("problems.jl")
include("options.jl")
include("integration.jl")
include("earthkernels.jl")

# Line-parameter results and their protocols
include("lineparameters/lineparameters.jl")
include("lineparameters/quantities.jl")
include("lineparameters/resolution.jl")
include("lineparameters/benchmark.jl")

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
include("earthreturn.jl")
include("impedance.jl")
include("admittance.jl")
include("lineparameters.jl")
include("cableconstants.jl")
include("observed_inputs.jl")
include("lineparameters/observations.jl")

# Line-parameter protocols and observation publication
include("lineparameters/base.jl")
include("textdisplay.jl")

# Submodule `ModalAnalysis`
include("modalanalysis/ModalAnalysis.jl")
using .ModalAnalysis: ModalAnalysis

public completion_details, completed_inputs, completed_formulation, retain_gridpoint
public selectdetails
public AbstractModalOperators
public SpectralIntegral, integrate
public AirVoltageSpectrum, earth_spectral_term, earth_spectral_value,
       earth_spectral_points!, earth_contour_angle, earth_direct,
       outgoing_root, bessel_i0m1, bessel_current_ratio, special_besselix
public earth!, materials!, homogenize!,
       same_physical_state, layer_index, computation_type
public has_uncertainty_type, numerical_magnitude
public resolution_available
public observation_assumptions
public domain, LineParamsDomain, PhaseDomain, ModalDomain, line_coordinates
public internal_shunt_response, blueprint_dependencies
public InternalImpedanceFormulation, InsulationImpedanceFormulation,
       PipeImpedanceFormulation,
       EarthImpedanceFormulation, InsulationAdmittanceFormulation,
       SemiconAdmittanceFormulation,
       EarthAdmittanceFormulation, ShuntModelFormulation
public layer_admittance
public CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten, lineinput,
       earth_pairs

end # module Engine
