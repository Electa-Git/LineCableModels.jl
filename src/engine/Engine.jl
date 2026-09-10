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
export AbstractFormulationBackend, AbstractFormulationOptions
export LineCableModelsCoaxial, LineCableModelsFEM, LineCableModelsFEMOptions,
       LineCableModelsFEMError, LineParametersWorkspace
export constitutive, formula_id, EarthPair
export verbosity
export InternalImpedance, InsulationImpedance, EarthImpedance, PipeImpedance
export InsulationAdmittance, SemiconAdmittance, EarthAdmittance

export compute

# Module-specific dependencies
using LinearAlgebra: svd, svd!, eigvals, eigvals!, Diagonal, I, checksquare, cond, diag, ldiv!, lu, lu!,
                     mul!, qr, ColumnNorm
import LinearAlgebra: norm
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, build, R, L, C,
                          resistance, inductance, capacitance
import ..LineCableModels: nominal
import ..LineCableModels: constitutive, formula, formula_id,
                          FormulaMethod
import ..LineCableModels: parameterize
#! explicit-imports: off
import ..LineCableModels: description
#! explicit-imports: on
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  AbstractProblemResult, AbstractCoreResult,
                  FormulationOptions, ComputationOptions,
                  ComputationDetails,
                  formulation_options, computation_options, computation_details, details,
                  compute, observe, observables,
                  observation_request, observation_indices,
                  publication_table

using ..Units
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
using DoubleExponentialFormulas: QuadDE

include("interfaces.jl")
include("formulations.jl")
include("earthkernels.jl")

# Problem and coaxial formulation definitions
include("problems.jl")
include("options.jl")
include("integration.jl")
include("spectralsampling.jl")
include("earthbounds.jl")
include("compleximages.jl")

# Line-parameter results and their protocols
include("lineparameters/lineparameters.jl")
include("lineparameters/quantities.jl")
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

public hooks, SpectralIntegral, integrate
public has_uncertainty_type, spectral_magnitude
public reduce_primitive_matrices, potential_to_admittance
public layer_admittance
public ConsoleVerbosityLogger
public CableBlueprint, BlueprintConductor, BlueprintDielectric, flatten, lineinput,
       earth_pairs

end # module Engine
