"""
    LineCableModels.Engine

Calculate cable constants and frequency-dependent line-parameter matrices from
materialised cable systems.

# Overview

- Define materialised problems, formulations, and core results.
- Calculate conductor, insulation, and earth-return impedance and admittance.
- Assemble phase-domain series-impedance and shunt-admittance matrices.
- Apply bundle reduction, Kron elimination, and ideal transposition.
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
export AbstractFormulationBackend, AbstractFormulationOptions
export LineCableModelsCoaxial, LineCableModelsFEM, LineCableModelsFEMOptions,
       LineCableModelsFEMError, LineParametersWorkspace
export constitutive, formula_id, EarthPair
export verbosity
export InternalImpedance, InsulationImpedance, EarthImpedance
export InsulationAdmittance, SemiconAdmittance, EarthAdmittance

export compute, plot

# Module-specific dependencies
using LinearAlgebra: I, checksquare, cond, diag, ldiv!, lu, lu!, mul!, norm
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, build, R, L, C,
                          resistance, inductance, capacitance
import ..LineCableModels: nominal, uncertainty
import ..LineCableModels: constitutive, formula, formula_id, FormulaSpec
#! explicit-imports: off
import ..LineCableModels: description
#! explicit-imports: on
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  AbstractProblemResult, AbstractCoreResult,
                  FormulationOptions, ComputationOptions,
                  ComputationDetails,
                  formulation_options, computation_options, computation_details, details,
                  compute, observe, @observe, observables, validate_observables,
                  unit_targets,
                  observation_request, observation_indices,
                  request_identity, request_quantity, request_indices,
                  computation_owner, publication_table

using ..Units
using ..PlotBuilder
using ..Materials
import ..EarthProps
using ..EarthProps: EarthMaterial, EarthModel, EHEM
using ..DataModel: CableDesign, LineCableSystem, ncables, nphases
import ..DataModel
import ..TextDisplay
import ..LineCableModels: validate
import ..Validation
import Logging
using Logging: AbstractLogger, ConsoleLogger, with_logger
import SpecialFunctions
using QuadGK: alloc_segbuf, quadgk

include("interfaces.jl")
include("formulations.jl")
include("earthkernels.jl")

# Problem and coaxial formulation definitions
include("problems.jl")
include("options.jl")

# Line-parameter results and their protocols
include("lineparameters/lineparameters.jl")
include("lineparameters/quantities.jl")
include("lineparameters/benchmark.jl")
include("matrixops.jl")

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

# Submodule `SemiconAdmittance`
include("semiconadmittance/SemiconAdmittance.jl")
using .SemiconAdmittance: SemiconAdmittance

# Submodule `EarthAdmittance`
include("earthadmittance/EarthAdmittance.jl")
using .EarthAdmittance: EarthAdmittance

# Native workspace and numerical action
include("input.jl")
include("logging.jl")
include("earthreturn.jl")
include("impedance.jl")
include("admittance.jl")
include("compute.jl")
include("reduction.jl")

# Line-parameter protocols and renderer-independent plot definitions
include("lineparameters/base.jl")
include("lineparameters/publication.jl")
include("lineparameters/plot.jl")
include("lineparameters/plotdefinition.jl")
include("lineparameters/comparisonplot.jl")
include("textdisplay.jl")

public has_uncertainty_type
public reduce_primitive_matrices, potential_to_admittance
public ConsoleVerbosityLogger

end # module Engine
