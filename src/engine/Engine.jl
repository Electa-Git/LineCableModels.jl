"""
    LineCableModels.Engine

Calculate cable constants and frequency-dependent line-parameter matrices from
materialised cable systems.

# Overview

- Define materialised problems, formulations, and primitive results.
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
export CableConstantsProblem, LineParametersProblem,
       LineParameters, SeriesImpedance, ShuntAdmittance,
       RMSError, LineParametersBenchmark, compare,
       absolute_error, relative_error,
       Z, Y, R, X, L, G, B, C,
       series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance,
       frequencies, nconductors, nfrequencies, basis,
       kronify
export AbstractFormulation, AnalyticalFormulation, Formulation
export EarthProperties, CPEarth, LineParametersTrace, verbosity
export InternalImpedance, InsulationImpedance, EarthImpedance
export InsulationAdmittance, EarthAdmittance, EHEM, Transforms

export compute, plot

# Module-specific dependencies
using LinearAlgebra
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES,
                           METHODLIST, FUNCTIONNAME
import ..LineCableModels: retired_fem_sector, _RETIRED_FEM
import ..LineCableModels: basis, R, L, C, resistance, inductance, capacitance
import ..LineCableModels: nominal, standard_uncertainty
import ..Grammar: AbstractProblemDefinition, AbstractFormulation,
                  AbstractProblemResult, FormulationOptions, ComputationOptions,
                  formulation_options, computation_options, compute, observe, observables

using ..Units
using ..PlotBuilder
using ..Materials
using ..EarthProps: EarthModel
using ..DataModel: CableDesign, CableConstants, LineCableSystem, ncables, nphases
import ..DataModel
import ..LineCableModels: validate
import ..Validation
using Logging
using SpecialFunctions
using QuadGK: quadgk

const FEM = _RETIRED_FEM

include("interfaces.jl")
include("formulations.jl")
include("earthkernels.jl")

include("earthproperties/EarthProperties.jl")
using .EarthProperties: CPEarth

# Problem and analytical formulation definitions
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

# Immutable solver input and numerical action
include("input.jl")
include("lineparameters/trace.jl")
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

end # module Engine
