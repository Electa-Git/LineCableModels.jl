module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, validate, description, maxfill, set_backend!
export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions
export formulation_options, computation_options
export compute, observables, primitives, preprocess
export basis, domain, frequencies, nconductors, nfrequencies, ncables, nphases
export Z, Y, R, X, L, G, B, C
export series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance
# High-level modeling grammar:
export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace
export has_uncertainty, nominal, standard_uncertainty
export @gridspace, @relax
export Combinatorial, LinearError, MonteCarlo, ParametricProblem
export ParametricResult, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity, RLCG
export result, statistics, samples, histograms, uncertain_value
export Material, MaterialsLibrary, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export make_stranded, make_screened, WireEstimate

# Materialized results, reusable designs, and presentation:
export CableConstants, CableConstantsProblem, LineParameters, CablesLibrary, preview

# Engine:
export Formulation, AnalyticalFormulation, SeriesImpedance, ShuntAdmittance, kronify,
       LineParameters, PhaseDomain, ModalDomain
export Fortescue
export LineParametersTrace

# Import/Export:
export export_data, import_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions

include("docstrings.jl")
include("interfaces.jl")
include("retired.jl")

# Package-local shared calculation grammar.
include("Grammar.jl")
using .Grammar:
                AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                AbstractParametricResult, AbstractUncertaintyResult,
                FormulationOptions, ComputationOptions,
                formulation_options, computation_options,
                compute, observables, primitives, preprocess,
                nominal, standard_uncertainty

# Submodule `UnitHandler`
include("unithandler/UnitHandler.jl")

# Submodule `Validation`
include("validation/Validation.jl")
using .Validation: validate

# Submodule `PlotBuilder`
include("plotbuilder/PlotBuilder.jl")
using .PlotBuilder.BackendHandler: set_backend!
using .PlotBuilder: UIPlot, export_svg
export UIPlot, export_svg

# Submodule `Materials`
include("materials/Materials.jl")
import .Materials
using .Materials: Material, MaterialsLibrary

# Submodule `EarthProps`
include("earthprops/EarthProps.jl")
import .EarthProps

# Submodule `DataModel`
include("datamodel/DataModel.jl")
using .DataModel: CableConstants, CablesLibrary, preview

# Submodule `Engine`
include("engine/Engine.jl")
using .Engine: LineParameters, SeriesImpedance,
               ShuntAdmittance, kronify, Formulation, CableConstantsProblem,
               AnalyticalFormulation, LineParametersTrace
using .Engine.Transforms: Fortescue

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")
using .ParametricBuilder:
                          Grid, AbsoluteError, DeterministicGrid, RelativeGrid,
                          AbsoluteGrid, AbstractGrid, AbstractUncertainGrid,
                          UncertainValue, Gridspace, has_uncertainty,
                          @gridspace, @relax,
                          Combinatorial, ParametricProblem, ParametricResult,
                          result,
                          Conductor, Insulator, CableBuilder,
                          at, trifoil, hflat, vflat, Earth, SystemBuilder,
                          WireEstimate, make_stranded, make_screened

# Submodule `UQ`
include("uq/UQ.jl")
using .UQ:
           LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult,
           SampleSummary, HistogramDensity, RLCG,
           statistics, samples, histograms, uncertain_value

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, import_data, load!, save

end
