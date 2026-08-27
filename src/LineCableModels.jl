"""
    LineCableModels

Calculate electrical parameters for overhead and underground cable systems.

The public API constructs materialised or finite parametric cable models,
selects numerical formulations, evaluates cable constants and line-parameter
matrices, propagates declared uncertainty, and describes plots for optional
Makie renderers.
"""
module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, validate, description, maxfill, set_backend!
export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractCoreResult, AbstractResultSpace
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables, project
export quantity, native_unit, display_unit, scale_factor, label, symbol
export basis, domain, frequencies, nconductors, nfrequencies, ncables, nphases
export Z, Y, R, X, L, G, B, C
export series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance
# High-level modelling grammar:
export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace
export has_uncertainty, nominal, standard_uncertainty
export @gridspace, @relax
export Combinatorial, LinearError, MonteCarlo, ParametricProblem
export ParametricResult, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity
export result, statistics, samples, histograms, uncertain_value
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution
export report, TableReportDefinition, XLSXReportDefinition, ReportArtifact
export Material, MaterialsLibrary, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export make_stranded, make_screened, WireEstimate

# Materialised results, reusable designs, and presentation:
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

# Submodule `Units`
include("units/Units.jl")
using .Units: quantity, native_unit, display_unit, scale_factor, label, symbol

# Package-local shared calculation grammar.
include("grammar/Grammar.jl")
using .Grammar:
                AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                AbstractCoreResult, AbstractResultSpace,
                AbstractParametricResult, AbstractUncertaintyResult,
                FormulationOptions, ComputationOptions, ComputationDetails,
                formulation_options, computation_options, computation_details, details,
                compute, observe, @observe, observables,
                nominal, standard_uncertainty

# Submodule `Validation`
include("validation/Validation.jl")
using .Validation: validate

# Submodule `PlotBuilder`
include("plotbuilder/PlotBuilder.jl")
using .PlotBuilder: UIPlot, export_svg, set_backend!
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
using .DataModel: CableConstants, CablesLibrary, preview, ncables, nphases

# Submodule `Engine`
include("engine/Engine.jl")
using .Engine: LineParameters, SeriesImpedance,
               ShuntAdmittance, kronify, Formulation, CableConstantsProblem,
               AnalyticalFormulation, LineParametersTrace,
               description, domain, frequencies, nconductors, nfrequencies,
               Z, Y, X, G, B, series_impedance, shunt_admittance,
               reactance, conductance, susceptance,
               LineParamsDomain, PhaseDomain, ModalDomain
using .Engine.Transforms: Fortescue

public LineParamsDomain

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")
using .ParametricBuilder:
                          Grid, AbsoluteError, DeterministicGrid, RelativeGrid,
                          AbsoluteGrid, AbstractGrid, AbstractUncertainGrid,
                          UncertainValue, Gridspace, has_uncertainty,
                          @gridspace, @relax,
                          Combinatorial, ParametricProblem, ParametricResult,
                          result, project,
                          Conductor, Insulator, CableBuilder,
                          at, trifoil, hflat, vflat, Earth, SystemBuilder,
                          WireEstimate, make_stranded, make_screened

# Submodule `UQ`
include("uq/UQ.jl")
using .UQ:
           LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult,
           SampleSummary, HistogramDensity,
           statistics, samples, histograms, uncertain_value,
           root_seed, point_seed, trial_count,
           confidence, cdf_tolerance, sampling_distribution

# Submodule `ReportBuilder`
include("reportbuilder/ReportBuilder.jl")
using .ReportBuilder:
                      report, TableReportDefinition, XLSXReportDefinition, ReportArtifact

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, import_data, load!, save

end
