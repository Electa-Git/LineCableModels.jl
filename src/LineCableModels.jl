module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, validate, description, maxfill, set_backend!
export basis, domain, frequencies, nconductors, nfrequencies, ncables, nphases
export Z, Y, R, X, L, G, B, C
export series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance
# High-level modeling grammar:
export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export AbstractSpec, Gridspace, Configuration, configurations, materialize
export has_uncertainty, configuration_manifest, nominal, standard_uncertainty
export Material, MaterialsLibrary, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export make_stranded, make_screened, WireEstimate

# Materialized results, reusable designs, and presentation:
export CableConstants, CableConstantsProblem, LineParameters, CablesLibrary, preview

# Engine:
export Formulation, compute!, SeriesImpedance, ShuntAdmittance, kronify,
       LineParameters, PhaseDomain, ModalDomain
export Fortescue
export EMTTrace
export FullParametric, MonteCarlo, FullParametricResult, MonteCarloResult
export CalculationManifest, ConfigurationFailure, SampleSummary, HistogramPDF, RLCG
export result, statistics, samples, histograms, uncertain_value, manifest

# Import/Export:
export export_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions

include("docstrings.jl")
include("interfaces.jl")
include("retired.jl")

# Submodule `UnitHandler`
include("unithandler/UnitHandler.jl")

# Submodule `Validation`
include("validation/Validation.jl")

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
using .Engine: compute!, LineParameters, SeriesImpedance,
               ShuntAdmittance, kronify, Formulation, CableConstantsProblem,
               EMTTrace
using .Engine.Transforms: Fortescue

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")
using .ParametricBuilder:
                          Grid, AbsoluteError, DeterministicGrid, RelativeGrid,
                          AbsoluteGrid,
                          AbstractGrid, AbstractUncertainGrid, UncertainValue,
                          AbstractSpec, Gridspace, Configuration, configurations,
                          materialize,
                          has_uncertainty, configuration_manifest,
                          Conductor, Insulator, CableBuilder,
                          at, trifoil, hflat, vflat, Earth, SystemBuilder,
                          WireEstimate, make_stranded, make_screened

# Unified Gridspace computation policies and typed result grammar.
include("computation/Computation.jl")
using .Computation:
                    FullParametric, MonteCarlo, FullParametricResult, MonteCarloResult,
                    CalculationManifest, ConfigurationFailure, SampleSummary, HistogramPDF,
                    RLCG,
                    result, statistics, samples, histograms, uncertain_value, manifest

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, load!, save

end
