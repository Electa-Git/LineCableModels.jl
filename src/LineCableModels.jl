module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, set_verbosity!, set_backend!
export basis, domain, frequencies, nconductors, nfrequencies
export Z, Y, R, X, L, G, B, C
export series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance
export ParametricSweep, cases, results, ncases

# Materials:
export Material, MaterialsLibrary

# Data model (design + system):
export Thickness, Diameter, WireArray, Strip, Tubular, Semicon, Insulator, Sector,
       SectorParams, SectorInsulator
export ConductorGroup, InsulatorGroup
export CableComponent, CableDesign, CableConstants, NominalData
export CablesLibrary
export CablePosition, LineCableSystem
export trifoil_formation, flat_formation, preview, equivalent, MaxFill

# Earth properties:
export EarthModel

# Engine:
export LineParametersProblem,
       FormulationSet, compute!, SeriesImpedance, ShuntAdmittance, kronify,
       LineParameters, PhaseDomain, ModalDomain

# Parametric builder:
# export make_stranded, make_screened
# export conductor, insulator

# Import/Export:
export export_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions

# Submodule `Commons`
include("commons/Commons.jl")
using .Commons: IMPORTS, EXPORTS, add!, PhaseDomain, ModalDomain, domain,
                basis, Z, Y, R, X, L, G, B, C,
                series_impedance, shunt_admittance,
                resistance, reactance, inductance,
                conductance, susceptance, capacitance,
                frequencies, nconductors, nfrequencies
# Submodule `UncertainBessels`
include("uncertainbessels/UncertainBessels.jl")

# Submodule `Utils`
include("utils/Utils.jl")
using .Utils: set_verbosity!

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
using .Materials: Material, MaterialsLibrary

# Submodule `EarthProps`
include("earthprops/EarthProps.jl")
using .EarthProps: EarthModel

# Submodule `DataModel`
include("datamodel/DataModel.jl")
using .DataModel: Thickness, Diameter, CircStrands, RectStrands, Strip, Tubular, Semicon,
                  Insulator, ConductorGroup, InsulatorGroup, CableComponent, CableDesign,
                  CableConstants,
                  NominalData,
                  CablesLibrary, CablePosition, LineCableSystem, trifoil_formation,
                  flat_formation,
                  preview, equivalent, MaxFill, Sector, SectorParams,
                  SectorInsulator

# Submodule `Engine`
include("engine/Engine.jl")
using .Engine: LineParametersProblem, compute!, LineParameters, SeriesImpedance,
               ShuntAdmittance, kronify, FormulationSet

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")
using .ParametricBuilder: ParametricSweep, cases, results, ncases

# Submodule `UQ`
include("uq/UQ.jl")
using .UQ: SampleSummary, RLCG, HistogramPDF, CableConstantsMC, LineParametersMC,
           sample, trial, mc, statistics, has_samples, samples,
           has_distributions, distribution, surrogate, ntrials, confidence,
           mean, std, quantile
export SampleSummary, RLCG, HistogramPDF, CableConstantsMC, LineParametersMC,
       sample, trial, mc, statistics, has_samples, samples,
       has_distributions, distribution, surrogate, ntrials, confidence,
       mean, std, quantile

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, load!, save

# Aliases for backward compatibility
const WireArray = CircStrands  # alias for now
export WireArray  # export aliases

end
