"""
    LineCableModels

Calculate electrical parameters for overhead and underground cable systems.

The public API constructs materialized or finite parametric cable models,
evaluates cable constants and line-parameter matrices with selected numerical formulations. It propagates declared uncertainty and provides plotting methods when Makie is loaded.
"""
module LineCableModels

## Public API
# -------------------------------------------------------------------------
# Core generics:
export add!, build, homogenize, validate, description, constitutive
export formula, formula_id
export AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult
export AbstractCoreResult, AbstractResultSpace
export AbstractParametricResult, AbstractUncertaintyResult
export FormulationOptions, ComputationOptions, ComputationDetails
export formulation_options, computation_options, computation_details, details
export compute, observe, @observe, observables
export ObservedResult, kron_reduce
export quantity, native_unit, display_unit, scale_factor, label, symbol
export basis, line_length, domain, frequencies, nconductors, nfrequencies, ncables, nphases
export Z, Y, R, X, L, G, B, C
export series_impedance, shunt_admittance,
       resistance, reactance, inductance,
       conductance, susceptance, capacitance
# High-level modelling grammar:
export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace
export has_uncertainty, nominal, uncertainty
export @gridspace
export Combinatorial, LinearError, MonteCarlo, ParametricProblem
export ParametricResult, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity
export statistics, samples, histograms, uncertain
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution
export report, TableReportDefinition, XLSXReportDefinition, ReportArtifact
export AbstractMaterial, Material, RadialDielectric, MaterialsLibrary
export AbstractShape, AbstractPrimitive
export Disk, Rectangle, Ellipse, Sector, Annulus, Polygon, Shell
export Pose2
export EmptyBoundary, resolve, boundary, area, perimeter, centroid, support, tessellate
export r_in, r_ex, thickness, outer_radius
export AbstractCablePart, Region, Stack
export Group, Assembly
export Enclosure
export Ring, Polar, Fill, Lattice, placements
export capacity, FillFactor
export LayRatio, Pitch, LayAngle, Helix, pitch, angle, overlength
export at, trefoil, hflat, vflat, layer, homogeneous
export AbstractEarthModel, AbstractEarthLayer, AbstractEarthMaterial, EarthMaterial,
       EarthLayer, EarthModel
export terminal, core, stranded, milliken, rope, cores, tape, insulation, screen, sheath
export armor, bedding, jacket, filler, pipe, duct
export solid, shell, wires, layers, assembly
export @cable, @system, @earth, @terminal, @assembly, @pipe, @duct
export @at, @hflat, @vflat, @trefoil
export @distribute
export estimate_stranding, estimate_screen, WireEstimate

public Gridpoint

# Materialized results, reusable designs, and presentation:
export CableDesign, LineCableSystem, DatasheetInfo, datasheet
export CableGeometry, PlacedRegion
export CableConstants, CableConstantsProblem, CableConstantsFormulation,
       LineParametersProblem, LineParameters, CablesLibrary
export preview, show_material_scale

# Engine:
export Formulation, LineParametersFormulation, CableConstantsFormulation,
       LineCableModelsCoaxial,
       LineCableModelsFEM, LineCableModelsFEMError, BoundarySolveError,
       SeriesImpedance, ShuntAdmittance,
       LineParameters, PhaseDomain, ModalDomain
export ModalAnalysisProblem, ModalAnalysisFormulation,
       LineCableModelsModal, ModalOperators, operators, Tv, Ti, gamma, alpha, beta, velocity, Zc, Yc,
       PropagationParameters, H, transform
export ModalAnalysis

# Import/Export:
export export_data, import_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions
using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random
import Logging

include("docstrings.jl")
include("interfaces.jl")

public FormulaDefinition, FormulaMethod

# Submodule `Units`
include("units/Units.jl")
using .Units: quantity, native_unit, display_unit, scale_factor, label, symbol

# Package-local shared calculation grammar.
include("commons/Commons.jl")
using .Commons:
                AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                AbstractCoreResult, AbstractResultSpace,
                AbstractParametricResult, AbstractUncertaintyResult,
                FormulationOptions, ComputationOptions, ComputationDetails,
                formulation_options, computation_options, computation_details, details,
                observe, @observe, observables, ObservedResult, kron_reduce
import .Commons: compute, validate
using .Commons: FormulaDefinition, FormulaMethod
include("logging.jl")
include("formulas.jl")

# Bounded text formatting consumes the completed shared declaration grammar.
include("textdisplay/TextDisplay.jl")
import .Commons: nominal, uncertainty

# Root-owned finite parameter grammar. These primitives are loaded before the
# domain modules so every constructor enters the same scalar-or-Gridspace path
# without a late-loaded bridge through ParametricBuilder.
include("grid.jl")
include("gridspace.jl")

public parameterize, materialize, sample_uncertainty

# Thin native plotting handles and optional-extension entry points.
include("plotbuilder/PlotBuilder.jl")
using .PlotBuilder:
                    UIPlot, plot, preview, show_material_scale, export_svg,
                    figurelegend!, panellegend!, figuretitle!, paneltitle!,
                    figurecolorbars!, axisscale!, resetview!, addwidget!, removewidget!,
                    plotwindow, materialcolors, materialscale!
export UIPlot, export_svg, figurelegend!, panellegend!, figuretitle!, paneltitle!
export figurecolorbars!, axisscale!, resetview!, addwidget!, removewidget!
export materialcolors, materialscale!
public PlotBuilder, plot, plotwindow

# Submodule `Materials`
include("materials/Materials.jl")
import .Materials
using .Materials: AbstractMaterial, Material, RadialDielectric, MaterialsLibrary

# Submodule `Earth`
include("earth/Earth.jl")
import .Earth
public Earth
using .Earth: AbstractEarthModel, AbstractEarthLayer, AbstractEarthMaterial,
              EarthMaterial, EarthLayer, EarthModel, layer, homogeneous

# Submodule `DataModel`
include("datamodel/DataModel.jl")
using .DataModel: CableDesign, CableGeometry, PlacedRegion,
                  LineCableSystem, DatasheetInfo, datasheet,
                  CablesLibrary, ncables, nphases,
                  AbstractShape, AbstractPrimitive,
                  Disk, Rectangle, Ellipse, Sector, Annulus, Polygon, Shell,
                  Pose2,
                  EmptyBoundary, resolve, boundary, area, perimeter, centroid, support,
                  tessellate,
                  r_in, r_ex, thickness, outer_radius,
                  AbstractCablePart, Region, Stack,
                  Group, Assembly, Enclosure
using .DataModel: Ring, Polar, Fill, Lattice, capacity, placements,
                  FillFactor,
                  LayRatio, Pitch, LayAngle, Helix, pitch, angle, overlength

# Submodule `Engine`
include("engine/Engine.jl")
using .Engine: LineParameters, LineParametersProblem, CableConstants,
               CableConstantsProblem, CableConstantsFormulation, SeriesImpedance,
               ShuntAdmittance, Formulation,
               LineParametersFormulation, LineCableModelsCoaxial,
               LineCableModelsFEM,
               LineCableModelsFEMError, BoundarySolveError,
               domain, frequencies, nconductors, nfrequencies,
               Z, Y, X, G, B, series_impedance, shunt_admittance,
               reactance, conductance, susceptance,
               LineParamsDomain, PhaseDomain, ModalDomain

public LineParamsDomain

# Submodule `ModalAnalysis`
include("modalanalysis/ModalAnalysis.jl")
using .ModalAnalysis: ModalAnalysisProblem, ModalAnalysisFormulation,
                   LineCableModelsModal, ModalOperators, operators,
                   Tv, Ti, gamma, alpha, beta, velocity, Zc, Yc, PropagationParameters, H, transform

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")
using .ParametricBuilder:
                          @gridspace,
                          Combinatorial, ParametricProblem, ParametricResult,
                          terminal, core, stranded, milliken, rope, cores, tape,
                          insulation, screen, sheath, armor, bedding, jacket,
                          filler, pipe, duct, solid, shell, wires, layers,
                          assembly,
                          at, trefoil, hflat, vflat,
                          WireEstimate, estimate_stranding, estimate_screen
using .ParametricBuilder: @cable, @system, @earth, @terminal, @assembly, @pipe,
                          @duct, @at, @hflat, @vflat, @trefoil, @distribute

# Submodule `UQ`
include("uq/UQ.jl")
using .UQ:
           LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult,
           SampleSummary, HistogramDensity,
           statistics, samples, histograms, uncertain,
           root_seed, point_seed, trial_count,
           confidence, cdf_tolerance, sampling_distribution

# Completed-result measurement projections.
include("performance.jl")
public benchmark

# Submodule `ReportBuilder`
include("reportbuilder/ReportBuilder.jl")
using .ReportBuilder:
                      report, TableReportDefinition, XLSXReportDefinition, ReportArtifact

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, import_data, load!, save

# External-tool integration. Native execution is deferred until compute.
include("pscad/PSCAD.jl")
export PSCAD

end
