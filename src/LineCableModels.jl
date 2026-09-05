"""
    LineCableModels

Calculate electrical parameters for overhead and underground cable systems.

The public API constructs materialised or finite parametric cable models,
selects numerical formulations, evaluates cable constants and line-parameter
matrices, propagates declared uncertainty, and adds a small high-level plotting
surface when Makie is loaded.
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
export has_uncertainty, nominal, uncertainty
export @gridspace
export Combinatorial, LinearError, MonteCarlo, ParametricProblem
export ParametricResult, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity
export statistics, samples, histograms, uncertain
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution
export report, TableReportDefinition, XLSXReportDefinition, ReportArtifact
export AbstractMaterial, Material, MaterialsLibrary, Conductor, Insulator, Semiconductor
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
export AbstractEarthModel, EarthLayer, EarthModel
export terminal, core, stranded, milliken, rope, cores, tape, insulation, screen, sheath
export armor, bedding, jacket, filler, pipe, duct
export solid, shell, wires, layers, assembly
export @cable, @system, @earth, @terminal, @assembly, @pipe, @duct
export @at, @hflat, @vflat, @trefoil
export @distribute
export make_stranded, make_screened, WireEstimate

public Gridpoint

# Materialised results, reusable designs, and presentation:
export CableDesign, LineCableSystem, DatasheetInfo, catalogue
export CableGeometry, PlacedRegion
export CableConstants, CableConstantsProblem, CableConstantsFormulation,
       LineParametersProblem, LineParameters, CablesLibrary
export preview, show_material_scale

# Engine:
export Formulation, LineParametersFormulation, CableConstantsFormulation,
       LineCableModelsCoaxial,
       LineCableModelsFEM, LineCableModelsFEMOptions, LineCableModelsFEMError,
       SeriesImpedance, ShuntAdmittance, kronify,
       LineParameters, PhaseDomain, ModalDomain
export ModalTransformationProblem, ModalTransformationFormulation,
       LineCableModelsModal, ModalOperators, operators
export Transforms

# Import/Export:
export export_data, import_data, save, load!
# -------------------------------------------------------------------------

import DocStringExtensions: DocStringExtensions
using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random

include("docstrings.jl")
include("interfaces.jl")
include("formulas.jl")

public FormulaMethod

# Submodule `Units`
include("units/Units.jl")
using .Units: quantity, native_unit, display_unit, scale_factor, label, symbol

# Package-local bounded text formatting. Domain modules extend Base display
# methods and call this owner through qualified operations.
include("textdisplay/TextDisplay.jl")

# Package-local shared calculation grammar.
include("grammar/Grammar.jl")
using .Grammar:
                AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                AbstractCoreResult, AbstractResultSpace,
                AbstractParametricResult, AbstractUncertaintyResult,
                FormulationOptions, ComputationOptions, ComputationDetails,
                formulation_options, computation_options, computation_details, details,
                compute, observe, @observe, observables
import .Grammar: nominal, uncertainty

# Submodule `InputValidation`
include("inputvalidation/InputValidation.jl")
using .InputValidation: validate

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
                    plotwindow, materialcolors, materialscale!
export UIPlot, export_svg, figurelegend!, panellegend!, figuretitle!, paneltitle!
export materialcolors, materialscale!
public PlotBuilder, plot, plotwindow

# Submodule `Materials`
include("materials/Materials.jl")
import .Materials
using .Materials: AbstractMaterial, Material, MaterialsLibrary

# Submodule `EarthProps`
include("earthprops/EarthProps.jl")
import .EarthProps
using .EarthProps: AbstractEarthModel, EarthLayer, EarthModel, layer, homogeneous

# Submodule `DataModel`
include("datamodel/DataModel.jl")
using .DataModel: CableDesign, CableGeometry, PlacedRegion,
                  LineCableSystem, DatasheetInfo, catalogue,
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
               ShuntAdmittance, kronify, Formulation,
               LineParametersFormulation, LineCableModelsCoaxial,
               LineCableModelsFEM, LineCableModelsFEMOptions,
               LineCableModelsFEMError,
               domain, frequencies, nconductors, nfrequencies,
               Z, Y, X, G, B, series_impedance, shunt_admittance,
               reactance, conductance, susceptance,
               LineParamsDomain, PhaseDomain, ModalDomain

public LineParamsDomain

# Submodule `Transforms`
include("transforms/Transforms.jl")
using .Transforms: ModalTransformationProblem, ModalTransformationFormulation,
                   LineCableModelsModal, ModalOperators, operators

# Submodule `ParametricBuilder`
include("parametricbuilder/ParametricBuilder.jl")
using .ParametricBuilder:
                          @gridspace,
                          Combinatorial, ParametricProblem, ParametricResult,
                          Conductor, Insulator,
                          terminal, core, stranded, milliken, rope, cores, tape,
                          insulation, screen, sheath, armor, bedding, jacket,
                          filler, pipe, duct, solid, shell, wires, layers,
                          assembly,
                          at, trefoil, hflat, vflat,
                          WireEstimate, make_stranded, make_screened
using .ParametricBuilder: Semiconductor
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

# Submodule `ReportBuilder`
include("reportbuilder/ReportBuilder.jl")
using .ReportBuilder:
                      report, TableReportDefinition, XLSXReportDefinition, ReportArtifact

# Submodule `ImportExport`
include("importexport/ImportExport.jl")
using .ImportExport: export_data, import_data, load!, save

end
