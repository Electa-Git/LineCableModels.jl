"""
    LineCableModels.DataModel

Define materialised cable geometry and line arrangements.

# Overview

- Define conductive and dielectric cable parts.
- Assemble parts into [`CableDesign`](@ref) values.
- Place cables in a [`LineCableSystem`](@ref).
- Calculate base-state resistance, inductance, capacitance, and conductance.
- Describe cable and system previews for PlotBuilder.
- Store cable designs in [`CablesLibrary`](@ref).

# Dependencies

$(IMPORTS)

"""
module DataModel

# Export public API
export CircStrands, RectStrands, Strip, Tubular  # Conductor types
export Semicon, Insulator  # Insulator types
export ConductorGroup, InsulatorGroup  # Group types
export CableComponent, CableDesign, CableConstants  # Cable design types
export CablePosition, LineCableSystem  # System types
export CablesLibrary, NominalData  # Support types
export trifoil_formation, flat_formation, outer_radius, maxfill  # Geometry
export AbstractPrimitive, AbstractShape
export Disk, Rectangle, Ellipse, Sector, Annulus, Shell, Polygon, Pose2
export EmptyBoundary, DiskShape, RectangleShape, EllipseShape, SectorShape
export AnnulusShape, PolygonShape, PlacedShape
export resolve, boundary, area, centroid, support, r_in, r_ex, thickness
export ncables, nphases
export preview, equivalent

# Module-specific dependencies
#! explicit-imports: off
# IMPORTS is expanded in the module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES, FUNCTIONNAME
import ..PlotBuilder
import ..Units
import ..LineCableModels: add!, validate, maxfill, nominal
import ..LineCableModels: basis, R, L, C, resistance, inductance, capacitance
import ..Grammar: AbstractCoreResult, observe, observables
using ..Materials: Material
import ..Validation
using ..Validation: IntegerField, Positive, Finite, Nonnegative, OneOf, Less
using Colors: HSL, RGB, RGBA, alpha, blue, green, red
import GeometryBasics
using GeometryBasics: Point2f
using Printf: @sprintf
using Statistics: mean

# Abstract types and interfaces
include("interfaces.jl")
include("types.jl")
include("geometry/primitives.jl")
include("geometry/shapes.jl")
include("geometry/pose.jl")
include("geometry/resolve.jl")

# Submodule `BaseParams`
include("baseparams/BaseParams.jl")
using .BaseParams

# Constructors
# Conductors
include("packing.jl")
include("circstrands.jl")
include("rectstrands.jl")
include("strip.jl")
include("tubular.jl")
include("conductorgroup.jl")

# Insulators
include("insulator.jl")
include("semicon.jl")
include("insulatorgroup.jl")

# Groups
include("nominaldata.jl")
include("cablecomponent/cablecomponent.jl")
include("cabledesign/cabledesign.jl")

# Library
include("cableslibrary/cableslibrary.jl")
include("linecablesystem/linecablesystem.jl")

# Geometry and language protocols
include("geometry.jl")
include("base.jl")
include("preview/definitions.jl")
include("preview/materials.jl")
include("preview/cable.jl")
include("preview/cables.jl")
include("preview/system.jl")
include("preview/materialscale.jl")

public base_parameters, compute_cable_constants
public preview_shapes, preview_materials
public PreviewPolygon, PreviewReferenceLine, PreviewPayload

end # module DataModel
