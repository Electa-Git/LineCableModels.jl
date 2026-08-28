"""
    LineCableModels.DataModel

Define the physical cable object model and completed line arrangements.

# Overview

- Declare primitives, material-bearing regions, and resolved geometry.
- Compose regions through [`Stack`](@ref), [`Group`](@ref),
  [`Assembly`](@ref), and [`Enclosure`](@ref).
- Build one physical declaration as a completed [`CableDesign`](@ref).
- Place completed designs and resolve connections as a [`LineCableSystem`](@ref).
- Describe cable and system previews for PlotBuilder.
- Store cable designs in [`CablesLibrary`](@ref).

# Dependencies

$(IMPORTS)

"""
module DataModel

# Export public API
export CableDesign, CableGeometry, PlacedRegion, LineCableSystem, CableConstants
export build
export CablesLibrary, NominalData
export trefoil_formation, flat_formation, outer_radius
export AbstractPrimitive, AbstractShape
export AbstractCablePart, Region, Stack
export Group, Assembly
export Enclosure
export Disk, Rectangle, Ellipse, Sector, Annulus, Shell, Polygon, Pose2
export EmptyBoundary, DiskShape, RectangleShape, EllipseShape, SectorShape
export AnnulusShape, PolygonShape, PlacedShape
export resolve, boundary, area, centroid, support, r_in, r_ex, thickness
export Ring, Polar, Lattice, DiameterFactor, placements
export LayRatio, Pitch, LayAngle, Helix, pitch, angle, overlength
export ncables, nphases
export preview

# Module-specific dependencies
#! explicit-imports: off
# IMPORTS is expanded in the module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES, FUNCTIONNAME
import ..PlotBuilder
import ..Units
import ..LineCableModels: add!, build, validate, nominal
import ..LineCableModels: basis, R, L, C, resistance, inductance, capacitance
import ..Grammar: AbstractCoreResult, observe, observables
using ..Materials: Material
import ..Validation
using Colors: HSL, RGB, RGBA, alpha, blue, green, red
import GeometryBasics
using GeometryBasics: Point2f
using Printf: @sprintf
using Statistics: mean
import Base: angle

# Abstract types and interfaces
include("interfaces.jl")
include("types.jl")
include("geometry/primitives.jl")
include("geometry/shapes.jl")
include("geometry/pose.jl")
include("geometry/resolve.jl")
include("design/region.jl")
include("design/stack.jl")
include("placement/patterns.jl")
include("placement/paths.jl")
include("placement/compaction.jl")
include("design/group.jl")
include("design/assembly.jl")
include("design/enclosure.jl")

# Submodule `BaseParams`
include("baseparams/BaseParams.jl")
using .BaseParams

include("nominaldata.jl")
include("design/cabledesign.jl")
include("cabledesign/cableconstants.jl")
include("cabledesign/dataframe.jl")

# Library
include("cableslibrary/cableslibrary.jl")
include("linecablesystem/linecablesystem.jl")

# Geometry and language protocols
include("geometry.jl")
include("preview/definitions.jl")
include("preview/materials.jl")
include("preview/cable.jl")
include("preview/cables.jl")
include("preview/system.jl")
include("preview/materialscale.jl")

public preview_shapes, preview_materials
public PreviewPolygon, PreviewReferenceLine, PreviewPayload

end # module DataModel
