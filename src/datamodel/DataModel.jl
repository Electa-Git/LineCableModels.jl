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
export build, homogenize
export CablesLibrary, DatasheetInfo, catalogue
export trefoil_formation, flat_formation, outer_radius
export AbstractShape, AbstractPrimitive
export AbstractCablePart, Region, Stack
export Group, Assembly
export Enclosure
export Disk, Rectangle, Ellipse, Sector, Annulus, Polygon, RoundedSector, Shell
export Pose2
export EmptyBoundary
export resolve, boundary, area, perimeter, centroid, support, r_in, r_ex, thickness
export tessellate
export Ring, Polar, Fill, Lattice, capacity, placements
export FillFactor, DiameterFactor, TabulatedCompaction, AffineCompaction
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
import ..TextDisplay
import ..LineCableModels: add!, build, homogenize, validate, nominal
import ..LineCableModels: _construction
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
include("geometry/pose.jl")
include("geometry/primitives.jl")
include("geometry/shell.jl")
include("geometry/roundedsector.jl")
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

include("design/cabledesign.jl")
include("flatten.jl")
include("cabledesign/cableconstants.jl")

# Library
include("cableslibrary/datasheetinfo.jl")
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

# Bounded human-readable representations for the completed physical grammar.
include("textdisplay.jl")

public preview_shapes, preview_materials
public PreviewPolygon, PreviewReferenceLine, PreviewPayload
public RoundedSectorShape, ShellShape
public flatten

end # module DataModel
