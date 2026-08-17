"""
    LineCableModels.DataModel

The [`DataModel`](@ref) module provides data structures, constructors and utilities for modeling power cables within the [`LineCableModels.jl`](index.md) package. This module includes definitions for various cable components, and visualization tools for cable designs.

# Overview

- Provides objects for detailed cable modeling with the [`CableDesign`](@ref) and supporting types: [`CircStrands`](@ref), [`Strip`](@ref), [`Tubular`](@ref), [`Semicon`](@ref), and [`Insulator`](@ref).
- Includes objects for cable **system** modeling with the [`LineCableSystem`](@ref) type, and multiple formation patterns like trifoil and flat arrangements.
- Contains functions for calculating the base electric properties of all elements within a [`CableDesign`](@ref), namely: resistance, inductance (via GMR), shunt capacitance, and shunt conductance (via loss factor).
- Offers visualization tools for previewing cable cross-sections and system layouts.
- Provides a library system for storing and retrieving cable designs.

# Dependencies

$(IMPORTS)

"""
module DataModel

# Export public API
export CircStrands, RectStrands, Strip, Tubular  # Conductor types
export Semicon, Insulator  # Insulator types
export SectorParams, Sector, SectorInsulator  # Removed API tombstones
export ConductorGroup, InsulatorGroup  # Group types
export CableComponent, CableDesign, CableConstants  # Cable design types
export CablePosition, LineCableSystem  # System types
export CablesLibrary, NominalData  # Support types
export trifoil_formation, flat_formation, outer_radius, maxfill  # Geometry
export ncables, nphases
export preview, equivalent

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES,
                           FUNCTIONNAME, METHODLIST
import ..PlotBuilder
import ..UnitHandler
import ..LineCableModels: add!, validate, maxfill, nominal, standard_uncertainty
import ..LineCableModels: basis, R, L, C, resistance, inductance, capacitance
import ..LineCableModels: ncables, nphases
import ..LineCableModels: retired_fem_sector
import ..LineCableModels: SectorParams, Sector, SectorInsulator
using ..Materials: Material
import ..Validation
using ..Validation: IntegerField, Positive, Finite, IsA, Nonnegative, OneOf,
                    Greater, Less, PhysicalFillLimit
using DataFrames
using Colors
using LinearAlgebra
using GeometryBasics: Point, Point2f, Polygon
using Printf: @sprintf
using Statistics: mean
# Abstract types & interfaces
include("types.jl")

# Submodule `BaseParams`
include("baseparams/BaseParams.jl")
using .BaseParams

# Constructors
# Conductors
include("strands_handler.jl")
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
include("cablecomponent.jl")
include("cabledesign.jl")

# Library
include("cableslibrary.jl")
include("linecablesystem.jl")

# Helpers & overrides
include("geometry.jl")
include("io.jl")
include("plotspecs.jl")

"""
$(TYPEDSIGNATURES)

Preview a cable design or cable system with a loaded Makie backend.

Load `CairoMakie`, `GLMakie`, or `WGLMakie` before calling this function.
"""
function preview end

function preview(args...; kwargs...)
    throw(
        ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling preview.",
    ),
    )
end

"""
$(TYPEDSIGNATURES)

Display the resistivity, permeability, and permittivity color scales used by
[`preview`](@ref). This internal helper supports preview development and visual
regression tests. Load a Makie backend before calling it.
"""
function show_material_scale end

function show_material_scale(args...; kwargs...)
    throw(
        ArgumentError(
        "Plotting is optional. Load CairoMakie, GLMakie, or WGLMakie before calling show_material_scale.",
    ),
    )
end

end # module DataModel
