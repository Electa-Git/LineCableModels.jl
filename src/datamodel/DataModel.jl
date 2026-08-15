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
export Thickness, Diameter  # Type definitions
export CircStrands, RectStrands, Strip, Tubular, SectorParams, Sector  # Conductor types
export Semicon, Insulator, SectorInsulator  # Insulator types
export ConductorGroup, InsulatorGroup  # Group types
export CableComponent, CableDesign, CableConstants  # Cable design types
export CablePosition, LineCableSystem  # System types
export CablesLibrary, NominalData  # Support types
export trifoil_formation, flat_formation, get_outer_radius, MaxFill  # Helpers
export preview, equivalent

# Module-specific dependencies
using ..Commons
import ..PlotBuilder
import ..UnitHandler
import ..Commons: add!
import ..Commons: basis, R, L, C, resistance, inductance, capacitance
using ..Utils:
               resolve_T, to_certain, to_nominal, is_headless,
               is_in_testset, to_lower, to_upper
import ..Utils: coerce_to_T, to_lower
using ..Materials: Material
import ..Validation: Validation, sanitize, validate!, has_radii, has_temperature,
                     extra_rules, IntegerField, Positive, Finite, Normalized, IsA,
                     required_fields,
                     coercive_fields, keyword_fields, keyword_defaults, _kwdefaults_nt,
                     is_radius_input,
                     Nonneg, OneOf, Greater, PhysicalFillLimit, Satisfies
using Measurements
using DataFrames
using Colors
using LinearAlgebra
using GeometryBasics: Point, Point2f, Polygon
using Printf: @sprintf
using Statistics: mean
# Abstract types & interfaces
include("types.jl")
include("radii.jl")

# Submodule `BaseParams`
include("baseparams/BaseParams.jl")
using .BaseParams

# Constructors
include("macros.jl")
include("validation.jl")

# Conductors
include("strands_handler.jl")
include("circstrands.jl")
include("rectstrands.jl")
include("strip.jl")
include("tubular.jl")
include("conductorgroup.jl")
include("sector.jl")

# Insulators
include("insulator.jl")
include("semicon.jl")
include("insulatorgroup.jl")
include("sectorinsulator.jl")

# Groups
include("nominaldata.jl")
include("cablecomponent.jl")
include("cabledesign.jl")

# Library
include("cableslibrary.jl")
include("linecablesystem.jl")

# Helpers & overrides
include("helpers.jl")
include("io.jl")
include("typecoercion.jl")
include("plotspecs.jl")

"""
    preview(object; kwargs...)

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
    show_material_scale(; kwargs...)

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

# Aliases for backward compatibility
const WireArray = CircStrands
export WireArray

end # module DataModel
