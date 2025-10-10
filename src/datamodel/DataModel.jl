"""
	LineCableModels.DataModel

The [`DataModel`](@ref) module provides data structures, constructors and utilities for modeling power cables within the [`LineCableModels.jl`](index.md) package. This module includes definitions for various cable components, and visualization tools for cable designs.

# Overview

- Provides objects for detailed cable modeling with the [`CableDesign`](@ref) and supporting types: [`WireArray`](@ref), [`Strip`](@ref), [`Tubular`](@ref), [`Semicon`](@ref), and [`Insulator`](@ref).
- Includes objects for cable **system** modeling with the [`LineCableSystem`](@ref) type, and multiple formation patterns like trifoil and flat arrangements.
- Contains functions for calculating the base electric properties of all elements within a [`CableDesign`](@ref), namely: resistance, inductance (via GMR), shunt capacitance, and shunt conductance (via loss factor).
- Offers visualization tools for previewing cable cross-sections and system layouts.
- Provides a library system for storing and retrieving cable designs.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module DataModel

# Export public API
export Thickness, Diameter  # Type definitions
export WireArray, Strip, Tubular, SectorParams, Sector  # Conductor types
export Semicon, Insulator, SectorInsulator  # Insulator types
export ConductorGroup, InsulatorGroup  # Group types
export CableComponent, CableDesign  # Cable design types
export CablePosition, LineCableSystem  # System types
export CablesLibrary, NominalData  # Support types
export trifoil_formation, flat_formation  # Formation helpers
export preview, equivalent

# Module-specific dependencies
using ..Commons
import ..Commons: add!
using ..Utils:
	resolve_T, to_certain, to_nominal, is_headless,
	is_in_testset, to_lower, to_upper
import ..Utils: coerce_to_T, to_lower
using ..Materials: Material
import ..BackendHandler: set_backend!, ensure_backend!, current_backend_symbol,
	backend_available, renderfig, next_fignum
import ..PlotUIComponents: gl_screen, with_icon, MI_REFRESH, MI_SAVE, ICON_TTF
import ..Validation: Validation, sanitize, validate!, has_radii, has_temperature,
	extra_rules, IntegerField, Positive, Finite, Normalized, IsA, required_fields,
	coercive_fields, keyword_fields, keyword_defaults, _kwdefaults_nt, is_radius_input,
	Nonneg, OneOf
using Measurements
using DataFrames
using Colors
using Plots
using DisplayAs: DisplayAs
using LinearAlgebra
using Makie: Point, Point2f # otherwise will require adding GeometryBasics as a dependency
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
include("wirearray.jl")
include("strip.jl")
include("tubular.jl")
include("conductorgroup.jl")
include("sector.jl")

# Insulators
include("insulator.jl")
include("semicon.jl")
include("insulatorgroup.jl")
include("SectorInsulator.jl")


# Groups
include("nominaldata.jl")
include("cablecomponent.jl")
include("cabledesign.jl")

# Library
include("cableslibrary.jl")
include("linecablesystem.jl")

# Helpers & overrides
include("helpers.jl")
include("preview.jl")
include("io.jl")
include("typecoercion.jl")

end # module DataModel