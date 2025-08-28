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
export WireArray, Strip, Tubular  # Conductor types
export Semicon, Insulator  # Insulator types
export ConductorGroup, InsulatorGroup  # Group types
export CableComponent, CableDesign  # Cable design types
export CablePosition, LineCableSystem  # System types
export CablesLibrary, NominalData  # Support types
export trifoil_formation, flat_formation  # Formation helpers

# Load common dependencies
using ..LineCableModels
include("utils/commondeps.jl")

# Module-specific dependencies
using Measurements
using DataFrames
using Colors
using Plots
using DisplayAs: DisplayAs
using ..Utils
using ..Materials
using ..EarthProps
using ..Validation
import ..Validation: sanitize, validate!, has_radii, has_temperature, extra_rules, IntegerField, Positive, Finite, Normalized, IsA, required_fields, coercive_fields, keyword_fields, keyword_defaults, _kwdefaults_nt

# Abstract types & constructors
include("datamodel/types.jl")
include("datamodel/macros.jl")
include("datamodel/validation.jl")
include("datamodel/radii.jl")

# Conductors
include("datamodel/wirearray.jl")
include("datamodel/strip.jl")
include("datamodel/tubular.jl")
include("datamodel/conductorgroup.jl")

# Insulators
include("datamodel/insulator.jl")
include("datamodel/semicon.jl")
include("datamodel/insulatorgroup.jl")

# Submodule `BaseParams`
include("datamodel/BaseParams.jl")
@force using .BaseParams
@reexport using .BaseParams

# Groups
include("datamodel/nominaldata.jl")
include("datamodel/cablecomponent.jl")
include("datamodel/cabledesign.jl")

# Library
include("datamodel/cableslibrary.jl")
include("datamodel/linecablesystem.jl")

# Helpers & overrides
include("datamodel/helpers.jl")
include("datamodel/preview.jl")
include("datamodel/io.jl")
include("datamodel/typecoercion.jl")

end # module DataModel