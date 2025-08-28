"""
    LineCableModels.Engine.FEM

The [`FEM`](@ref) module provides functionality for generating geometric meshes for cable cross-sections, assigning physical properties, and preparing the system for electromagnetic simulation within the [`LineCableModels.jl`](index.md) package.

# Overview

- Defines core types [`FEMFormulation`](@ref), and [`FEMWorkspace`](@ref) for managing simulation parameters and state.
- Implements a physical tag encoding system (CCOGYYYYY scheme for cable components, EPFXXXXX for domain regions).
- Provides primitive drawing functions for geometric elements.
- Creates a two-phase workflow: creation → fragmentation → identification.
- Maintains all state in a structured [`FEMWorkspace`](@ref) object.

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module FEM

# Export public API
export MeshTransition, calc_domain_size
export compute!, preview_results
export FormulationSet, Electrodynamics, Darwin

# Load common dependencies
using ...LineCableModels
include("../utils/commondeps.jl")

# Module-specific dependencies
using Measurements
using LinearAlgebra
using Colors
using ...Utils
using ...Materials
using ...EarthProps
using ...DataModel
using ...Engine
import ...DataModel: AbstractCablePart, AbstractConductorPart, AbstractInsulatorPart
import ..Engine: AbstractFormulationSet, AbstractFormulationOptions, AbstractImpedanceFormulation, AbstractAdmittanceFormulation, compute!
import ...Utils: display_path

# FEM specific dependencies
using Gmsh
using GetDP
using GetDP: Problem, get_getdp_executable, add!

"""
$(TYPEDEF)

Abstract base type for workspace containers in the FEM simulation framework.
Workspace containers maintain the complete state of a simulation, including 
intermediate data structures, identification mappings, and results.

Concrete implementations should provide state tracking for all phases of the 
simulation process from geometry creation through results analysis.
"""
abstract type AbstractWorkspace end

"""
$(TYPEDEF)

Abstract type for entity data to be stored within the FEMWorkspace.
"""
abstract type AbstractEntityData end

"""
$(TYPEDEF)

Core entity data structure containing common properties for all entity types.

$(TYPEDFIELDS)
"""
struct CoreEntityData
    "Encoded physical tag \\[dimensionless\\]."
    physical_group_tag::Int
    "Name of the elementary surface."
    elementary_name::String
    "Target mesh size \\[m\\]."
    mesh_size::Float64
end

"""
$(TYPEDEF)

Entity data structure for cable parts.

$(TYPEDFIELDS)
"""
struct CablePartEntity{T<:AbstractCablePart} <: AbstractEntityData
    "Core entity data."
    core::CoreEntityData
    "Reference to original cable part."
    cable_part::T
end

"""
$(TYPEDEF)

Entity data structure for domain surfaces external to cable parts.

$(TYPEDFIELDS)
"""
struct SurfaceEntity <: AbstractEntityData
    "Core entity data."
    core::CoreEntityData
    "Material properties of the domain."
    material::Material
end

"""
$(TYPEDEF)

Entity data structure for domain curves (boundaries and layer interfaces).

$(TYPEDFIELDS)
"""
struct CurveEntity <: AbstractEntityData
    "Core entity data."
    core::CoreEntityData
    "Material properties of the domain."
    material::Material
end

"""
$(TYPEDEF)

Entity container that associates Gmsh entity with metadata.

$(TYPEDFIELDS)
"""
struct GmshObject{T<:AbstractEntityData}
    "Gmsh entity tag (will be defined after boolean fragmentation)."
    tag::Int32
    "Entity-specific data."
    data::T
end

"""
$(TYPEDSIGNATURES)

Constructs a [`GmshObject`](@ref) instance with automatic type conversion.

# Arguments

- `tag`: Gmsh entity tag (will be converted to Int32)
- `data`: Entity-specific data conforming to [`AbstractEntityData`](@ref)

# Returns

- A [`GmshObject`](@ref) instance with the specified tag and data.

# Notes

This constructor automatically converts any integer tag to Int32 for compatibility with the Gmsh C API, which uses 32-bit integers for entity tags.

# Examples

```julia
# Create domain entity with tag and data
core_data = CoreEntityData([0.0, 0.0, 0.0])
domain_data = SurfaceEntity(core_data, material)
entity = $(FUNCTIONNAME)(1, domain_data)
```
"""
function GmshObject(tag::Integer, data::T) where {T<:AbstractEntityData}
    return GmshObject{T}(Int32(tag), data)
end

mutable struct Darwin <: AbstractImpedanceFormulation
    problem::GetDP.Problem
    resolution_name::String

    function Darwin()
        return new(GetDP.Problem(), "Darwin")
    end

end

mutable struct Electrodynamics <: AbstractAdmittanceFormulation
    problem::GetDP.Problem
    resolution_name::String

    function Electrodynamics()
        return new(GetDP.Problem(), "Electrodynamics")
    end
end

# Include auxiliary files
include("fem/meshtransitions.jl") # Mesh transition objects
include("fem/problemdefs.jl")     # Problem definitions
include("fem/workspace.jl")       # Workspace functions
include("fem/encoding.jl")        # Tag encoding schemes
include("fem/drawing.jl")         # Primitive drawing functions
include("fem/identification.jl")  # Entity identification
include("fem/mesh.jl")            # Mesh generation
include("fem/materials.jl")       # Material handling
include("fem/helpers.jl")         # Various utilities
include("fem/visualization.jl")   # Visualization functions
include("fem/space.jl")           # Domain creation functions
include("fem/cable.jl")           # Cable geometry creation functions
include("fem/solver.jl")          # Solver functions

end # module FEM
