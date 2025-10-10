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

# Module-specific dependencies
using ...Commons
using ...Materials
using ...EarthProps
using ...DataModel
using ...Engine
import ...Engine: kronify, reorder_M, reorder_indices, merge_bundles!, FormulationSet,
	AbstractFormulationSet, AbstractImpedanceFormulation, AbstractAdmittanceFormulation,
	compute!
import ...Engine: AbstractFormulationOptions, LineParamOptions, build_options, _COMMON_SYMS
import ...DataModel: AbstractCablePart, AbstractConductorPart, AbstractInsulatorPart
using ...Utils:
	display_path, set_verbosity!, is_headless, to_nominal, symtrans!, symtrans,
	line_transpose!
using Measurements
using LinearAlgebra
using Colors
using Makie: Point, Point2f # otherwise will require adding GeometryBasics as a dependency
# FEM specific dependencies
using Gmsh
using GetDP
using GetDP: Problem, get_getdp_executable, add!


include("types.jl")
include("lineparamopts.jl")   # Line parameter options

# Include auxiliary files
include("meshtransitions.jl") # Mesh transition objects
include("problemdefs.jl")     # Problem definitions
include("workspace.jl")       # Workspace functions
include("encoding.jl")        # Tag encoding schemes
include("drawing.jl")         # Primitive drawing functions
include("identification.jl")  # Entity identification
include("mesh.jl")            # Mesh generation
include("materialprops.jl")       # Material handling
include("helpers.jl")         # Various utilities
include("visualization.jl")   # Visualization functions
include("space.jl")           # Domain creation functions
include("cable.jl")           # Cable geometry creation functions
include("solver.jl")          # Solver functions
include("base.jl")            # Base namespace extensions

end # module FEM
