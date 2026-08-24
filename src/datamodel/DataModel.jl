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
import ..Units
import ..LineCableModels: add!, validate, maxfill, nominal, standard_uncertainty
import ..LineCableModels: basis, R, L, C, resistance, inductance, capacitance
import ..Grammar: AbstractProblemResult, observables
import ..LineCableModels: retired_fem_sector
import ..LineCableModels: SectorParams, Sector, SectorInsulator
using ..Materials: Material
import ..Validation
using ..Validation: IntegerField, Positive, Finite, IsA, Nonnegative, OneOf,
                    Greater, Less
using DataFrames
using Colors
using LinearAlgebra
using GeometryBasics: Point, Point2f, Polygon
using Printf: @sprintf
using Statistics: mean

# Abstract types and interfaces
include("interfaces.jl")
include("types.jl")

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
include("preview/system.jl")
include("preview/materialscale.jl")

end # module DataModel
