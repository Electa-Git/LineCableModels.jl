"""
    LineCableModels.ParametricBuilder

Construct finite parameter spaces and materialise cable problems from explicit
`Grid` inputs.

# Overview

- Define deterministic and uncertainty-bearing finite sources.
- Compose sources with product or zip semantics.
- Materialise materials, cable parts, cable designs, positions, earth models,
  and line-parameter problems.
- Evaluate every materialised problem with `Combinatorial`.
- Project completed result spaces into finite spaces of downstream problems.
- Estimate stranded-conductor and wire-screen patterns.
"""
module ParametricBuilder

export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace
export has_uncertainty, nominal, uncertainty
export @gridspace, @relax
export build

export Combinatorial, ParametricProblem, ParametricResult
export result, project

export Material, Conductor, Insulator, Semiconductor, Filler
export CableDesign, LineCableSystem
export DiskDefinition, RectangleDefinition, EllipseDefinition
export SectorDefinition, AnnulusDefinition, ShellDefinition, PolygonDefinition, Pose2
export Region, Stack
export Group, Assembly
export Enclosure
export at, trefoil, hflat, vflat, Earth
export WireEstimate, make_stranded, make_screened

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random
using RequiredInterfaces: @required
using ..Grammar: @orchestrator
import ..LineCableModels: add!, build
import ..Grammar
import ..Grammar: compute, computation_options, computation_details, details,
                  nominal, uncertainty, check_core_result,
                  computation_owner
using ..Grammar:
                 AbstractProblemDefinition, AbstractFormulation, AbstractResultSpace,
                 AbstractParametricResult,
                 ComputationOptions, ComputationDetails
import ..Materials
import ..Materials: Material
import ..DataModel
import ..DataModel: DiskDefinition, RectangleDefinition, EllipseDefinition
import ..DataModel: SectorDefinition, AnnulusDefinition, ShellDefinition
import ..DataModel: PolygonDefinition, Pose2
import ..DataModel: Region, Stack
import ..DataModel: Group, Assembly
import ..DataModel: Enclosure
import ..DataModel: CableDesign, LineCableSystem
import ..EarthProps
import ..Engine

include("grid.jl")
include("gridspace.jl")
include("macros.jl")
include("results.jl")
include("project.jl")

include("material.jl")
include("geometry.jl")
include("physicaltree.jl")
include("conductor/Conductor.jl")
include("insulator/Insulator.jl")
include("semiconductor/Semiconductor.jl")
include("filler/Filler.jl")
include("positions.jl")
include("system.jl")

include("engine/cableconstants.jl")
include("traversal.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, make_stranded, make_screened

public AbstractProjectionDefinition, entitle, select, derive, materialize
public traverse, sample_uncertainty

end
