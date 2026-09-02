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
- Transport completed result spaces into target-bearing downstream problem
  spaces.
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
export result

export Material, Conductor, Insulator, Semiconductor
export CableDesign, LineCableSystem
export Disk, Rectangle, Ellipse, Sector, Annulus, Polygon, RoundedSector, Shell, Pose2
export Region, Stack
export Group, Assembly
export Enclosure
export terminal, core, stranded, rope, cores, tape, insulation, screen, sheath
export armor, bedding, jacket, filler, pipe, duct
export solid, shell, wires, layers, assembly
export capacity, Hexa
export FillFactor, DiameterFactor, TabulatedCompaction, AffineCompaction
export at, trefoil, hflat, vflat, layer, homogeneous, EarthLayer, EarthModel
export @cable, @system, @earth, @terminal, @assembly, @pipe, @duct
export @at, @hflat, @vflat, @trefoil
export @distribute
export WireEstimate, make_stranded, make_screened

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random
import ..LineCableModels: add!, build, Gridpoint, validate
import ..LineCableModels: _construction, _construction_axis, _finite_construction
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
import ..DataModel: Disk, Rectangle, Ellipse, Sector, Annulus, Polygon
import ..DataModel: RoundedSector, Shell, Pose2
import ..DataModel: Region, Stack
import ..DataModel: Group, Assembly
import ..DataModel: Enclosure
import ..DataModel: capacity, Hexa, FillFactor, DiameterFactor
import ..DataModel: TabulatedCompaction, AffineCompaction
import ..DataModel: CableDesign, LineCableSystem
import ..EarthProps
using ..EarthProps: EarthLayer, EarthModel, layer, homogeneous
import ..Engine
import ..TextDisplay

include("grid.jl")
include("gridspace.jl")
include("macros.jl")
include("results.jl")

include("material.jl")
include("geometry.jl")
include("physicaltree.jl")
include("conductor/Conductor.jl")
include("insulator/Insulator.jl")
include("semiconductor/Semiconductor.jl")
include("positions.jl")
include("construction_macros.jl")
include("system.jl")

include("engine/cableconstants.jl")
include("traversal.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, make_stranded, make_screened

include("textdisplay.jl")

public materialize, traverse, sample_uncertainty

end
