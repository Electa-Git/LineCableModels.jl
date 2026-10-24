"""
    LineCableModels.ParametricBuilder

Construct finite parameter spaces and materialize cable problems from explicit
`Grid` inputs.

# Overview

- Define deterministic and uncertainty-bearing finite sources.
- Compose sources with product or zip semantics.
- Declare cable parts, placements and systems with the construction vocabulary
  and macros.
- Transport completed result spaces into target-bearing downstream problem
  spaces.
- Estimate stranded-conductor and wire-screen patterns.
"""
module ParametricBuilder
import ..LineCableModels

export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace
export has_uncertainty, nominal, uncertainty
export build

export Material
export CableDesign, LineCableSystem
export Disk, Rectangle, Ellipse, Sector, Annulus, Polygon, Shell, Pose2
export Region, Stack
export Group, Assembly
export Enclosure
export terminal, core, stranded, milliken, rope, cores, tape, insulation, screen, sheath
export armor, bedding, jacket, filler, pipe, duct
export solid, shell, wires, layers, assembly
export capacity
export FillFactor
export at, trefoil, hflat, vflat, layer, homogeneous, EarthLayer, EarthModel
export @cable, @system, @earth, @terminal, @assembly, @pipe, @duct
export @at, @hflat, @vflat, @trefoil
export @distribute
export WireEstimate, estimate_stranding, estimate_screen

import Random
import ..LineCableModels: build, validate
import ..LineCableModels: Grid, AbsoluteError, DeterministicGrid, RelativeGrid
import ..LineCableModels: AbsoluteGrid, AbstractGrid, AbstractUncertainGrid
import ..LineCableModels: UncertainValue, Gridspace, ScalarOrGrid, has_uncertainty
import ..LineCableModels: Builder, target
import ..Commons: nominal, uncertainty
import ..Materials
import ..Materials: Material
import ..DataModel
import ..DataModel: Disk, Rectangle, Ellipse, Sector, Annulus, Polygon
import ..DataModel: Shell, Pose2
import ..DataModel: Region, Stack
import ..DataModel: Group, Assembly
import ..DataModel: Enclosure
import ..DataModel: capacity, FillFactor
import ..DataModel: CableDesign, LineCableSystem
import ..DataModel: trefoil, hflat, vflat
import ..Earth
using ..Earth: EarthLayer, EarthModel, layer, homogeneous
import ..TextDisplay

include("physicaltree.jl")
include("positions.jl")
include("construction_macros.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, estimate_stranding, estimate_screen

include("textdisplay.jl")

end
