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
- Estimate stranded-conductor and wire-screen patterns.
"""
module ParametricBuilder

export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export Gridspace
export has_uncertainty, nominal, standard_uncertainty
export @gridspace, @relax

export Combinatorial, ParametricProblem, ParametricResult
export result

export Material, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export WireEstimate, make_stranded, make_screened

using DocStringExtensions: SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
using Random
import ..LineCableModels: add!, maxfill
import ..Grammar
import ..Grammar: compute, computation_options, computation_details, details,
                  nominal, standard_uncertainty, check_core_result
using ..Grammar:
                 AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                 AbstractResultSpace, AbstractParametricResult, AbstractUncertaintyResult,
                 ComputationOptions, ComputationDetails
import ..Materials
import ..Materials: Material
import ..DataModel
import ..EarthProps
import ..Engine

include("grid.jl")
include("gridspace.jl")
include("macros.jl")
include("results.jl")

include("material.jl")
include("parts.jl")
include("conductor/Conductor.jl")
include("insulator/Insulator.jl")
include("cablebuilder.jl")
include("positions.jl")
include("systembuilder.jl")

include("compute.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, make_stranded, make_screened

end
