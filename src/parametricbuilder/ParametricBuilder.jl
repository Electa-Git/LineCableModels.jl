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
import ..Grammar: compute, observables, nominal, standard_uncertainty
using ..Grammar:
                 AbstractProblemDefinition, AbstractFormulation, AbstractProblemResult,
                 AbstractParametricResult, AbstractUncertaintyResult
import ..Materials
import ..Materials: Material
import ..DataModel
import ..EarthProps
import ..Engine

include("grid.jl")
include("gridspace.jl")
include("macros.jl")
include("results.jl")

include("materialdefinition.jl")
include("cablebuilderdefinition.jl")
include("positiondefinition.jl")
include("systembuilderdefinition.jl")

include("compute.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: WireEstimate, make_stranded, make_screened

end
