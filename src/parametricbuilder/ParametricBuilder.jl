module ParametricBuilder

export Grid, AbsoluteError, DeterministicGrid, RelativeGrid, AbsoluteGrid
export AbstractGrid, AbstractUncertainGrid, UncertainValue
export AbstractSpec, Gridspace, Configuration, configurations, materialize
export has_uncertainty, configuration_manifest, nominal, standard_uncertainty
export @gridspace, @relax

export Material, Conductor, Insulator, CableBuilder
export at, trifoil, hflat, vflat, Earth, SystemBuilder
export make_stranded, make_screened

using ..Commons: f₀, T₀, SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS
import ..Commons: add!
import ..Materials
import ..Materials: Material
import ..DataModel
import ..EarthProps
import ..Engine

include("gridspace/grid.jl")
include("gridspace/gridspace.jl")
include("gridspace/macros.jl")

include("materialspec.jl")
include("cablebuilderspec.jl")
include("positionspec.jl")
include("systembuilderspec.jl")

include("wirepatterns/WirePatterns.jl")
using .WirePatterns: make_stranded, make_screened

end
