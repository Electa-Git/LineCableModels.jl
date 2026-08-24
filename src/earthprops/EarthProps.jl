"""
    LineCableModels.EarthProps

Define static physical descriptions of homogeneous and layered earth. Analysis
frequencies and frequency-dependent property formulations belong to
`LineCableModels.Engine`.

# Public actions

- Construct and validate [`EarthLayer`](@ref) and [`EarthModel`](@ref).
- Extend an earth model with `add!`.
- Present earth data through Base and DataFrames protocols.
"""
module EarthProps

export EarthLayer, EarthModel

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using DataFrames
import ..LineCableModels: add!, validate
import ..Validation

include("earthlayer.jl")
include("earthmodel.jl")
include("dataframe.jl")
include("base.jl")

end # module EarthProps
