"""
    LineCableModels.EarthProps

Define static homogeneous and layered-earth descriptions, measured
frequency-dependent material relations, and equivalent homogeneous-earth
reductions required by line-parameter formulations.

# Public actions

- Construct and validate [`EarthLayer`](@ref), [`EarthModel`](@ref), and the
  ephemeral [`EarthMaterial`](@ref).
- Select measured frequency dependence through [`FD`](@ref).
- Select equivalent homogeneous-earth reductions through [`EHEM`](@ref).
- Extend an earth model with `add!`.
- Present earth data through the Base display protocol.
"""
module EarthProps

export EarthMaterial, EarthLayer, EarthModel
export build
export FD, EHEM

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: add!, build, validate
import ..LineCableModels: _construction
using ..Materials: AbstractMaterial
import ..Validation
import ..TextDisplay

include("earthmaterial.jl")
include("earthlayer.jl")
include("earthmodel.jl")

include("fd/FD.jl")
include("ehem/EHEM.jl")

include("base.jl")

end # module EarthProps
