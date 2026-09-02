"""
    LineCableModels.EarthProps

Define static homogeneous and layered-earth descriptions, measured
frequency-dependent material relations, and equivalent homogeneous-earth
reductions required by line-parameter formulations.

# Public actions

- Declare earth descriptions with [`layer`](@ref) and [`homogeneous`](@ref),
  represented by [`EarthLayer`](@ref) and [`EarthModel`](@ref).
- Construct the ephemeral [`EarthMaterial`](@ref) used by the engine.
- Select measured frequency dependence through [`FD`](@ref).
- Select equivalent homogeneous-earth reductions through [`EHEM`](@ref).
- Build immutable earth models from complete ordered layer declarations.
- Present earth data through the Base display protocol.
"""
module EarthProps

export EarthMaterial, EarthLayer, EarthModel
export layer, homogeneous
export build
export FD, EHEM

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: build, validate
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
