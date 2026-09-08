"""
    LineCableModels.Earth

Define static homogeneous and layered-earth descriptions, measured
frequency-dependent material relations, and equivalent homogeneous-earth
reductions required by line-parameter formulations.

# Public actions

- Declare earth descriptions with [`layer`](@ref) and [`homogeneous`](@ref),
  represented by [`EarthLayer`](@ref) and [`EarthModel`](@ref).
- Construct the ephemeral [`EarthMaterial`](@ref) used by the engine.
- Select measured frequency dependence through [`FrequencyDependent`](@ref).
- Select equivalent homogeneous-earth reductions through [`EquivalentHomogeneous`](@ref).
- Build immutable earth models from complete ordered layer declarations.
- Present earth data through the Base display protocol.
"""
module Earth

export AbstractEarthModel, AbstractEarthLayer, AbstractEarthMaterial, EarthMaterial,
       EarthLayer, EarthModel
export layer, homogeneous
export build
export FrequencyDependent, EquivalentHomogeneous

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using RequiredInterfaces: @required
import ..LineCableModels: build, validate
import ..LineCableModels: parameterize
using ..Materials: AbstractMaterial
import ..TextDisplay

"Supertype for materialized static earth-layer and earth-model descriptions."
abstract type AbstractEarthModel end

"Supertype for one static earth layer; retains the abstract earth-model contract."
abstract type AbstractEarthLayer <: AbstractEarthModel end

"Supertype for frequency-evaluated earth constitutive properties."
abstract type AbstractEarthMaterial <: AbstractMaterial end

@required AbstractEarthModel begin
    validate(::AbstractEarthModel)
end

include("earthmaterial.jl")
include("earthlayer.jl")
include("earthmodel.jl")

include("frequencydependent/FrequencyDependent.jl")
include("equivalenthomogeneous/EquivalentHomogeneous.jl")

include("base.jl")

end # module Earth
