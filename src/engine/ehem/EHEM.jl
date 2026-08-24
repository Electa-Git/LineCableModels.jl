"""
    LineCableModels.Engine.EHEM

Define equivalent-homogeneous-earth reductions for evaluated layered-earth
properties.

# Dependencies

$(IMPORTS)

"""
module EHEM

# Export public API
export EnforceLayer

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using ...EarthProps: EarthModel
import ..Engine: AbstractEHEMFormulation, description

include("enforcelayer.jl")

end # module EHEM
