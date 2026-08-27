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
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
using ...EarthProps: EarthModel
import ..Engine: AbstractEHEMFormulation, description

include("enforcelayer.jl")

end # module EHEM
