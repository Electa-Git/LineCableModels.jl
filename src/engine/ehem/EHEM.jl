"""
    LineCableModels.Engine.EHEM

# Dependencies

$(IMPORTS)

"""
module EHEM

# Export public API
export EnforceLayer

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using ...EarthProps: EarthModel
import ...LineCableModels: description
import ..Engine: AbstractEHEMFormulation

include("enforcelayer.jl")

end # module EHEM
