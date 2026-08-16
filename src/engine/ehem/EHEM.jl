"""
    LineCableModels.Engine.EHEM

# Dependencies

$(IMPORTS)

"""
module EHEM

# Export public API
export EnforceLayer

# Module-specific dependencies
using ...Commons
using ...EarthProps: EarthModel
import ...Commons: get_description
import ..Engine: AbstractEHEMFormulation

include("enforcelayer.jl")

end # module EHEM
