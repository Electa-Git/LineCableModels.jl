"""
	LineCableModels.Engine.EHEM

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module EHEM

# Export public API
export EnforceLayer

# Module-specific dependencies
using ...Commons
using ...EarthProps: EarthModel
import ...Commons: get_description
import ..Engine: AbstractEHEMFormulation
using Measurements

include("enforcelayer.jl")

end # module EHEM