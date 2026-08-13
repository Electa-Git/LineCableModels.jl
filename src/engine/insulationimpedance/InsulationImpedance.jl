"""
	LineCableModels.Engine.InsulationImpedance

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module InsulationImpedance

# Export public API
export Lossless


# Module-specific dependencies
using ...Commons
import ...Commons: get_description
import ..Engine: InsulationImpedanceFormulation
using Measurements

include("lossless.jl")

end # module InsulationImpedance
