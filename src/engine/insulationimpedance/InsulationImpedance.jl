"""
    LineCableModels.Engine.InsulationImpedance

# Dependencies

$(IMPORTS)

"""
module InsulationImpedance

# Export public API
export Lossless

# Module-specific dependencies
using ...Commons
import ...Commons: get_description
import ..Engine: InsulationImpedanceFormulation

include("lossless.jl")

end # module InsulationImpedance
