"""
    LineCableModels.Engine.InsulationImpedance

# Dependencies

$(IMPORTS)

"""
module InsulationImpedance

# Export public API
export Lossless

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ...LineCableModels: description
import ..Engine: InsulationImpedanceFormulation

include("lossless.jl")

end # module InsulationImpedance
