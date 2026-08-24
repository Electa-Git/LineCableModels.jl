"""
    LineCableModels.Engine.InsulationImpedance

Define series-impedance formulations for cable insulation.

# Dependencies

$(IMPORTS)

"""
module InsulationImpedance

# Export public API
export Lossless

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ..Engine: InsulationImpedanceFormulation, description

include("lossless.jl")

end # module InsulationImpedance
