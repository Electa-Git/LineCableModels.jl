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
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InsulationImpedanceFormulation, description

include("lossless.jl")

end # module InsulationImpedance
