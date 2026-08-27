"""
    LineCableModels.Engine.InsulationAdmittance

Define shunt-admittance formulations for cable insulation. `Lossless` retains
capacitance only. `ParallelRC` retains capacitance and dielectric conductance.

# Dependencies

$(IMPORTS)

"""
module InsulationAdmittance

# Export public API
export Lossless, ParallelRC

# Module-specific dependencies
#! explicit-imports: off
# IMPORTS is expanded in this module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: description
import ..Engine: conductivity
import ..Engine: InsulationAdmittanceFormulation

include("lossless.jl")
include("parallelrc.jl")

end # module InsulationAdmittance
