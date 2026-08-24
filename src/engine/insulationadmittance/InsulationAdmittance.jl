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
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDSIGNATURES
import ..Engine: description
import ..Engine: conductivity
import ..Engine: InsulationAdmittanceFormulation

include("lossless.jl")
include("parallelrc.jl")

end # module InsulationAdmittance
