"""
    LineCableModels.Engine.InsulationAdmittance

# Dependencies

$(IMPORTS)

"""
module InsulationAdmittance

# Export public API
export Lossless, ParallelRC

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDSIGNATURES
import ...LineCableModels: description
import ..Engine: conductivity
import ..Engine: InsulationAdmittanceFormulation

include("lossless.jl")
include("parallelrc.jl")

end # module InsulationAdmittance
