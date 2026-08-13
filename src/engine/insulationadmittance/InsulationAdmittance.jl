"""
    LineCableModels.Engine.InsulationAdmittance

# Dependencies

$(IMPORTS)

"""
module InsulationAdmittance

# Export public API
export Lossless, ParallelRC

# Module-specific dependencies
using ...Commons
import ...Commons: get_description
using ...Utils: _to_σ
import ..Engine: InsulationAdmittanceFormulation
using Measurements

include("lossless.jl")
include("parallelrc.jl")

end # module InsulationAdmittance
