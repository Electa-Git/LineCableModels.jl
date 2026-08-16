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
import ...Commons: get_description, _typed_ε₀
using ...Utils: _to_σ
import ..Engine: InsulationAdmittanceFormulation

include("lossless.jl")
include("parallelrc.jl")

end # module InsulationAdmittance
