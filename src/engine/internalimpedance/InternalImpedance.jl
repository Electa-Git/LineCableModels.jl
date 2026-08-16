"""
    LineCableModels.Engine.InternalImpedance

# Dependencies

$(IMPORTS)

"""
module InternalImpedance

# Export public API
export ScaledBessel

# Module-specific dependencies
using ...Commons
import ...Commons: get_description
import ..Engine: InternalImpedanceFormulation
using LinearAlgebra
using SpecialFunctions: besselix, besselkx
using ...Utils: _to_σ

include("scaledbessel.jl")

end # module InternalImpedance
