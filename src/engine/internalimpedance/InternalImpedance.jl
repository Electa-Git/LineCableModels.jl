"""
	LineCableModels.Engine.InternalImpedance

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module InternalImpedance

# Export public API
export ScaledBessel

# Module-specific dependencies
using ...Commons
import ...Commons: get_description
import ..Engine: InternalImpedanceFormulation
using Measurements
using LinearAlgebra
using ...UncertainBessels: besselix, besselkx
using ...Utils: _to_σ

include("scaledbessel.jl")

end # module InternalImpedance

