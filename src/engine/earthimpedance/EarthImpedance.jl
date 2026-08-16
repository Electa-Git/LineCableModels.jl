"""
    LineCableModels.Engine.EarthImpedance

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Papadopoulos

# Module-specific dependencies
using ...Commons
import ...Commons: get_description
import ..Engine: EarthImpedanceFormulation
using QuadGK: quadgk
using ...Utils: _to_σ, _bessel_diff, to_nominal

include("homogeneous.jl")
include("base.jl")

end # module EarthImpedance
