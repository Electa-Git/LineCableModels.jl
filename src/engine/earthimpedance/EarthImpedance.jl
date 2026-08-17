"""
    LineCableModels.Engine.EarthImpedance

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Papadopoulos

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ...LineCableModels: description, nominal
import ..Engine: EarthImpedanceFormulation
import ..Engine: conductivity, bessel_difference
using QuadGK: quadgk

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("homogeneous.jl")
include("base.jl")

end # module EarthImpedance
