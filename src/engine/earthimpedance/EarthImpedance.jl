"""
    LineCableModels.Engine.EarthImpedance

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Papadopoulos
export Ametani, Deri, DirectNumericalIntegration, Lucca, Saad, Wedepohl

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDSIGNATURES
import ...LineCableModels: description, nominal
import ..Engine: EarthImpedanceFormulation
import ..Engine: conductivity, bessel_difference
using QuadGK: quadgk

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("homogeneous.jl")
include("reference.jl")
include("base.jl")

end # module EarthImpedance
