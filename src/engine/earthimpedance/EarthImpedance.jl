"""
    LineCableModels.Engine.EarthImpedance

Define earth-return impedance formulations. Engine evaluates `Papadopoulos`
directly. External implementations may support the formula identities `Deri`,
`Wedepohl`, `Saad`, `Ametani`, and `Lucca`.

# Dependencies

$(IMPORTS)

"""
module EarthImpedance

# Export public API
export Papadopoulos
export Ametani, Deri, Lucca, Saad, Wedepohl

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDSIGNATURES
import ...LineCableModels: nominal
import ..Engine: EarthImpedanceFormulation, description
import ..Engine: conductivity, bessel_difference
using QuadGK: quadgk

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("homogeneous.jl")
include("formulations.jl")
include("base.jl")

end # module EarthImpedance
