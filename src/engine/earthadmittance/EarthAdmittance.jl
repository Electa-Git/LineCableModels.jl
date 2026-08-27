"""
    LineCableModels.Engine.EarthAdmittance

Define formulations for the earth contribution to shunt admittance.

`IdealGround` sets the environmental potential coefficient to zero.
`Papadopoulos` evaluates the homogeneous-earth integral.

# Dependencies

$(IMPORTS)

"""
module EarthAdmittance

# Export public API
export IdealGround, Papadopoulos

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDSIGNATURES
#! explicit-imports: on
import ...LineCableModels: nominal
import ..Engine: EarthAdmittanceFormulation, description
import ..Engine: conductivity, bessel_difference
using QuadGK: quadgk

vacuum_permittivity(value) = one(value) * 88541878128 * (one(value) * 10)^(-22)
vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("homogeneous.jl")
include("idealground.jl")
include("base.jl")

end # module EarthAdmittance
