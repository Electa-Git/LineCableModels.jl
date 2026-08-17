"""
    LineCableModels.Engine.EarthAdmittance

# Dependencies

$(IMPORTS)

"""
module EarthAdmittance

# Export public API
export Papadopoulos

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ...LineCableModels: description, nominal
import ..Engine: EarthAdmittanceFormulation
import ..Engine: conductivity, bessel_difference
using QuadGK: quadgk

vacuum_permittivity(value) = one(value) * 88541878128 * (one(value) * 10)^(-22)
vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("homogeneous.jl")
include("base.jl")

end # module EarthAdmittance
