"""
    LineCableModels.Engine.InternalImpedance

Define conductor internal-impedance formulations based on scaled modified
Bessel functions.

# Dependencies

$(IMPORTS)

"""
module InternalImpedance

# Export public API
export ScaledBessel

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
#! explicit-imports: on
import ..Engine: InternalImpedanceFormulation, description
import ..Engine: conductivity
import ..Engine: special_besselix, special_besselkx

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("scaledbessel.jl")

end # module InternalImpedance
