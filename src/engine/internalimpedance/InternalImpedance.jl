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
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ..Engine: InternalImpedanceFormulation, description
import ..Engine: conductivity
using LinearAlgebra
import ..Engine: special_besselix, special_besselkx

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("scaledbessel.jl")

end # module InternalImpedance
