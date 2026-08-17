"""
    LineCableModels.Engine.InternalImpedance

# Dependencies

$(IMPORTS)

"""
module InternalImpedance

# Export public API
export ScaledBessel

# Module-specific dependencies
using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ...LineCableModels: description
import ..Engine: InternalImpedanceFormulation
import ..Engine: conductivity
using LinearAlgebra
import ..Engine: special_besselix, special_besselkx

vacuum_permeability(value) = one(value) * 4 * (one(value) * π) * (one(value) * 10)^(-7)

include("scaledbessel.jl")

end # module InternalImpedance
