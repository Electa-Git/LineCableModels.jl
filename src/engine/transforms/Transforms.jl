"""
    LineCableModels.Engine.Transforms

# Dependencies

$(IMPORTS)

"""
module Transforms

# Export public API
export Fortescue

# Module-specific dependencies
using ...Commons
import ...Commons: get_description, PhaseDomain, ModalDomain
import ...Commons: basis
import ...Utils: symtrans, symtrans!, offdiag_ratio, to_nominal
import ..Engine:
                 AbstractTransformFormulation, LineParameters, SeriesImpedance,
                 ShuntAdmittance

#
using Measurements
using LinearAlgebra
# using GenericLinearAlgebra
using NLsolve

include("fortescue.jl")
include("eiglevenberg.jl")

function (F::AbstractTransformFormulation)(
        lp::LineParameters{
        Tc, U, ModalDomain, Basis},
) where {Tc <: COMPLEXSCALAR, U <: REALSCALAR, Basis}
    throw(
        ErrorException(
        "Not yet implemented: inverse $(nameof(typeof(F)))( ::LineParameters{<:COMPLEXSCALAR,<:REALSCALAR,ModalDomain} )",
    ),
    )
end

end # module Transforms
