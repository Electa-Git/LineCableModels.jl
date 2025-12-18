"""
	LineCableModels.Engine.Transforms

# Dependencies

$(IMPORTS)

# Exports

$(EXPORTS)
"""
module Transforms

# Export public API
export Fortescue

# Module-specific dependencies
using ...Commons
import ...Commons: get_description
import ...Utils: symtrans, symtrans!, offdiag_ratio, to_nominal
import ..Engine:
	AbstractTransformFormulation, LineParameters, SeriesImpedance, ShuntAdmittance,
	PhaseDomain, ModalDomain
#
using Measurements
using LinearAlgebra
# using GenericLinearAlgebra
using NLsolve


include("fortescue.jl")
include("eiglevenberg.jl")

function (F::AbstractTransformFormulation)(
	lp::LineParameters{Tc, U, ModalDomain},
) where {Tc <: COMPLEXSCALAR, U <: REALSCALAR}
	throw(
		ErrorException(
			"Not yet implemented: inverse $(nameof(typeof(F)))( ::LineParameters{<:COMPLEXSCALAR,<:REALSCALAR,ModalDomain} )",
		),
	)
end

end # module Transforms
