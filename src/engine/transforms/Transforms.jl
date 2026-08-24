"""
    LineCableModels.Engine.Transforms

Define phase/modal transformations and matrix normalisation operations used by
analytical line-parameter calculations.

# Public actions

- Transform [`LineParameters`](@ref) with an owned formulation.
- Enforce reciprocity and ideal transposition on square parameter matrices.

# Dependencies

$(IMPORTS)
"""
module Transforms

export Fortescue

using DocStringExtensions: IMPORTS, TYPEDSIGNATURES
import ...LineCableModels: basis, nominal
import ..Engine:
                 AbstractTransformFormulation, LineParameters, SeriesImpedance,
                 ShuntAdmittance, description, PhaseDomain, ModalDomain
using LinearAlgebra
using NLsolve

include("matrices.jl")
include("fortescue.jl")
include("eiglevenberg.jl")
include("transform.jl")

end # module Transforms
