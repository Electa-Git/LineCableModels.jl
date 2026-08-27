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

#! explicit-imports: off
# IMPORTS is expanded in the module docstring rather than called as Julia code.
using DocStringExtensions: IMPORTS
#! explicit-imports: on
using DocStringExtensions: TYPEDSIGNATURES
import ...LineCableModels: basis, nominal
import ..Engine:
                 AbstractTransformFormulation, LineParameters, description,
                 PhaseDomain, ModalDomain
using LinearAlgebra: Diagonal, I, checksquare, diag, eigen
using NLsolve: converged, nlsolve

include("matrices.jl")
include("fortescue.jl")
include("eiglevenberg.jl")
include("transform.jl")

end # module Transforms
