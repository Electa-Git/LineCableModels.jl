"This package's own implementation of a formulation, which computes it in process."
struct LCMBackend <: AbstractBackend end

"Backend that runs the finite-element computation of [`FEM`](@ref), in the Gmsh extension."
struct GmshBackend <: AbstractBackend end

"""
    backend(F)

Return the backend that runs the formulations of type `F`.
"""
function backend end

backend(::Type{<:CoaxialFormulation}) = LCMBackend()
backend(::Type{<:FEM}) = GmshBackend()

defaults(::Type{LCMBackend}) = (trace = false,)

"""
$(TYPEDSIGNATURES)

Compute the line parameters of `problem` by `formulation`, on the backend that
[`backend`](@ref) gives for its type. The keywords are the options of that backend.
"""
function compute(problem::LineParametersProblem, formulation::AbstractFormulation; kw...)
    return compute(backend(typeof(formulation)), problem, formulation; kw...)
end

# The Gmsh extension defines the finite-element computation.
function compute(::GmshBackend, ::LineParametersProblem, ::AbstractFormulation; kw...)
    throw(ArgumentError(
        "the finite-element backend needs the Gmsh extension: run `using Gmsh` first"))
end
