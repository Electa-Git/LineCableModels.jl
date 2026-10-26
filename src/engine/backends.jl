"""
$(TYPEDEF)

Represent the backend that executes the coaxial formulation in process.
"""
struct LCMBackend <: AbstractBackend end

"""
$(TYPEDEF)

Represent the backend that executes [`FEM`](@ref) through the Gmsh extension.
"""
struct GmshBackend <: AbstractBackend end

"""
    backend(F)

Return the default backend of a formulation type. Pass a backend first to
[`compute`](@ref) to choose another.

# Arguments

- `F`: formulation type.

# Returns

- A backend singleton.
"""
function backend end

backend(::Type{<:CoaxialFormulation}) = LCMBackend()
backend(::Type{<:FEM}) = GmshBackend()

defaults(::Type{LCMBackend}) = (trace = false,)

"""
$(TYPEDSIGNATURES)

Compute line parameters with the default backend from [`backend`](@ref).

# Arguments

- `problem`: cable system, earth, temperature, and frequencies.
- `formulation`: formulas and formulation controls.

# Keywords

- `kw...`: computation options of the default backend.

# Returns

- A [`LineParameters`](@ref) record.

# Errors

- Throws `ArgumentError` for an unknown computation option or an unavailable backend
  implementation.
"""
function compute(problem::LineParametersProblem, formulation::AbstractFormulation; kw...)
    return compute(backend(typeof(formulation)), problem, formulation; kw...)
end

"""
$(TYPEDSIGNATURES)

Build the computation options once and compute line parameters on the supplied backend.

# Arguments

- `backend`: backend singleton that executes the computation.
- `problem`: cable system, earth, temperature, and frequencies.
- `formulation`: formulas and formulation controls.

# Keywords

- `kw...`: computation options of `backend`.

# Returns

- A [`LineParameters`](@ref) record.

# Errors

- Throws `ArgumentError` for an unknown computation option or an unsupported formulation.
"""
function compute(backend::AbstractBackend, problem::LineParametersProblem,
        formulation::AbstractFormulation; kw...)
    return compute(backend, problem, formulation, ComputationOptions(typeof(backend); kw...))
end

# A backend implements the formulations it has a method for.
function compute(backend::AbstractBackend, ::LineParametersProblem,
        formulation::AbstractFormulation, ::ComputationOptions)
    throw(ArgumentError("$(nameof(typeof(backend))) does not implement $(nameof(typeof(formulation)))"))
end

# The Gmsh extension defines the finite-element computation.
function compute(::GmshBackend, ::LineParametersProblem, ::AbstractFormulation; kw...)
    throw(ArgumentError(
        "the Gmsh backend implements FEM and needs the Gmsh extension: run `using Gmsh` first"))
end
