"""
    LineCableModelsGmshExt

Declare and validate the computation options of the Gmsh backend for finite elements.
Build a cable cross-section with shared contact points in the caller's Gmsh session.
"""
module LineCableModelsGmshExt

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using Gmsh: Gmsh, gmsh
using Printf: @sprintf
import LineCableModels.Commons: build, compute, defaults, select, validate
using LineCableModels.Commons: ComputationOptions, resolve
using LineCableModels.Materials: Material
import LineCableModels.DataModel: area
using LineCableModels.DataModel: AbstractShape, AbstractPrimitive, Circle, Annulus,
                                 Rectangle, Offset, AnnularSector, DifferenceShape,
                                 PlacedRegion, Placement, support, r_in, r_ex
import LineCableModels.CAD: draw
using LineCableModels.CAD: Segment, Arc, Curve, Contour
using LineCableModels.Engine: GmshBackend, FEM, LineParametersProblem, admittivity,
                              resistivity, reductions, properties, vacuum_permittivity,
                              vacuum_permeability

public GetdpArtifact

include("options.jl")
include("model.jl")
include("geometry.jl")
include("mesh.jl")
include("getdp.jl")

"""
$(TYPEDSIGNATURES)

Return the declared computation options of the Gmsh backend.

# Arguments

- `GmshBackend`: backend type.

# Returns

- `mesh_size_factor`: factor that scales every native resolution target
  \\[dimensionless\\]. Larger values coarsen the mesh.
- `domain_size_factor`: multiplier of the earth length scale for the domain's half-width
  \\[dimensionless\\].
- `pml_reflection`: nominal attenuation target of the PML in the open interval (0, 1)
  \\[dimensionless\\].
- `workers`: number of Gmsh and GetDP processes at once.
- `solver_threads`: number of BLAS and OpenMP threads of each GetDP process.
- `keep_run_directory`: whether to keep the run directory.
- `trace`: whether to keep the primitive matrices `Z` \\[Ω/m\\] and `P` \\[Ω·m\\] by
  terminal in the result.
- `getdp_executable`: GetDP executable as a path or [`GetdpArtifact`](@ref).
- `gmsh_verbosity`: Gmsh's `-v` level, from 0 to 5.
- `getdp_verbosity`: GetDP's `-v` level, from 0 to 5.
- `overrides`: native expert controls by native name, each at its default in
  [`OVERRIDES`](@ref).
"""
function defaults(::Type{GmshBackend})
    return (mesh_size_factor = 1.0,
        domain_size_factor = 2.0,
        pml_reflection = 1e-3,
        workers = clamp(Sys.CPU_THREADS ÷ 4, 1, 8),
        solver_threads = 1,
        keep_run_directory = false,
        trace = false,
        getdp_executable = GetdpArtifact(),
        gmsh_verbosity = 2,
        getdp_verbosity = 2,
        overrides = map(entry -> entry.default, OVERRIDES))
end

"""
$(TYPEDSIGNATURES)

Build and validate the computation options once, then pass them with the problem and
formulation to the backend's computation method.

# Arguments

- `backend`: Gmsh backend singleton.
- `problem`: cable system, earth, temperature, and frequencies.
- `formulation`: finite-element formulation and its material laws.

# Keywords

- `kw...`: computation options declared by [`defaults`](@ref).

# Errors

- Throws `ArgumentError` for an unknown or invalid computation option or an unsupported
  formulation.
"""
function compute(backend::GmshBackend, problem::LineParametersProblem,
        formulation::FEM; kw...)
    options = ComputationOptions(GmshBackend; kw...)
    return compute(backend, problem, formulation, options)
end

end
