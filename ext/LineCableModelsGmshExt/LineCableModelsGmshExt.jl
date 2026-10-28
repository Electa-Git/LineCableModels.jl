"""
    LineCableModelsGmshExt

Compute finite-element series impedances and shunt admittances with Gmsh and GetDP.
Build and mesh the cross-section, solve every frequency, and reduce the terminal matrices.
"""
module LineCableModelsGmshExt

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using Gmsh: Gmsh, gmsh
using Printf: @sprintf
using LinearAlgebra: I, ldiv!, lu!
using Pkg.Artifacts: ensure_artifact_installed
import LineCableModels.Commons: build, compute, defaults, select, validate
using LineCableModels.Commons: ComputationOptions, resolve
using LineCableModels.Materials: Material
import LineCableModels.DataModel: area
using LineCableModels.DataModel: AbstractShape, AbstractPrimitive, Circle, Annulus,
                                 Rectangle, Offset, AnnularSector, DifferenceShape,
                                 PlacedRegion, Placement, support, r_in, r_ex
import LineCableModels.CAD: draw
using LineCableModels.CAD: Segment, Arc, Curve, Contour
import LineCableModels.Engine: solve!
using LineCableModels.Engine: GmshBackend, FEM, LineParametersProblem, LineParameters,
                              bundle!, eliminate!, transposition!, admittivity,
                              resistivity, reductions, properties, nominal,
                              vacuum_permittivity,
                              vacuum_permeability

public GetdpArtifact

include("options.jl")
include("model.jl")
include("geometry.jl")
include("mesh.jl")
include("getdp.jl")
include("result.jl")

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

"""
$(TYPEDSIGNATURES)

Validate the problem, build and mesh its cross-section, solve the finite-element problem at
every frequency, and apply the selected conductor reductions.

# Arguments

- `GmshBackend`: backend that runs Gmsh and GetDP.
- `problem`: cables, earth, temperature \\[°C\\], and frequencies \\[Hz\\].
- `formulation`: finite-element material laws and reduction controls.
- `options`: mesh, solver, trace, and run-directory controls.

# Returns

- A `LineParameters` record with series impedances \\[Ω/m\\], shunt admittances
  \\[S/m\\], terminal labels, and primitive matrices when traced.

# Errors

- Throws for invalid problem data or unavailable native executables.
- Retains the run directory on failure. Native process errors identify the frequency and
  log, with its last twenty lines.
"""
function compute(::GmshBackend, problem::LineParametersProblem{T}, formulation::FEM,
        options::ComputationOptions{GmshBackend}) where {T <: Real}
    mapped = map(x -> Float64(nominal(x)), problem)
    validate(mapped, formulation)
    model = FieldModel(mapped, formulation)
    path = directory(mapped)
    for name in (joinpath("input", "getdp"), "mesh", joinpath("results", "jobs"),
        "logs", "work")
        mkpath(joinpath(path, name))
    end
    for source in readdir(joinpath(@__DIR__, "getdp"); join = true)
        cp(source, joinpath(path, "input", "getdp", basename(source)))
    end
    write(joinpath(path, "input", "model_data.pro"), model, options)
    entry = "Include \"model_data.pro\";\n" *
            "Include \"getdp/parameters.pro\";\n" *
            "Include \"physical.geo\";\n" *
            "Include \"getdp/geometry.geo\";\n" *
            "Include \"getdp/mesh.geo\";\n"
    write(joinpath(path, "input", "model.geo"), entry)
    owned = !Bool(gmsh.isInitialized())
    owned && gmsh.initialize(String[], false, false)
    current = gmsh.model.getCurrent()
    models = Set(gmsh.model.list())
    terminal = gmsh.option.getNumber("General.Terminal")
    verbosity = gmsh.option.getNumber("General.Verbosity")
    result = try
        gmsh.option.setNumber("General.Terminal", options.gmsh_verbosity > 0 ? 1 : 0)
        gmsh.option.setNumber("General.Verbosity", options.gmsh_verbosity)
        geometry = build(Geometry, model, options)
        write(joinpath(path, "input", "physical.geo"), geometry)
        mesh(path, model, options)
        plan = reductions(mapped.system, formulation.options)
        n = length(model.label)
        m = length(plan.kept)
        e = length(plan.eliminated)
        nf = length(mapped.frequencies)
        Z = Array{ComplexF64, 3}(undef, m, m, nf)
        Y = Array{ComplexF64, 3}(undef, m, m, nf)
        trace = options.trace ?
                (Z = Array{ComplexF64, 3}(undef, n, n, nf),
            P = Array{ComplexF64, 3}(undef, n, n, nf)) :
                (Z = Array{ComplexF64, 3}(undef, 0, 0, 0),
            P = Array{ComplexF64, 3}(undef, 0, 0, 0))
        buffers = (factor = Matrix{ComplexF64}(undef, e, e),
            coupling = Matrix{ComplexF64}(undef, m, e),
            rhs = Matrix{ComplexF64}(undef, e, m),
            reduced = Matrix{ComplexF64}(undef, m, m),
            identity = Matrix{ComplexF64}(I, m, m))
        solve!(Z, Y, trace, path, model, formulation, buffers, options)
        LineParameters{T}(Array{Complex{T}, 3}(Z), Array{Complex{T}, 3}(Y),
            problem.frequencies, plan.terminals, convert(T, problem.system.line_length),
            formulation, map(array -> Array{Complex{T}, 3}(array), trace))
    finally
        if owned
            gmsh.finalize()
        else
            for name in setdiff(Set(gmsh.model.list()), models)
                gmsh.model.setCurrent(name)
                gmsh.model.remove()
            end
            isempty(current) || gmsh.model.setCurrent(current)
            gmsh.option.setNumber("General.Terminal", terminal)
            gmsh.option.setNumber("General.Verbosity", verbosity)
        end
    end
    if options.keep_run_directory
        @info "FEM run directory" directory=path
    else
        rm(path; recursive = true)
        isempty(readdir(dirname(path))) && rm(dirname(path))
    end
    return result
end

end
