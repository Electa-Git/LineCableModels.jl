"Native expert controls, each with its default and allowed range."
const OVERRIDES = (
    PmlSideThicknessFactor = (default = 1.0, range = Interval{:open, :open}(0.0, Inf)),
    PmlTopThicknessFactor = (default = 1.0, range = Interval{:open, :open}(0.0, Inf)),
    PmlBottomThicknessFactor = (default = 1.0, range = Interval{:open, :open}(0.0, Inf)),
    PmlSideLayers = (default = 16, range = 1:(typemax(Cint) - 1)),
    PmlTopLayers = (default = 16, range = 1:(typemax(Cint) - 1)),
    PmlBottomLayers = (default = 16, range = 1:(typemax(Cint) - 1)),
    PmlSideGrading = (default = log(20), range = Interval{:closed, :open}(0.0, Inf)),
    PmlTopGrading = (default = log(20), range = Interval{:closed, :open}(0.0, Inf)),
    PmlBottomGrading = (default = log(20), range = Interval{:closed, :open}(0.0, Inf)),
    ExteriorMeshSizeFactor = (default = 1.0, range = Interval{:closed, :open}(1.0, Inf)),
    InterfaceRefinementFactor = (default = 1.0, range = Interval{:closed, :open}(1.0, Inf)),
    ConductorGeometryTolerance = (default = 1e-3, range = Interval{:open, :open}(0.0, Inf)),
    ConductorSkinDepthElements = (default = 6.0, range = Interval{:open, :open}(0.0, Inf)),
    ConductorMeshGrowth = (
        default = sqrt(1.25), range = Interval{:closed, :open}(1.0, Inf)),
    ConductorSkinDepths = (default = 5.0, range = Interval{:open, :open}(0.0, Inf)),
    ConductorThicknessElements = (default = 4, range = 1:typemax(Cint)),
    VolumeQuadrature = (default = 12, range = (4, 7, 12, 13)),
    PhysicalVolumeQuadrature = (default = 3, range = (3, 4, 7, 12, 13)),
    PmlQuadrangles = (default = 1, range = (0, 1)),
    PmlQuadrature = (default = 4, range = (4, 9, 16)),
    MumpsOrdering = (default = 0, range = (-1, 0, 2, 3, 4, 5, 6, 7)),
    PetscPrealloc = (default = 0, range = 0:typemax(Cint)),
    MumpsErrorAnalysis = (default = 2, range = (0, 1, 2)),
    MumpsRefinementMax = (default = 2, range = 0:typemax(Cint)),
    MumpsBackwardErrorTolerance = (
        default = 1e-12, range = Interval{:open, :open}(0.0, Inf)),
    ThinLayerStrips = (default = 1, range = (0, 1)))

"""
$(TYPEDSIGNATURES)

Return the segments per turn whose inscribed polygon bounds the relative error
of a circle's area by `τ` \\[dimensionless\\]:

``n=12⋅2^{max(0,⌈log₂(√(2π²/(3τ))/12)⌉)}``.

The count is an `Int`.
"""
function segments(τ::Real)
    return 12 * 2^max(0, ceil(Int, log2(sqrt(2π^2 / (3τ)) / 12)))
end

"""
$(TYPEDSIGNATURES)

Check the Gmsh backend's factors, process counts, verbosity levels, and native controls.

# Arguments

- `options`: computation options of the Gmsh backend.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a nonfinite or nonpositive `mesh_size_factor`.
- Throws `ArgumentError` for a nonfinite or nonpositive `domain_size_factor`.
- Throws `ArgumentError` for `pml_reflection` outside the open interval (0, 1).
- Throws `ArgumentError` for `workers` less than 1.
- Throws `ArgumentError` for `solver_threads` less than 1.
- Throws `ArgumentError` for `gmsh_verbosity` outside `0:5`.
- Throws `ArgumentError` for `getdp_verbosity` outside `0:5`.
- Throws `ArgumentError` for a native control in `overrides` outside its range in
  [`OVERRIDES`](@ref).
"""
function validate(options::ComputationOptions{GmshBackend})
    isfinite(options.mesh_size_factor) && options.mesh_size_factor > 0 ||
        throw(ArgumentError(
            "mesh_size_factor must be finite and positive; received $(options.mesh_size_factor)"))
    isfinite(options.domain_size_factor) && options.domain_size_factor > 0 ||
        throw(ArgumentError(
            "domain_size_factor must be finite and positive; received $(options.domain_size_factor)"))
    0 < options.pml_reflection < 1 || throw(ArgumentError(
        "pml_reflection must be strictly between 0 and 1; received $(options.pml_reflection)"))
    options.workers >= 1 || throw(ArgumentError(
        "workers must be at least 1; received $(options.workers)"))
    options.solver_threads >= 1 || throw(ArgumentError(
        "solver_threads must be at least 1; received $(options.solver_threads)"))
    options.gmsh_verbosity in 0:5 || throw(ArgumentError(
        "gmsh_verbosity must be in 0:5; received $(options.gmsh_verbosity)"))
    options.getdp_verbosity in 0:5 || throw(ArgumentError(
        "getdp_verbosity must be in 0:5; received $(options.getdp_verbosity)"))
    for (name, entry) in pairs(OVERRIDES)
        value = options.overrides[name]
        value in entry.range || throw(ArgumentError(
            "overrides.$name must be in $(entry.range); received $value"))
    end
    return options
end
