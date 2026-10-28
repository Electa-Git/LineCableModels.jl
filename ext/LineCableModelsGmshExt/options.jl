"""
$(TYPEDEF)

Represent the GetDP 3.5.0 executable of the package's artifact.
"""
struct GetdpArtifact end

select(::GetdpArtifact, path::AbstractString) = String(path)

"""
$(TYPEDEF)

Represent a real interval for a native control. The low bound is inside when `closed`
is true. The high bound is always outside.

$(TYPEDFIELDS)
"""
struct Interval{T <: Real}
    "Low bound."
    low::T
    "High bound."
    high::T
    "Whether the low bound is inside."
    closed::Bool
    function Interval(low::T, high::T, closed::Bool) where {T <: Real}
        return new{T}(low, high, closed)
    end
end

function Base.in(x::Real, interval::Interval)
    return (interval.closed ? interval.low <= x : interval.low < x) && x < interval.high
end

function Base.show(io::IO, interval::Interval)
    print(io, interval.closed ? '[' : '(',
        interval.low, ", ",
        isinf(interval.high) ? "∞" : interval.high, ')')
end

"Native expert controls, each with its default and allowed range."
const OVERRIDES = (
    PmlSideThicknessFactor = (default = 1.0, range = Interval(0.0, Inf, false)),
    PmlTopThicknessFactor = (default = 1.0, range = Interval(0.0, Inf, false)),
    PmlBottomThicknessFactor = (default = 1.0, range = Interval(0.0, Inf, false)),
    PmlSideLayers = (default = 16, range = 1:(typemax(Cint) - 1)),
    PmlTopLayers = (default = 16, range = 1:(typemax(Cint) - 1)),
    PmlBottomLayers = (default = 16, range = 1:(typemax(Cint) - 1)),
    PmlSideGrading = (default = log(20), range = Interval(0.0, Inf, true)),
    PmlTopGrading = (default = log(20), range = Interval(0.0, Inf, true)),
    PmlBottomGrading = (default = log(20), range = Interval(0.0, Inf, true)),
    ExteriorMeshSizeFactor = (default = 1.0, range = Interval(1.0, Inf, true)),
    InterfaceRefinementFactor = (default = 1.0, range = Interval(1.0, Inf, true)),
    ConductorGeometryTolerance = (default = 1e-3, range = Interval(0.0, Inf, false)),
    ConductorSkinDepthElements = (default = 6.0, range = Interval(0.0, Inf, false)),
    ConductorMeshGrowth = (default = sqrt(1.25), range = Interval(1.0, Inf, true)),
    ConductorSkinDepths = (default = 5.0, range = Interval(0.0, Inf, false)),
    ConductorThicknessElements = (default = 4, range = 1:typemax(Cint)),
    VolumeQuadrature = (default = 12, range = (4, 7, 12, 13)),
    PhysicalVolumeQuadrature = (default = 3, range = (3, 4, 7, 12, 13)),
    PmlQuadrangles = (default = 1, range = (0, 1)),
    PmlQuadrature = (default = 4, range = (4, 9, 16)),
    MumpsOrdering = (default = 0, range = (-1, 0, 2, 3, 4, 5, 6, 7)),
    PetscPrealloc = (default = 0, range = 0:typemax(Cint)),
    MumpsErrorAnalysis = (default = 2, range = (0, 1, 2)),
    MumpsRefinementMax = (default = 2, range = 0:typemax(Cint)),
    MumpsBackwardErrorTolerance = (default = 1e-12, range = Interval(0.0, Inf, false)),
    ThinLayerStrips = (default = 1, range = (0, 1)))

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
