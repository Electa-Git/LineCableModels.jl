module LineCableModelsGmshExt

using Gmsh
using LineCableModels

const _GETDP_DEPRECATION_MESSAGE = "The bundled GetDP frontend is a temporary compatibility layer and will be retired with the legacy FEM interface."

function _warn_getdp_frontend()
    Base.depwarn(_GETDP_DEPRECATION_MESSAGE, :GetDPFrontend)
    return nothing
end

module FEMImplementation

using LineCableModels.Commons
using LineCableModels.Materials
using LineCableModels.EarthProps
using LineCableModels.DataModel
using LineCableModels.Engine
import LineCableModels.Engine: kronify, reorder_M, reorder_indices, merge_bundles!,
                               AbstractFormulationSet, AbstractImpedanceFormulation,
                               AbstractAdmittanceFormulation,
                               compute!
import LineCableModels.Commons: PhaseDomain, ModalDomain, LineParamsDomain, domain
import LineCableModels.Engine: AbstractFormulationOptions, LineParamOptions, build_options,
                               _COMMON_SYMS
import LineCableModels.DataModel: AbstractCablePart, AbstractConductorPart,
                                  AbstractInsulatorPart
using LineCableModels.Utils: display_path, set_verbosity!, is_headless, to_nominal,
                             symtrans!, symtrans, line_transpose!
using Measurements
using LinearAlgebra
using Colors
using GeometryBasics: Point, Point2f
using Gmsh

include("fem/getdp_frontend/GetDPFrontend.jl")
using .GetDPFrontend
using .GetDPFrontend: Problem, get_getdp_executable, add!
const GetDP = GetDPFrontend

const FEM_SOURCE = joinpath(@__DIR__, "..", "src", "engine", "fem")

include(joinpath(FEM_SOURCE, "types.jl"))
include(joinpath(FEM_SOURCE, "lineparamopts.jl"))
include(joinpath(FEM_SOURCE, "meshtransitions.jl"))
include(joinpath(FEM_SOURCE, "problemdefs.jl"))
include(joinpath(FEM_SOURCE, "workspace.jl"))
include(joinpath(FEM_SOURCE, "encoding.jl"))
include(joinpath(FEM_SOURCE, "drawing.jl"))
include(joinpath(FEM_SOURCE, "identification.jl"))
include(joinpath(FEM_SOURCE, "mesh.jl"))
include(joinpath(FEM_SOURCE, "materialprops.jl"))
include(joinpath(FEM_SOURCE, "helpers.jl"))
include(joinpath(FEM_SOURCE, "visualization.jl"))
include(joinpath(FEM_SOURCE, "space.jl"))
include(joinpath(FEM_SOURCE, "cable.jl"))
include(joinpath(FEM_SOURCE, "solver.jl"))
include(joinpath(FEM_SOURCE, "base.jl"))

end # module FEMImplementation

function Darwin(args...; kwargs...)
    _warn_getdp_frontend()
    return FEMImplementation.Darwin(args...; kwargs...)
end

function Electrodynamics(args...; kwargs...)
    _warn_getdp_frontend()
    return FEMImplementation.Electrodynamics(args...; kwargs...)
end

MeshTransition(args...; kwargs...) = FEMImplementation.MeshTransition(args...; kwargs...)
function calc_domain_size(args...; kwargs...)
    FEMImplementation.calc_domain_size(args...; kwargs...)
end
preview_results(args...; kwargs...) = FEMImplementation.preview_results(args...; kwargs...)

function formulation_set(; kwargs...)
    _warn_getdp_frontend()
    return FEMImplementation.formulation_set(; kwargs...)
end

end # module LineCableModelsGmshExt
