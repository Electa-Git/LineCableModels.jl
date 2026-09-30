module LineCableModelsGmshExt

using Gmsh: gmsh
import Gmsh
import JSON3
import Logging
using Base.BinaryPlatforms: HostPlatform, triplet
using Pkg.Artifacts: artifact_hash, ensure_artifact_installed
using Logging: AbstractLogger, SimpleLogger, @debug, @info,
               @warn, with_logger
using Printf: @sprintf
using SHA: sha256
using LinearAlgebra: I, cond, lu, norm
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import LineCableModels
import LineCableModels: compute
import LineCableModels.DataModel
import LineCableModels.Earth
import LineCableModels.Grammar
import LineCableModels.Engine
import LineCableModels.ImportExport
using LineCableModels.Grammar: computation_options, ComputationOptions, ComputationDetails
using LineCableModels: LineCableModelsFEM, LineCableModelsFEMError,
                       LineParametersProblem, LineParameters,
                       SeriesImpedance, ShuntAdmittance, PhaseDomain

include("model.jl")
include("pml_mesh.jl")
include("formulations.jl")
include("geometry.jl")
include("mesh.jl")
include("getdp.jl")
include("workers.jl")
include("results.jl")
include("compute.jl")
include("import.jl")
include("export.jl")

end
