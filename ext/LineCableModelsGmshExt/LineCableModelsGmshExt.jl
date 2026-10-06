module LineCableModelsGmshExt

using Gmsh: gmsh
import Gmsh
import JSON3
import Logging
using Base.BinaryPlatforms: HostPlatform, triplet
using Pkg.Artifacts: artifact_hash, ensure_artifact_installed
using Logging: AbstractLogger, SimpleLogger, @debug, @info,
               @warn, with_logger
using Dates: Dates, now
using Random: RandomDevice, rand
using Printf: @sprintf
using SHA: sha256
using LinearAlgebra: I, cond, lu, norm
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import LineCableModels
import LineCableModels: compute, Formulation, description, formula_id
using LineCableModels: formula, parameterize
import LineCableModels.TextDisplay
import LineCableModels.DataModel
import LineCableModels.Earth
import LineCableModels.Grammar
import LineCableModels.Engine
import LineCableModels.ImportExport
import LineCableModels.Grammar: computation_options, formulation_options, computation_details
using LineCableModels.Grammar: AbstractFormulation, FormulationOptions, ComputationOptions, ComputationDetails, details
using LineCableModels.Engine: LineParametersFormulation, LineCableModelsCoaxial, InsulationAdmittance, SemiconAdmittance
using LineCableModels.Materials: TemperatureDependent
using LineCableModels: LineParametersProblem, LineParameters,
                       SeriesImpedance, ShuntAdmittance, PhaseDomain

public LineCableModelsFEM, LineCableModelsFEMError, FEMMesh, FEMElementBlock, FEMFieldMap, FEMFieldBlock

include("formulations.jl")
include("options.jl")
include("data.jl")
include("textdisplay.jl")
include("model.jl")
include("geometry.jl")
include("mesh.jl")
include("getdp.jl")
include("workers.jl")
include("results.jl")
include("compute.jl")
include("import.jl")
include("export.jl")

end
