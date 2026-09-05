module LineCableModelsGmshExt

using Gmsh: gmsh
import Gmsh
import JSON3
import Logging
using Logging: AbstractLogger, ConsoleLogger, SimpleLogger, @debug, @info,
               @warn, with_logger
using Printf: @sprintf
using SHA: sha256

import LineCableModels
import LineCableModels: compute
import LineCableModels.DataModel
import LineCableModels.EarthProps
import LineCableModels.Engine
import LineCableModels.ImportExport
using LineCableModels.Grammar: computation_options
using LineCableModels: LineCableModelsFEM, LineCableModelsFEMError,
                       LineParametersProblem, LineParameters,
                       SeriesImpedance, ShuntAdmittance, PhaseDomain

include("model.jl")
include("geometry.jl")
include("mesh.jl")
include("onelab.jl")
include("getdp.jl")
include("results.jl")
include("compute.jl")

end
