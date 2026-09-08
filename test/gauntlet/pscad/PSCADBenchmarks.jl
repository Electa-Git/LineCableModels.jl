module PSCADBenchmarks

using Base64: base64encode
using SHA: sha256
import TOML
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.Earth: EarthModel
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description, parameterize, computation_details
import LineCableModels.Engine: AbstractAdmittanceFormulation,
                               AbstractImpedanceFormulation,
                               EarthAdmittanceFormulation,
                               EarthImpedanceFormulation,
                               Formulation,
                               InsulationAdmittanceFormulation,
                               LineParametersProblem,
                               verbosity
import LineCableModels.Grammar: AbstractFormulation, ComputationOptions,
                                FormulationOptions, computation_options, compute,
                                formulation_options
import ..GauntletSupport: GAUNTLET_ROOT, WORK_ROOT,
                          benchmark_metadata, formulation_record, semantic_sha256

export PSCADFormulation, RemoteConfig,
       read_pscad_result, remote_command, run_remote_pscad, formulas

include("formulations.jl")
include("outputs.jl")
include("remote/files.jl")
include("remote.jl")

end
