module PSCADBenchmarks

using Base64: base64encode
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.EarthProps: EarthModel
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description
import LineCableModels.Engine: AbstractAdmittanceFormulation,
                               AbstractImpedanceFormulation,
                               EarthAdmittanceFormulation,
                               EarthImpedanceFormulation,
                               Formulation,
                               InsulationAdmittanceFormulation,
                               LineParametersProblem,
                               computation_options, verbosity
import LineCableModels.Grammar: AbstractFormulation, compute
import ..GauntletSupport: GAUNTLET_ROOT, WORK_ROOT,
                          benchmark_metadata, formulation_record

export NativeEarthAdmittance, NativeInsulationAdmittance,
       PSCADFormulation, RemoteConfig,
       read_pscad_result, remote_command, run_remote_pscad

include("formulations.jl")
include("outputs.jl")
include("remote/files.jl")
include("remote.jl")

end
