module PSCADBenchmarks

using Base64: base64encode
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.EarthProps: EarthModel
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description
import LineCableModels.Engine: AbstractAdmittanceFormulation,
                               AbstractFormulation,
                               AbstractFormulationOptions,
                               AbstractImpedanceFormulation,
                               ComputeOptions,
                               EarthAdmittanceFormulation,
                               EarthImpedanceFormulation,
                               Formulation,
                               InsulationAdmittanceFormulation,
                               LineParametersProblem,
                               ProblemDefinition,
                               compute!, verbosity
import ..GauntletSupport: GAUNTLET_ROOT, WORK_ROOT,
                          benchmark_metadata, formulation_record

export NativeEarthAdmittance, NativeInsulationAdmittance,
       PSCADFormulation, PSCADOptions, RemoteConfig,
       read_pscad_result, remote_command, run_remote_pscad

include("formulations.jl")
include("outputs.jl")
include("remote/files.jl")
include("remote.jl")

end
