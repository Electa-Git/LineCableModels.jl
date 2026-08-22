module PSCADHarness

using Base64: base64encode
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.EarthProps: EarthModel
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description
import ..GauntletSupport: GAUNTLET_ROOT, GauntletCase,
                          benchmark_local, comparison_passes,
                          load_prior_snapshot, performance_comparison,
                          persist_snapshot, validate_structure, work_path

export RemoteConfig, PSCAD_FORMULATIONS, formulation_spec,
       read_pscad_result, remote_command, run_live, run_remote_pscad

include("formulations.jl")
include("outputs.jl")
include("remote.jl")

end
