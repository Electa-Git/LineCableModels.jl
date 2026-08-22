module PSCADHarness

using Base64: base64encode
using BenchmarkTools
using Dates: datetime2unix, now
using EzXML
using Statistics
using LineCableModels
using LineCableModels.DataModel: LineCableSystem
using LineCableModels.EarthProps: EarthModel
using LineCableModels.Engine
using LineCableModels.ImportExport
import LineCableModels: description
import ..GauntletSupport: GAUNTLET_ROOT, GauntletCase,
                          snapshot_path, validate_structure,
                          work_path, write_snapshot

export RemoteConfig, PSCAD_FORMULATIONS, formulation_spec,
       read_pscad_result, remote_command, run_live, run_remote_pscad

include("formulations.jl")
include("outputs.jl")
include("../benchmark.jl")
include("remote.jl")
include("legacy_export.jl")

end
