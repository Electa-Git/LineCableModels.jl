"""
    Gauntlet

Provide a typed, datasource-agnostic validation application for comparing
native LineCableModels results with independently produced evidence.
"""
module Gauntlet

using BenchmarkTools
using CodecZlib
using Dates
using Downloads
using DocStringExtensions: FUNCTIONNAME, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using EzXML
import JLD2
using JSON3
using LineCableModels
using LinearAlgebra
using Printf
using SHA
using Statistics
using TOML
using Tables
using Tar
using ZipFile

import Base: getindex, iterate, keys, length, show
import LineCableModels: Z, Y, basis, compute!, domain, frequencies

const DM = LineCableModels.DataModel
const EP = LineCableModels.EarthProps
const EN = LineCableModels.Engine

export Approximate,
       AllCases,
       Assumption,
       Case,
       Check,
       Coax,
       Comparison,
       Datasource,
       Dataset,
       Deferred,
       Diagnostic,
       Exact,
       ExactOnly,
       Family,
       Fail,
       Fidelity,
       Fit,
       FitCheck,
       FEM,
       MatrixCheck,
       Metrics,
       Mixed,
       ModalCheck,
       Modes,
       NoReduction,
       Overhead,
       PSCAD,
       Pass,
       Perf,
       PerformanceCheck,
       PhysicalCheck,
       Pipe,
       Port,
       Provenance,
       Policy,
       Reduced,
       Reference,
       ReferenceRejected,
       Reduction,
       Rejected,
       Report,
       Retained,
       Suite,
       Tolerance,
       Terminal,
       Trial,
       Unavailable,
       Verdict,
       assumptions,
       compare,
       datasource,
       decode,
       evaluate,
       gauntlet,
       ingest,
       load,
       performance_baseline,
       provenance,
       pending,
       reference,
       stage_artifacts,
       stage_smoke,
       verify,
       write_report

include("types.jl")
include("fits.jl")
include("dataset.jl")
include("compare.jl")
include("modal.jl")
include("calibration.jl")
include("fixtures.jl")
include("performance.jl")
include("reports.jl")
include("artifacts.jl")
include("pscad/io.jl")
include("pscad/materialize.jl")
include("pscad/load.jl")
include("pscad/ingest.jl")
include("pscad/calibration.jl")
include("pscad/artifacts.jl")
include("pscad/harvest.jl")
include("pscad/harvest_cli.jl")

function (@main)(arguments)
    exit(harvest_main(arguments))
end

end
