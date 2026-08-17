"""
    Gauntlet

Provide a typed validation application for comparing native LineCableModels
results with independently harvested PSCAD evidence.
"""
module Gauntlet

using BenchmarkTools
using CodecZlib
using Dates
using DocStringExtensions: FUNCTIONNAME, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using JLD2
using JSON3
using LineCableModels
using LinearAlgebra
using Printf
using SHA
using Statistics
using TOML
using Tables
using Tar

import Base: getindex, iterate, keys, length, show
import LineCableModels: Z, Y, basis, compute!, domain, frequencies

export Approximate,
       AllCases,
       Assumption,
       Case,
       Check,
       Coax,
       Comparison,
       Corpus,
       Diagnostic,
       Exact,
       ExactOnly,
       Family,
       Fail,
       Fidelity,
       Fit,
       FitCheck,
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
       Source,
       Suite,
       Tolerance,
       Terminal,
       Trial,
       Unavailable,
       Verdict,
       assumptions,
       compare,
       evaluate,
       gauntlet,
       ingest,
       performance_baseline,
       provenance,
       reference,
       stage_artifacts,
       stage_smoke,
       write_report

include("types.jl")
include("fits.jl")
include("pscad.jl")
include("materialize.jl")
include("corpus.jl")
include("compare.jl")
include("modal.jl")
include("calibration.jl")
include("fixtures.jl")
include("performance.jl")
include("reports.jl")
include("ingest.jl")
include("artifacts.jl")

end
