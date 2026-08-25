"""
    LineCableModels.ReportBuilder

Build human-facing tables and optional plot artifacts from published scientific
observations.
"""
module ReportBuilder

export AbstractReportDefinition, ReportArtifact, TableReport, XLSXReport, report
export entitle, select, tabulate, illustrate, encode, write, finish

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using LinearAlgebra
using Printf
using Tables
using XLSX
import DataFrames: DataFrame, metadata, metadatakeys, metadata!
import ..Grammar: observables
import ..Units
import ..PlotBuilder
import ..DataModel
import ..Engine
import ..UQ
import ..LineCableModels: basis, frequencies, Z, Y, R, X, L, G, B, C

include("grammar.jl")
include("tables.jl")
include("montecarlo.jl")
include("xlsx.jl")

end
