"""
    LineCableModels.ReportBuilder

Build human-facing tables and optional plot artifacts from published scientific
observations.
"""
module ReportBuilder

export AbstractReportDefinition, ReportArtifact
export TableReportDefinition, CableConstantsTableDefinition
export LineParametersTableDefinition, BenchmarkTableDefinition
export MonteCarloTableDefinition, XLSXReportDefinition, report
export entitle, select, tabulate, illustrate, encode, write, finish

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
using LinearAlgebra
using Printf
import DataFrames: DataFrame, metadata, metadatakeys, metadata!
import ..Grammar: observables, validate_observables, unit_targets
import ..Units
import ..PlotBuilder
import ..DataModel
import ..Engine
import ..Engine: line_requests, line_parent
import ..UQ
import ..LineCableModels: basis, frequencies, Z, Y, R, X, L, G, B, C

include("grammar.jl")
include("tables.jl")
include("montecarlo.jl")
include("xlsx.jl")

public clip, encode_cell

end
