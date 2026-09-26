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
export select, tabulate, illustrate, encode, write

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import Statistics
import DataFrames
import DataFrames: DataFrame, metadata, metadata!, Not
import ..Grammar: observables
import ..Grammar: ObservedResult
import ..Grammar
import ..Units
import ..DataModel
import ..Engine
import ..Grammar: AbstractUncertaintyResult, request_identity, request_quantity, request_indices
import ..LineCableModels: validate, description
import ..LineCableModels
import ..PlotBuilder
import ..UQ
import ..TextDisplay
import ..LineCableModels: Z, Y, R, X, L, G, B, C

include("grammar.jl")
include("tables.jl")
include("comparisons.jl")
include("montecarlo.jl")
include("performance.jl")
include("xlsx.jl")
include("textdisplay.jl")

public observation_columns, encode_cell, XLSXSheet, XLSXWorkbook

end
