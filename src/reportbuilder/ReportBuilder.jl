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
using Printf: @sprintf
import Statistics
import DataFrames
using RequiredInterfaces: @required
import DataFrames: DataFrame, metadata, metadata!, select!, Not
import ..Grammar: observables
import ..Grammar: ObservationPublication
using ..Grammar: @observe
import ..Units
import ..DataModel
import ..Engine
import ..Engine: domain
import ..ParametricBuilder: ParametricResult
import ..Grammar: AbstractCoreResult, AbstractUncertaintyResult, details, observe, detach,
                  request_identity, request_quantity, request_indices, observation_request,
                  materialize_observation
import ..LineCableModels: validate, description, formula_id
import ..LineCableModels
import ..PlotBuilder
import ..UQ
import ..TextDisplay
import ..LineCableModels: basis, Z, Y, R, X, L, G, B, C

include("grammar.jl")
include("tables.jl")
include("comparisons.jl")
include("montecarlo.jl")
include("performance.jl")
include("xlsx.jl")
include("textdisplay.jl")

public observation_columns, encode_cell, XLSXSheet, XLSXWorkbook

end
