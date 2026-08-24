"""
    LineCableModels.UQ

Own uncertainty formulations, evaluation, retained statistical products, and
uncertainty-specific presentation.
"""
module UQ

export LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity, RLCG
export result, statistics, samples, histograms, uncertain_value, manifest

using Random
using Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import DataFrames: DataFrame, metadata!
import ..LineCableModels: PhaseDomain, basis, domain, frequencies
import ..DataModel
import ..Engine
import ..Grammar
import ..Grammar: compute, observables
import ..ParametricBuilder
import ..ParametricBuilder: result, manifest
import ..PlotBuilder
import ..UnitHandler
using ..Grammar:
                 AbstractFormulation, AbstractProblemResult, AbstractUncertaintyResult
using ..ParametricBuilder:
                          ParametricProblem, ConfigurationFailure,
                          configurations, configuration_manifest, materialize

include("types.jl")
include("compute.jl")
include("dataframe.jl")
include("plotspecs.jl")

end
