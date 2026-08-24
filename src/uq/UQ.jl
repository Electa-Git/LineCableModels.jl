"""
    LineCableModels.UQ

Define direct linear propagation, conditional Monte Carlo sampling, retained
statistics, and uncertainty-result presentation.
"""
module UQ

export LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity, RLCG
export result, statistics, samples, histograms, uncertain_value

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
import ..ParametricBuilder: result
import ..PlotBuilder
import ..QuantityUnits
using ..Grammar:
                 AbstractFormulation, AbstractProblemResult, AbstractUncertaintyResult
using ..ParametricBuilder:
                           ParametricProblem

include("formulations.jl")
include("statistics.jl")
include("results.jl")
include("base.jl")
include("observables.jl")
include("linearerror.jl")
include("montecarlo/compute.jl")
include("dataframe.jl")
include("montecarlo/plot.jl")

end
