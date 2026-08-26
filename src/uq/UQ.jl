"""
    LineCableModels.UQ

Define direct linear propagation, conditional Monte Carlo sampling, retained
statistics, and uncertainty-result presentation.
"""
module UQ

export LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity, RLCG
export result, statistics, samples, histograms, uncertain_value
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution

using Random
using Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: PhaseDomain, basis, domain, frequencies, R, L, C
import ..DataModel
import ..Engine
import ..Grammar
import ..Grammar: compute, computation_options, computation_details, details,
                  observe, observables, check_core_result
import ..ParametricBuilder
import ..ParametricBuilder: result
import ..PlotBuilder
import ..Units
using ..Grammar:
                 AbstractFormulation, AbstractProblemResult, AbstractResultSpace,
                 AbstractUncertaintyResult,
                 ComputationOptions, ComputationDetails
using ..ParametricBuilder:
                           ParametricProblem

include("formulations.jl")
include("statistics.jl")
include("results.jl")
include("base.jl")
include("observations.jl")
include("linearerror.jl")
include("montecarlo/compute.jl")
include("montecarlo/plot.jl")

end
