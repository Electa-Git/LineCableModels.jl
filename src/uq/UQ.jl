"""
    LineCableModels.UQ

Define direct linear propagation, conditional Monte Carlo sampling, retained
statistics, and uncertainty-result presentation.
"""
module UQ

export LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity
export result, statistics, samples, histograms, uncertain_value
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution
export cumulative_probability, quantile_pairs

import Random
import Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, domain, frequencies, R, L, C
import ..DataModel
import ..Engine
import ..Grammar: compute, computation_options, computation_details, details,
                  observe, observables, check_core_result, computation_owner,
                  detach
import ..ParametricBuilder
import ..ParametricBuilder: result, traverse
import ..Units
using ..Grammar:
                 AbstractFormulation, AbstractUncertaintyResult,
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

end
