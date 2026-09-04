"""
    LineCableModels.UQ

Define direct linear propagation, conditional Monte Carlo sampling, retained
statistics, and uncertainty-result presentation.
"""
module UQ

export LinearError, MonteCarlo, PolynomialChaos
export LinearErrorResult, MonteCarloResult, PolynomialChaosResult
export SampleSummary, HistogramDensity
export statistics, samples, histograms, expansions, validation, uncertain
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution
export cumulative_probability, quantile_pairs

import Random
import Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, frequencies, R, L, C
import ..LineCableModels: points, realize, realize_arguments
import ..DataModel
import ..Engine
import ..Grammar: compute, computation_options, computation_details, details,
                  observe, observables, check_core_result,
                  detach, publication_table, request_identity, request_indices,
                  observation_indices
import ..ParametricBuilder
import ..ParametricBuilder: traverse
import ..Units
import ..TextDisplay
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
include("publication.jl")
include("linearerror.jl")
include("montecarlo/compute.jl")
include("textdisplay.jl")

end
