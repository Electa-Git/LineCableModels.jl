"""
    LineCableModels.UQ

Define direct linear propagation, conditional Monte Carlo sampling, retained
statistics, and uncertainty-result presentation.
"""
module UQ

export LinearError, MonteCarlo, Sensitivity
export LinearErrorResult, MonteCarloResult, SensitivityResult
export SampleSummary, HistogramDensity
export result, statistics, samples, histograms, uncertain
export root_seed, point_seed, trial_count
export first_order, total_order, second_order
export confidence, cdf_tolerance, sampling_distribution
export cumulative_probability, quantile_pairs

import Random
import Statistics
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, frequencies, R, L, C
import ..DataModel
import ..Engine
import ..Grammar: compute, computation_options, computation_details, details,
                  observe, observables, check_core_result, computation_owner,
                  detach, publication_table, request_identity, request_indices,
                  observation_indices
import ..ParametricBuilder
import ..ParametricBuilder: result, traverse
import ..Units
import ..TextDisplay
using ..Grammar:
                 AbstractFormulation, AbstractProblemResult,
                 AbstractUncertaintyResult,
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
