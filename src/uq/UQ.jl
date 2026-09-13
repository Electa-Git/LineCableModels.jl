"""
    LineCableModels.UQ

Define direct linear propagation, conditional Monte Carlo sampling, retained
statistics, and uncertainty-result presentation.
"""
module UQ
import ..LineCableModels: description, formula_id
import ..Grammar: formulation_options

export LinearError, MonteCarlo, LinearErrorResult, MonteCarloResult
export SampleSummary, HistogramDensity
export statistics, samples, histograms, uncertain
export root_seed, point_seed, trial_count
export confidence, cdf_tolerance, sampling_distribution
export cumulative_probability, quantile_pairs

import Random
import Statistics
import ..LineCableModels
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: basis, frequencies, R, L, C, nominal, uncertainty
import ..LineCableModels: points, realize, realize_arguments, Gridpoint
import ..LineCableModels: progress_receiver, report_progress, with_progress_scope, with_scan_progress
import ..DataModel
import ..Engine
import ..Grammar: compute, computation_options, computation_details, details,
                  observe, observables, check_core_result,
                  detach, publication_table, request_identity, request_indices,
                  observation_indices, observation_resolution, observation_request
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
include("comparisons.jl")
include("linearerror.jl")
include("montecarlo/compute.jl")
include("textdisplay.jl")

end
