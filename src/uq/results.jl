"""
$(TYPEDEF)

Store ordered core results from a [`LinearError`](@ref) calculation.

$(TYPEDFIELDS)
"""
struct LinearErrorResult{T, F, D <: ComputationDetails} <: AbstractUncertaintyResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Uncertainty-bearing core results in Gridspace traversal order."
    values::Vector{T}
    "Typed supplemental output retained by the propagation."
    details::D

    function LinearErrorResult(
            formulation::F,
            values::Vector{T},
            details::D
    ) where {T, F, D <: ComputationDetails}
        check_core_result(T)
        isempty(details) || keys(details) == (:points,) ||
            throw(ArgumentError(
                "LinearErrorResult details must be empty or contain only points",
            ))
        isempty(details) || length(details.points) == length(values) ||
            throw(DimensionMismatch(
                "retained details must contain one entry per core result",
            ))
        return new{T, F, D}(formulation, values, details)
    end
end

LinearErrorResult(formulation, values) = LinearErrorResult(formulation, values, (;))

"""
$(TYPEDEF)

Store core results reconstructed from sample means and Monte Carlo summaries.

$(TYPEDFIELDS)
"""
struct MonteCarloResult{T, F, ST <: AbstractVector, S, H, D <: ComputationDetails} <:
       AbstractUncertaintyResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Core results assembled from sample means in Gridspace traversal order."
    values::Vector{T}
    "Per-observable sample summaries."
    stats::ST
    "Retained joint samples, or `nothing` when not requested."
    sample_values::S
    "Retained marginal histograms, or `nothing` when not requested."
    histogram_values::H
    "Resolved root random seed."
    root_seed::UInt64
    "Deterministic random seed for each Gridspace point."
    point_seeds::Vector{UInt64}
    "Trial count used for each Gridspace point."
    trial_counts::Vector{Int}
    "Typed supplemental output retained by the propagation."
    details::D

    function MonteCarloResult(
            formulation::F,
            values::Vector{T},
            stats::ST,
            sample_values::S,
            histogram_values::H,
            root_seed::UInt64,
            point_seeds::Vector{UInt64},
            trial_counts::Vector{Int},
            details::D
    ) where {T, F, ST <: AbstractVector, S, H, D <: ComputationDetails}
        check_core_result(T)
        length(stats) == length(values) || throw(DimensionMismatch(
            "Monte Carlo statistics must contain one entry per core result",
        ))
        length(point_seeds) == length(values) || throw(DimensionMismatch(
            "Monte Carlo seeds must contain one entry per core result",
        ))
        length(trial_counts) == length(values) || throw(DimensionMismatch(
            "Monte Carlo trial counts must contain one entry per core result",
        ))
        sample_values === nothing || length(sample_values) == length(values) ||
            throw(DimensionMismatch(
                "retained samples must contain one entry per core result",
            ))
        histogram_values === nothing || length(histogram_values) == length(values) ||
            throw(DimensionMismatch(
                "retained histograms must contain one entry per core result",
            ))
        isempty(details) || keys(details) == (:trials,) ||
            throw(ArgumentError(
                "MonteCarloResult details must be empty or contain only trials",
            ))
        if !isempty(details)
            length(details.trials) == length(values) || throw(DimensionMismatch(
                "retained details must contain one entry per Gridspace point",
            ))
            all(length(records) == count
            for (records, count) in zip(details.trials, trial_counts)) ||
                throw(DimensionMismatch(
                    "retained details must contain one entry per Monte Carlo trial",
                ))
        end
        return new{T, F, ST, S, H, D}(
            formulation,
            values,
            stats,
            sample_values,
            histogram_values,
            root_seed,
            point_seeds,
            trial_counts,
            details
        )
    end
end

function MonteCarloResult(
        formulation,
        values,
        stats,
        sample_values,
        histogram_values,
        root_seed,
        point_seeds,
        trial_counts
)
    return MonteCarloResult(
        formulation,
        values,
        stats,
        sample_values,
        histogram_values,
        root_seed,
        point_seeds,
        trial_counts,
        (;)
    )
end

details(value::LinearErrorResult) = value.details
details(value::MonteCarloResult) = value.details
