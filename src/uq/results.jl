"""
$(TYPEDEF)

Store ordered primitive results from a [`LinearError`](@ref) calculation.

$(TYPEDFIELDS)
"""
struct LinearErrorResult{T, F} <: AbstractUncertaintyResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Uncertainty-bearing primitive results in Gridspace traversal order."
    values::Vector{T}

    function LinearErrorResult(formulation::F, values::Vector{T}) where {T, F}
        ParametricBuilder._primitive_result_type(T)
        return new{T, F}(formulation, values)
    end
end

"""
$(TYPEDEF)

Store primitive sample means and real-valued Monte Carlo summaries.

$(TYPEDFIELDS)
"""
struct MonteCarloResult{T, F, ST <: AbstractVector, S, H} <:
       AbstractUncertaintyResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Primitive results assembled from sample means in Gridspace traversal order."
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

    function MonteCarloResult(
            formulation::F,
            values::Vector{T},
            stats::ST,
            sample_values::S,
            histogram_values::H,
            root_seed::UInt64,
            point_seeds::Vector{UInt64},
            trial_counts::Vector{Int}
    ) where {T, F, ST <: AbstractVector, S, H}
        ParametricBuilder._primitive_result_type(T)
        length(stats) == length(values) || throw(DimensionMismatch(
            "Monte Carlo statistics must contain one entry per primitive result",
        ))
        length(point_seeds) == length(values) || throw(DimensionMismatch(
            "Monte Carlo seeds must contain one entry per primitive result",
        ))
        length(trial_counts) == length(values) || throw(DimensionMismatch(
            "Monte Carlo trial counts must contain one entry per primitive result",
        ))
        sample_values === nothing || length(sample_values) == length(values) ||
            throw(DimensionMismatch(
                "retained samples must contain one entry per primitive result",
            ))
        histogram_values === nothing || length(histogram_values) == length(values) ||
            throw(DimensionMismatch(
                "retained histograms must contain one entry per primitive result",
            ))
        return new{T, F, ST, S, H}(
            formulation,
            values,
            stats,
            sample_values,
            histogram_values,
            root_seed,
            point_seeds,
            trial_counts
        )
    end
end
