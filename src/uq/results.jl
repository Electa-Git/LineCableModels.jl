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

_monte_carlo_keys(::Type{<:DataModel.CableConstants}) = (:R, :L, :C)
_monte_carlo_keys(::Type{<:Engine.LineParameters}) = (:R, :L, :C, :G)

function _validated_product(product, expected_keys::Tuple, name::AbstractString)
    product isa NamedTuple || throw(ArgumentError(
        "Monte Carlo $name must be a named tuple",
    ))
    keys(product) == expected_keys || throw(ArgumentError(
        "Monte Carlo $name must contain exactly $(join(expected_keys, ", "))",
    ))
    return product
end

function _validate_cable_products(
        statistics_product,
        sample_product,
        histogram_product,
        trials::Int
)
    all(summary -> summary isa SampleSummary && summary.n == trials,
        Base.values(statistics_product)) || throw(DimensionMismatch(
        "cable-constant summaries must contain the Monte Carlo trial count",
    ))
    if sample_product !== nothing
        all(sample -> sample isa AbstractVector && length(sample) == trials,
            Base.values(sample_product)) || throw(DimensionMismatch(
            "cable-constant samples must share the Monte Carlo trial dimension",
        ))
    end
    if histogram_product !== nothing
        all(histogram -> histogram isa HistogramDensity,
            Base.values(histogram_product)) || throw(ArgumentError(
            "cable-constant histograms must contain HistogramDensity values",
        ))
    end
    return nothing
end

function _line_product_shape(product, expected_shape, name::AbstractString)
    all(field -> field isa AbstractArray && size(field) == expected_shape,
        Base.values(product)) || throw(DimensionMismatch(
        "line-parameter $name must match the core matrix and frequency dimensions",
    ))
    return nothing
end

function _validate_line_products(
        value::Engine.LineParameters,
        statistics_product,
        sample_product,
        histogram_product,
        trials::Int
)
    core_shape = size(observe(value, Engine.Z))
    _line_product_shape(statistics_product, core_shape, "summaries")
    all(summary -> summary isa SampleSummary && summary.n == trials,
        Iterators.flatten(Base.values(statistics_product))) || throw(DimensionMismatch(
        "line-parameter summaries must contain the Monte Carlo trial count",
    ))
    if sample_product !== nothing
        sample_shape = (core_shape..., trials)
        _line_product_shape(sample_product, sample_shape, "samples")
    end
    if histogram_product !== nothing
        _line_product_shape(histogram_product, core_shape, "histograms")
        all(histogram -> histogram isa HistogramDensity,
            Iterators.flatten(Base.values(histogram_product))) || throw(ArgumentError(
            "line-parameter histograms must contain HistogramDensity values",
        ))
    end
    return nothing
end

function _validate_monte_carlo_products(
        values::Vector{T},
        statistics_products,
        sample_products,
        histogram_products,
        trial_counts
) where {T <: Union{DataModel.CableConstants, Engine.LineParameters}}
    expected_keys = _monte_carlo_keys(T)
    for point in eachindex(values)
        trials = trial_counts[point]
        trials > 0 || throw(ArgumentError("Monte Carlo trial counts must be positive"))
        statistics_product = _validated_product(
            statistics_products[point], expected_keys, "statistics products")
        sample_product = sample_products === nothing ? nothing : _validated_product(
            sample_products[point], expected_keys, "sample products")
        histogram_product = histogram_products === nothing ? nothing : _validated_product(
            histogram_products[point], expected_keys, "histogram products")
        if values[point] isa DataModel.CableConstants
            _validate_cable_products(
                statistics_product, sample_product, histogram_product, trials)
        else
            _validate_line_products(
                values[point], statistics_product, sample_product, histogram_product, trials)
        end
    end
    return nothing
end

function _validate_monte_carlo_products(
        ::Vector,
        statistics_products,
        sample_products,
        histogram_products,
        trial_counts
)
    all(>(0), trial_counts) || throw(ArgumentError(
        "Monte Carlo trial counts must be positive",
    ))
    return nothing
end

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
        isempty(values) && throw(ArgumentError(
            "MonteCarloResult requires at least one core result",
        ))
        isconcretetype(eltype(stats)) || throw(ArgumentError(
            "Monte Carlo statistics must use a concrete product type",
        ))
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
        sample_values === nothing || sample_values isa AbstractVector || throw(
            ArgumentError("retained samples must be stored in a vector"),
        )
        histogram_values === nothing || histogram_values isa AbstractVector || throw(
            ArgumentError("retained histograms must be stored in a vector"),
        )
        sample_values === nothing || isconcretetype(eltype(sample_values)) || throw(
            ArgumentError("retained samples must use a concrete product type"),
        )
        histogram_values === nothing || isconcretetype(eltype(histogram_values)) || throw(
            ArgumentError("retained histograms must use a concrete product type"),
        )
        _validate_monte_carlo_products(
            values,
            stats,
            sample_values,
            histogram_values,
            trial_counts
        )
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
