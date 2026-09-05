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

_monte_carlo_keys(::Type{<:Engine.CableConstants}) = (:R, :L, :C, :G)
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
        value::Engine.CableConstants,
        statistics_product,
        sample_product,
        histogram_product,
        trials::Int
)
    count = length(value)
    all(field -> field isa AbstractVector && length(field) == count,
        Base.values(statistics_product)) || throw(DimensionMismatch(
        "cable-constant summaries must match the assembly count",
    ))
    all(summary -> summary isa SampleSummary && summary.n == trials,
        Iterators.flatten(Base.values(statistics_product))) || throw(DimensionMismatch(
        "cable-constant summaries must contain the Monte Carlo trial count",
    ))
    if sample_product !== nothing
        all(sample -> sample isa AbstractMatrix && size(sample) == (count, trials),
            Base.values(sample_product)) || throw(DimensionMismatch(
            "cable-constant samples must contain assembly and trial dimensions",
        ))
    end
    if histogram_product !== nothing
        all(field -> field isa AbstractVector && length(field) == count,
            Base.values(histogram_product)) || throw(DimensionMismatch(
            "cable-constant histograms must match the assembly count",
        ))
        all(histogram -> histogram isa HistogramDensity,
            Iterators.flatten(Base.values(histogram_product))) || throw(ArgumentError(
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
) where {T <: Union{Engine.CableConstants, Engine.LineParameters}}
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
        if values[point] isa Engine.CableConstants
            _validate_cable_products(
                values[point], statistics_product, sample_product,
                histogram_product, trials)
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

function _validate_failure_record(record)
    record isa NamedTuple &&
        keys(record) == (:attempt, :target_trial, :stage, :sample, :error) ||
        throw(ArgumentError(
        "Monte Carlo failure records must contain attempt, target_trial, stage, sample, and error",
    ))
    record.attempt isa Int && record.attempt > 0 || throw(ArgumentError(
        "Monte Carlo failure attempts must be positive integers",
    ))
    record.target_trial isa Int && record.target_trial > 0 || throw(ArgumentError(
        "Monte Carlo failure target trials must be positive integers",
    ))
    record.stage in (:sample, :build, :compute) || throw(ArgumentError(
        "Monte Carlo failure stages must be :sample, :build, or :compute",
    ))
    error = record.error
    error isa NamedTuple && keys(error) == (:type, :message, :stack) || throw(
        ArgumentError(
            "Monte Carlo failure errors must contain type, message, and stack",
        ),
    )
    error.type isa String && error.message isa String || throw(ArgumentError(
        "Monte Carlo failure type and message summaries must be strings",
    ))
    error.stack isa AbstractVector || throw(ArgumentError(
        "Monte Carlo failure stacks must be vectors",
    ))
    all(error.stack) do frame
        frame isa NamedTuple &&
            keys(frame) == (:function_name, :file, :line) &&
            frame.function_name isa String && frame.file isa String && frame.line isa Int
    end || throw(ArgumentError(
        "Monte Carlo failure stack frames must contain function_name, file, and line",
    ))
    return nothing
end

function _validate_failure_summary(summary, failures, accepted::Int)
    summary isa NamedTuple && keys(summary) == (
        :attempts, :accepted, :failed, :acceptance_rate, :by_type, :by_stage
    ) || throw(ArgumentError(
        "Monte Carlo failure summaries have an invalid schema",
    ))
    summary.attempts isa Int && summary.attempts > 0 || throw(ArgumentError(
        "Monte Carlo failure-summary attempts must be positive integers",
    ))
    summary.accepted == accepted || throw(DimensionMismatch(
        "Monte Carlo failure-summary accepted counts must match trial counts",
    ))
    summary.failed == length(failures) || throw(DimensionMismatch(
        "Monte Carlo failure-summary failed counts must match retained failures",
    ))
    summary.attempts == summary.accepted + summary.failed || throw(DimensionMismatch(
        "Monte Carlo failure-summary attempts must equal accepted plus failed trials",
    ))
    summary.acceptance_rate == summary.accepted / summary.attempts || throw(
        ArgumentError(
            "Monte Carlo failure-summary acceptance rates are inconsistent",
        ),
    )
    return nothing
end

function _validate_monte_carlo_details(details, values, trial_counts)
    isempty(details) && return nothing
    keys(details) == (:trials, :failures, :failure_summary) || throw(ArgumentError(
        "MonteCarloResult details must contain trials, failures, and failure_summary",
    ))
    point_count = length(values)
    all(length(product) == point_count
    for product in (details.trials, details.failures, details.failure_summary)) ||
        throw(DimensionMismatch(
        "retained Monte Carlo details must contain one entry per Gridspace point",
    ))
    for point in eachindex(values)
        records = details.trials[point]
        failures = details.failures[point]
        summary = details.failure_summary[point]
        length(records) == trial_counts[point] || throw(DimensionMismatch(
            "retained details must contain one entry per accepted Monte Carlo trial",
        ))
        failures isa AbstractVector || throw(ArgumentError(
            "retained Monte Carlo failures must be stored in vectors",
        ))
        foreach(_validate_failure_record, failures)
        _validate_failure_summary(summary, failures, trial_counts[point])
    end
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
    "Accepted-trial count used for each Gridspace point."
    trial_counts::Vector{Int}
    "Typed accepted-trial details and failure diagnostics retained by the propagation."
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
        _validate_monte_carlo_details(details, values, trial_counts)
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

"""
$(TYPEDSIGNATURES)

Transport every uncertainty-bearing linear-error result into one
target-bearing downstream problem point while preserving source cardinality
and order.
"""
function ParametricBuilder.Gridspace{Target}(
        source::LinearErrorResult
) where {Target}
    return ParametricBuilder.Gridspace{Target}(Target, (source,))
end

function ParametricBuilder.Gridspace{Target}(
        source::MonteCarloResult{T}
) where {Target, T}
    target_name = nameof(Target)
    throw(ArgumentError(
        "Gridspace transport from MonteCarloResult to problem $target_name requires " *
        "a reconstruction for result type $T. For built-in cable and line results, " *
        "load Measurements.jl with `using Measurements`. See `?Gridspace` and " *
        "https://electa-git.github.io/LineCableModels.jl/dev/gridspace/" *
        "#Transporting-completed-result-spaces",
    ))
end
