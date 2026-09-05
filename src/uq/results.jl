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

const _POLYNOMIAL_CHAOS_KEYS = (:R, :L, :C, :G)
const _POLYNOMIAL_CHAOS_VALIDATION_KEYS = (
    :stochastic_dimension,
    :polynomial_degree,
    :basis_terms,
    :collocation_evaluations,
    :validation_evaluations,
    :relative_rms_error,
    :relative_max_error,
    :worst,
)

_polynomial_chaos_shape(value::Engine.CableConstants) = size(value)
_polynomial_chaos_shape(value::Engine.LineParameters) = size(observe(value, R))

function _polynomial_chaos_product(product, name::AbstractString)
    product isa NamedTuple || throw(ArgumentError(
        "polynomial-chaos $name must be a named tuple",
    ))
    keys(product) == _POLYNOMIAL_CHAOS_KEYS || throw(ArgumentError(
        "polynomial-chaos $name must contain exactly R, L, C, and G",
    ))
    return product
end

function _validate_polynomial_chaos_expansion(expansion, shape::Tuple, terms::Int)
    expansion isa NamedTuple && keys(expansion) == (:basis, :coefficients) || throw(
        ArgumentError("polynomial-chaos expansions must contain basis and coefficients"),
    )
    coefficients = _polynomial_chaos_product(expansion.coefficients, "coefficients")
    expected_shape = (shape..., terms)
    for coefficient in Base.values(coefficients)
        coefficient isa AbstractArray && size(coefficient) == expected_shape || throw(
            DimensionMismatch(
                "polynomial-chaos coefficient arrays must have native dimensions followed by the basis-term dimension",
            ),
        )
        all(value -> value isa Real && isfinite(value), coefficient) || throw(
            ArgumentError("polynomial-chaos coefficients must be finite real values"),
        )
    end
    return nothing
end

function _validate_polynomial_chaos_statistics(stats, shape::Tuple)
    product = _polynomial_chaos_product(stats, "statistics")
    for summary in Base.values(product)
        summary isa NamedTuple && keys(summary) == (:mean, :std) || throw(
            ArgumentError("polynomial-chaos statistics must contain exactly mean and std"),
        )
        for field in Base.values(summary)
            field isa AbstractArray && size(field) == shape || throw(
                DimensionMismatch("polynomial-chaos statistics must match native result dimensions"),
            )
            all(value -> value isa Real && isfinite(value), field) || throw(
                ArgumentError("polynomial-chaos statistics must be finite real values"),
            )
        end
        all(>=(0), summary.std) || throw(ArgumentError(
            "polynomial-chaos standard deviations must be nonnegative",
        ))
    end
    return nothing
end

function _validate_polynomial_chaos_error_product(product, name::AbstractString)
    errors = _polynomial_chaos_product(product, name)
    all(Base.values(errors)) do value
        value isa Real && isfinite(value) && value >= 0
    end || throw(ArgumentError(
        "polynomial-chaos $name must contain finite nonnegative values",
    ))
    return nothing
end

function _validate_polynomial_chaos_record(record, formulation::PolynomialChaos)
    record isa NamedTuple && keys(record) == _POLYNOMIAL_CHAOS_VALIDATION_KEYS || throw(
        ArgumentError("polynomial-chaos validation records have an invalid schema"),
    )
    record.stochastic_dimension isa Int && record.stochastic_dimension >= 0 || throw(
        ArgumentError("stochastic dimensions must be nonnegative integers"),
    )
    record.polynomial_degree == formulation.degree || throw(DimensionMismatch(
        "validation polynomial degrees must match the formulation",
    ))
    record.basis_terms isa Int && record.basis_terms >= 1 || throw(ArgumentError(
        "basis-term counts must be positive integers",
    ))
    record.collocation_evaluations isa Int && record.collocation_evaluations >= 1 || throw(
        ArgumentError("collocation evaluation counts must be positive integers"),
    )
    record.validation_evaluations == formulation.validation_points || throw(
        DimensionMismatch("validation evaluation counts must match the formulation"),
    )
    _validate_polynomial_chaos_error_product(
        record.relative_rms_error, "relative RMS errors")
    _validate_polynomial_chaos_error_product(
        record.relative_max_error, "relative maximum errors")
    worst = record.worst
    worst isa NamedTuple && keys(worst) == (
        :quantity, :index, :frequency, :relative_rms_error, :relative_max_error
    ) || throw(ArgumentError(
        "worst validation locations must contain quantity, index, frequency, relative_rms_error, and relative_max_error",
    ))
    worst.quantity in _POLYNOMIAL_CHAOS_KEYS || throw(ArgumentError(
        "worst validation quantities must be R, L, C, or G",
    ))
    worst.index isa Tuple && all(index -> index isa Int, worst.index) || throw(
        ArgumentError("worst validation indices must be integer tuples"),
    )
    worst.frequency === nothing ||
        worst.frequency isa Real && isfinite(worst.frequency) || throw(
        ArgumentError("worst validation frequencies must be finite or nothing"),
    )
    all(field -> field isa Real && isfinite(field) && field >= 0,
        (worst.relative_rms_error, worst.relative_max_error)) || throw(
        ArgumentError("worst validation errors must be finite and nonnegative"),
    )
    return nothing
end

function _validate_polynomial_chaos_details(
        details::ComputationDetails,
        records::AbstractVector
)
    isempty(details) && return nothing
    keys(details) == (:points,) || throw(ArgumentError(
        "PolynomialChaosResult details must be empty or contain only points",
    ))
    length(details.points) == length(records) || throw(DimensionMismatch(
        "retained polynomial-chaos details must contain one entry per core result",
    ))
    for (point, record) in zip(details.points, records)
        point isa NamedTuple && keys(point) == (:collocation, :validation) || throw(
            ArgumentError(
                "retained polynomial-chaos point details must contain collocation and validation",
            ),
        )
        point.collocation isa AbstractVector && point.validation isa AbstractVector || throw(
            ArgumentError("retained polynomial-chaos solves must be stored in vectors"),
        )
        length(point.collocation) == record.collocation_evaluations || throw(
            DimensionMismatch("retained collocation solves must match validation records"),
        )
        length(point.validation) == record.validation_evaluations || throw(
            DimensionMismatch("retained validation solves must match validation records"),
        )
    end
    return nothing
end

"""
$(TYPEDEF)

Store core means, spectral expansions, coefficient-derived moments, and
independent validation diagnostics from a [`PolynomialChaos`](@ref)
calculation.

$(TYPEDFIELDS)
"""
struct PolynomialChaosResult{
    T,
    F,
    E <: AbstractVector,
    S <: AbstractVector,
    V <: AbstractVector,
    D <: ComputationDetails,
} <: AbstractUncertaintyResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Core results reconstructed from expansion means in Gridspace order."
    values::Vector{T}
    "Per-point basis and native coefficient products."
    expansion_values::E
    "Per-point native mean and standard-deviation products."
    stats::S
    "Per-point independent surrogate-validation diagnostics."
    validation_values::V
    "Typed collocation and validation solve details retained by the propagation."
    details::D

    function PolynomialChaosResult(
            formulation::F,
            values::Vector{T},
            expansion_values::E,
            stats::S,
            validation_values::V,
            details::D
    ) where {
            T,
            F,
            E <: AbstractVector,
            S <: AbstractVector,
            V <: AbstractVector,
            D <: ComputationDetails
    }
        formulation isa PolynomialChaos || throw(ArgumentError(
            "PolynomialChaosResult requires a PolynomialChaos formulation",
        ))
        check_core_result(T)
        T <: Union{Engine.CableConstants, Engine.LineParameters} || throw(
            ArgumentError(
                "PolynomialChaosResult supports CableConstants and LineParameters core results",
            ),
        )
        isempty(values) && throw(ArgumentError(
            "PolynomialChaosResult requires at least one core result",
        ))
        for (name, product) in (
                ("expansions", expansion_values),
                ("statistics", stats),
                ("validation records", validation_values))
            length(product) == length(values) || throw(DimensionMismatch(
                "polynomial-chaos $name must contain one entry per core result",
            ))
            isconcretetype(eltype(product)) || throw(ArgumentError(
                "polynomial-chaos $name must use a concrete product type",
            ))
        end
        for point in eachindex(values)
            record = validation_values[point]
            _validate_polynomial_chaos_record(record, formulation)
            shape = _polynomial_chaos_shape(values[point])
            _validate_polynomial_chaos_expansion(
                expansion_values[point], shape, record.basis_terms)
            _validate_polynomial_chaos_statistics(stats[point], shape)
        end
        _validate_polynomial_chaos_details(details, validation_values)
        return new{T, F, E, S, V, D}(
            formulation,
            values,
            expansion_values,
            stats,
            validation_values,
            details
        )
    end
end

function PolynomialChaosResult(
        formulation,
        values,
        expansion_values,
        stats,
        validation_values
)
    return PolynomialChaosResult(
        formulation, values, expansion_values, stats, validation_values, (;))
end

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
details(value::PolynomialChaosResult) = value.details

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
