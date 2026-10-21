"""
$(TYPEDEF)

Store ordered core results from a [`LinearError`](@ref) computation.

$(TYPEDFIELDS)
"""
struct LinearErrorResult{T, F} <: AbstractUncertaintyResult{T}
    "Higher-order formulation used for the computation."
    formulation::F
    "Uncertainty-bearing core results in Gridspace traversal order."
    values::Vector{T}
    "Supplemental output retained by the propagation."
    details::ComputationDetails

    function LinearErrorResult(
            formulation::F,
            values::Vector{T},
            details::ComputationDetails
    ) where {T, F}
        validate(T, LinearErrorResult)
        isempty(details.data) || keys(details.data) == (:points,) ||
            throw(ArgumentError(
                "LinearErrorResult details must be empty or contain only points",
            ))
        isempty(details.data) || length(details.data.points) == length(values) ||
            throw(DimensionMismatch(
                "retained details must contain one entry per core result",
            ))
        return new{T, F}(formulation, values, details)
    end
end

LinearErrorResult(formulation, values) = LinearErrorResult(formulation, values, ComputationDetails())

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

function _line_product_shape(product, expected_shape, name::AbstractString)
    all(field -> field isa AbstractArray && size(field) == expected_shape,
        Base.values(product)) || throw(DimensionMismatch(
        "line-parameter $name must match the core matrix and frequency dimensions",
    ))
    return nothing
end


"""
$(TYPEDEF)

Store uncertainty-bearing core results and their empirical Monte Carlo products.
The marginal representation is materialized during computation. Access reuses
the retained uncertainty-source identities.

$(TYPEDFIELDS)
"""
struct MonteCarloResult{T, F, ST <: AbstractVector, S, H} <: AbstractUncertaintyResult{T}
    "Higher-order formulation used for the computation."
    formulation::F
    "Stored mean ± sample standard deviation core results in Gridspace traversal order."
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
    "Accepted-trial details retained by the propagation."
    details::ComputationDetails

    function MonteCarloResult(
            formulation::F,
            values::Vector{T},
            stats::ST,
            sample_values::S,
            histogram_values::H,
            root_seed::UInt64,
            point_seeds::Vector{UInt64},
            trial_counts::Vector{Int},
            details::ComputationDetails
    ) where {T, F, ST <: AbstractVector, S, H}
        return validate(new{T, F, ST, S, H}(
            formulation,
            values,
            stats,
            sample_values,
            histogram_values,
            root_seed,
            point_seeds,
            trial_counts,
            details
        ))
    end
end

# Checks the stored cores, the empirical products of each point against its core result,
# and the retained trial and timing details.
function validate(result::MonteCarloResult{T}) where {T}
    (; values, stats, sample_values, histogram_values, point_seeds, trial_counts, details) =
        result
    validate(T, MonteCarloResult)
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
    if T <: Union{Engine.CableConstants, Engine.LineParameters}
        Engine.has_uncertainty_type(eltype(observe(first(values),R))) || throw(ArgumentError(
            "MonteCarloResult requires stored uncertainty-bearing cores; materialize the empirical summaries before construction"))
        expected_keys = _monte_carlo_keys(T)
        for point in eachindex(values)
            trials = trial_counts[point]
            trials > 0 || throw(ArgumentError("Monte Carlo trial counts must be positive"))
            statistics_product = _validated_product(
                stats[point], expected_keys, "statistics products")
            sample_product = sample_values === nothing ? nothing : _validated_product(
                sample_values[point], expected_keys, "sample products")
            histogram_product = histogram_values === nothing ? nothing : _validated_product(
                histogram_values[point], expected_keys, "histogram products")
            value = values[point]
            if value isa Engine.CableConstants
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
            else
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
            end
        end
    else
        all(>(0), trial_counts) || throw(ArgumentError(
            "Monte Carlo trial counts must be positive",
        ))
    end
    if !isempty(details.data)
        fields = Tuple(key for key in keys(details.data) if key !== :timing)
        fields in ((), (:trials,)) || throw(ArgumentError(
            "MonteCarloResult details must contain timing, trials, or both",
        ))
        point_count = length(values)
        if haskey(details.data, :timing)
            details.data.timing isa AbstractVector || throw(ArgumentError(
                "Monte Carlo timing must be a vector of population vectors"))
        end
        all(length(product) == point_count
        for product in Base.values(details.data)) ||
            throw(DimensionMismatch(
            "retained Monte Carlo details must contain one entry per Gridspace point",
        ))
        if haskey(details.data, :timing)
            timing = details.data.timing
            for point in eachindex(values)
                records = timing[point]
                records isa AbstractVector && all(record -> record isa NamedTuple, records) ||
                    throw(ArgumentError("Monte Carlo timing must contain vectors of named tuples"))
                length(records) == trial_counts[point] || throw(DimensionMismatch(
                    "Monte Carlo timing must contain one record per accepted trial"))
            end
        end
        if !isempty(fields)
            for point in eachindex(values)
                records = details.data.trials[point]
                length(records) == trial_counts[point] || throw(DimensionMismatch(
                    "retained details must contain one entry per accepted Monte Carlo trial",
                ))
            end
        end
    end
    return result
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
        ComputationDetails()
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
function Commons.Gridspace{Target}(
        source::LinearErrorResult
) where {Target}
    return Commons.Gridspace{Target}(Target, (source,))
end

function Commons.Gridspace{Target}(
        source::MonteCarloResult{T}
) where {Target, T}
    target_name = nameof(Target)
    throw(ArgumentError(
        "Gridspace transport from MonteCarloResult to problem $target_name requires " *
        "a transport for result type $T. See `?Gridspace` and " *
        "https://electa-git.github.io/LineCableModels.jl/dev/gridspace/" *
        "#Transporting-completed-result-spaces",
    ))
end

function Commons.Gridspace{Target}(
        source::MonteCarloResult{<:Union{Engine.CableConstants,Engine.LineParameters}}
) where {Target}
    return Commons.Gridspace{Target}(Target, (source,))
end

"""
$(TYPEDSIGNATURES)

Materialize retained R/L/C/G marginal summaries as an uncertainty-bearing core
result using `source`'s physical axes. The summaries include native units and
their standard deviations describe output spread, not uncertainty of the mean.
`details` defaults to the source's supplemental output. MC computation selects
only its point-wide metadata. Neither histogram bins nor sample retention
participate in this construction.
"""
function materialize(source::Engine.CableConstants,
        summaries::NamedTuple{(:R,:L,:C,:G)}; details::ComputationDetails=source.details)
    return Engine.CableConstants(source.terminals,
        materialize.(summaries.R), materialize.(summaries.L),
        materialize.(summaries.C), materialize.(summaries.G), source.frequency, details)
end

function materialize(source::Engine.LineParameters,
        summaries::NamedTuple{(:R,:L,:C,:G)}; details::ComputationDetails=source.details)
    angular = reshape(2π .* frequencies(source), 1, 1, :)
    return Engine.LineParameters(source.domain,
        complex.(materialize.(summaries.R), materialize.(summaries.L) .* angular),
        complex.(materialize.(summaries.G), materialize.(summaries.C) .* angular),
        frequencies(source); terminals=source.terminals,
        line_length=source.line_length, details)
end
