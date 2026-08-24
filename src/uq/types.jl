"""
$(TYPEDEF)

Apply established direct linear uncertainty propagation with `inner`.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation} <: AbstractFormulation
    "Formulation used for each materialized problem."
    inner::F
end

"""
$(TYPEDEF)

Select conditional Monte Carlo propagation over a
[`ParametricProblem`](@ref). Randomness is local and reproducible when `seed`
is supplied.

$(TYPEDFIELDS)
"""
struct MonteCarlo{F <: AbstractFormulation, D, S} <: AbstractFormulation
    "Formulation used for each sampled problem."
    inner::F
    "Requested trials, or `nothing` for DKW sizing."
    trials::Union{Nothing, Int}
    "Simultaneous empirical-CDF confidence [dimensionless]."
    confidence::Float64
    "Maximum empirical-CDF deviation used for DKW sizing [dimensionless]."
    cdf_tol::Float64
    "Sampling family or extension-provided univariate distribution."
    distribution::D
    "Optional root random seed."
    seed::S
    "Whether joint primitive samples are retained."
    return_samples::Bool
    "Whether marginal histogram densities are retained."
    return_histograms::Bool
    "Optional histogram bin count."
    bins::Union{Nothing, Int}
    function MonteCarlo(
            inner::F;
            trials::Union{Nothing, Integer} = nothing,
            confidence::Real = 0.95,
            cdf_tol::Real = 0.02,
            distribution = :normal,
            seed::Union{Nothing, Integer} = nothing,
            return_samples::Bool = false,
            return_histograms::Bool = false,
            bins::Union{Nothing, Integer} = nothing
    ) where {F <: AbstractFormulation}
        trials === nothing || trials > 0 || throw(ArgumentError("trials must be positive"))
        0 < confidence < 1 ||
            throw(ArgumentError("confidence must lie between zero and one"))
        0 < cdf_tol < 1 || throw(ArgumentError("cdf_tol must lie between zero and one"))
        bins === nothing || bins > 0 || throw(ArgumentError("bins must be positive"))
        distribution isa Symbol && distribution ∉ (:normal, :uniform) &&
            throw(ArgumentError(
                "unsupported distribution $(repr(distribution)); expected :normal, :uniform, a sampler function, or an extension-supported distribution",
            ))
        actual_seed = seed === nothing ? nothing : UInt64(seed)
        return new{F, typeof(distribution), typeof(actual_seed)}(
            inner,
            trials === nothing ? nothing : Int(trials),
            Float64(confidence),
            Float64(cdf_tol),
            distribution,
            actual_seed,
            return_samples,
            return_histograms,
            bins === nothing ? nothing : Int(bins)
        )
    end
end

"""
$(TYPEDEF)

Store scalar Monte Carlo statistics for one declared real observable.

$(TYPEDFIELDS)
"""
struct SampleSummary{T <: Real}
    "Arithmetic mean."
    mean::T
    "Sample standard deviation."
    std::T
    "Minimum."
    min::T
    "5th percentile."
    q05::T
    "Median."
    median::T
    "95th percentile."
    q95::T
    "Maximum."
    max::T
    "Number of successful samples."
    n::Int

    function SampleSummary(
            mean::T,
            std::T,
            min::T,
            q05::T,
            median::T,
            q95::T,
            max::T,
            n::Integer
    ) where {T <: Real}
        all(isfinite, (mean, std, min, q05, median, q95, max)) || throw(
            ArgumentError("sample statistics must be finite"),
        )
        std >= zero(T) || throw(ArgumentError(
            "sample standard deviation must be nonnegative",
        ))
        min <= q05 <= median <= q95 <= max || throw(ArgumentError(
            "sample quantiles must be ordered between the minimum and maximum",
        ))
        n > 0 || throw(ArgumentError("sample count must be positive"))
        return new{T}(mean, std, min, q05, median, q95, max, Int(n))
    end
end

function SampleSummary(mean::Real, std::Real, min::Real, q05::Real,
        median::Real, q95::Real, max::Real, n::Integer)
    promoted = promote(mean, std, min, q05, median, q95, max)
    return SampleSummary(promoted..., n)
end

function SampleSummary(values::AbstractVector{<:Real})
    isempty(values) && throw(ArgumentError("cannot summarize an empty sample"))
    all(isfinite, values) || throw(ArgumentError("sample values must be finite"))
    sigma = length(values) == 1 ? zero(float(first(values))) : Statistics.std(values)
    promoted = promote(
        Statistics.mean(values),
        sigma,
        minimum(values),
        Statistics.quantile(values, 0.05),
        Statistics.quantile(values, 0.50),
        Statistics.quantile(values, 0.95),
        maximum(values)
    )
    return SampleSummary(promoted..., length(values))
end

"""
$(TYPEDEF)

Store a normalized piecewise-constant probability density. `edges` use the
units of the sampled quantity and `density` uses their reciprocal.

$(TYPEDFIELDS)
"""
struct HistogramDensity{T <: AbstractFloat}
    "Strictly increasing bin edges."
    edges::Vector{T}
    "Normalized density in each bin."
    density::Vector{T}

    function HistogramDensity(edges::Vector{T}, density::Vector{T}) where {T <:
                                                                           AbstractFloat}
        length(edges) == length(density) + 1 || throw(ArgumentError(
            "histogram edges must contain one more value than density",
        ))
        isempty(density) && throw(ArgumentError("histogram density cannot be empty"))
        all(isfinite, edges) || throw(ArgumentError("histogram edges must be finite"))
        widths = diff(edges)
        all(>(zero(T)), widths) || throw(ArgumentError(
            "histogram edges must be strictly increasing",
        ))
        all(value -> isfinite(value) && value >= zero(T), density) ||
            throw(ArgumentError("histogram density must be finite and nonnegative"))
        area = sum(density .* widths)
        area > zero(T) || throw(ArgumentError("histogram density must have positive area"))
        return new{T}(copy(edges), density ./ area)
    end
end

function HistogramDensity(edges::AbstractVector{<:Real}, density::AbstractVector{<:Real})
    T = promote_type(float(eltype(edges)), float(eltype(density)))
    isconcretetype(T) || (T = Float64)
    return HistogramDensity(Vector{T}(edges), Vector{T}(density))
end

"""
$(TYPEDEF)

Group resistance, inductance, capacitance, and conductance representations.

$(TYPEDFIELDS)
"""
struct RLCG{T}
    "Resistance values."
    R::T
    "Inductance values."
    L::T
    "Capacitance values."
    C::T
    "Conductance values."
    G::T
end

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

Store primitive mean representations and real-valued Monte Carlo summaries.

$(TYPEDFIELDS)
"""
struct MonteCarloResult{T, F, ST <: AbstractVector, S, H} <:
       AbstractUncertaintyResult{T}
    "Higher-order formulation used for the calculation."
    formulation::F
    "Primitive mean representations in Gridspace traversal order."
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

for ResultType in (:LinearErrorResult, :MonteCarloResult)
    @eval begin
        Base.size(value::$ResultType) = size(value.values)
        Base.length(value::$ResultType) = length(value.values)
        Base.getindex(value::$ResultType, index::Integer) = value.values[index]
        Base.iterate(value::$ResultType, state...) = iterate(value.values, state...)
        Base.IndexStyle(::Type{<:$ResultType}) = IndexLinear()
    end
end

result(value::Union{LinearErrorResult, MonteCarloResult}) = value.values
statistics(value::MonteCarloResult) = value.stats
samples(value::MonteCarloResult) = value.sample_values
histograms(value::MonteCarloResult) = value.histogram_values
uncertain_value(value::LinearErrorResult) = value.values
observables(value::LinearErrorResult) = (result = result(value),)

function observables(
        value::MonteCarloResult{<:Any, <:Any, <:AbstractVector, Nothing, Nothing},
)
    return (result = result(value), statistics = statistics(value))
end
function observables(
        value::MonteCarloResult{<:Any, <:Any, <:AbstractVector, <:Any, Nothing},
)
    return (
        result = result(value),
        statistics = statistics(value),
        samples = samples(value)
    )
end
function observables(
        value::MonteCarloResult{<:Any, <:Any, <:AbstractVector, Nothing, <:Any},
)
    return (
        result = result(value),
        statistics = statistics(value),
        histograms = histograms(value)
    )
end
function observables(value::MonteCarloResult)
    return (
        result = result(value),
        statistics = statistics(value),
        samples = samples(value),
        histograms = histograms(value)
    )
end
