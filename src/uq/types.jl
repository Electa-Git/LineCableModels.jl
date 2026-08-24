"""
$(TYPEDEF)

Apply established direct linear uncertainty propagation with `inner`.

$(TYPEDFIELDS)
"""
struct LinearError{F <: AbstractFormulation} <: AbstractFormulation
    "Formulation used for each materialized problem."
    inner::F
    "Invalid-configuration policy: `:error` or `:skip`."
    invalid::Symbol

    function LinearError(inner::F; invalid::Symbol = :error) where {F <:
                                                                    AbstractFormulation}
        invalid in (:error, :skip) || throw(ArgumentError(
            "invalid must be :error or :skip; got $(repr(invalid))",
        ))
        return new{F}(inner, invalid)
    end
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
    "Invalid-configuration policy: `:error` or `:skip`."
    invalid::Symbol

    function MonteCarlo(
            inner::F;
            trials::Union{Nothing, Integer} = nothing,
            confidence::Real = 0.95,
            cdf_tol::Real = 0.02,
            distribution = :normal,
            seed::Union{Nothing, Integer} = nothing,
            return_samples::Bool = false,
            return_histograms::Bool = false,
            bins::Union{Nothing, Integer} = nothing,
            invalid::Symbol = :error
    ) where {F <: AbstractFormulation}
        trials === nothing || trials > 0 || throw(ArgumentError("trials must be positive"))
        0 < confidence < 1 ||
            throw(ArgumentError("confidence must lie between zero and one"))
        0 < cdf_tol < 1 || throw(ArgumentError("cdf_tol must lie between zero and one"))
        bins === nothing || bins > 0 || throw(ArgumentError("bins must be positive"))
        invalid in (:error, :skip) || throw(ArgumentError(
            "invalid must be :error or :skip; got $(repr(invalid))",
        ))
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
            bins === nothing ? nothing : Int(bins),
            invalid
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
struct LinearErrorResult{T, F, S, D <: Dict{Symbol, NamedTuple}} <:
       AbstractUncertaintyResult{T}
    "Resolved higher-order formulation."
    formulation::F
    "Successful uncertainty-bearing primitive results."
    values::Vector{T}
    "Source parameter space."
    space::S
    "Failures, retained data, replay data, and manifest entries."
    details::D

    function LinearErrorResult(formulation::F, values::Vector{T}, space::S,
            details::D) where {T, F, S, D <: Dict{Symbol, NamedTuple}}
        ParametricBuilder._primitive_result_type(T)
        ParametricBuilder._validate_details(details)
        return new{T, F, S, D}(formulation, values, space, details)
    end
end

"""
$(TYPEDEF)

Store primitive mean representations and real-valued Monte Carlo summaries.

$(TYPEDFIELDS)
"""
struct MonteCarloResult{T, F, S, ST <: AbstractVector, D <: Dict{Symbol, NamedTuple}} <:
       AbstractUncertaintyResult{T}
    "Resolved higher-order formulation."
    formulation::F
    "Successful primitive mean representations."
    values::Vector{T}
    "Source parameter space."
    space::S
    "Per-observable sample summaries."
    stats::ST
    "Failures, retained data, replay data, and manifest entries."
    details::D

    function MonteCarloResult(formulation::F, values::Vector{T}, space::S,
            stats::ST, details::D) where {
            T, F, S, ST <: AbstractVector, D <: Dict{Symbol, NamedTuple}}
        ParametricBuilder._primitive_result_type(T)
        ParametricBuilder._validate_details(details)
        length(stats) == length(values) || throw(DimensionMismatch(
            "Monte Carlo statistics must contain one entry per primitive result",
        ))
        return new{T, F, S, ST, D}(formulation, values, space, stats, details)
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
samples(value::MonteCarloResult) = value.details[:samples].values
histograms(value::MonteCarloResult) = value.details[:histograms].values
uncertain_value(value::LinearErrorResult) = value.values
manifest(value::Union{LinearErrorResult, MonteCarloResult}) = value.details[:manifest].value

function observables(value::LinearErrorResult)
    return (
        result = result(value),
        details = value.details,
        manifest = manifest(value)
    )
end

function observables(value::MonteCarloResult)
    return (
        result = result(value),
        statistics = statistics(value),
        samples = samples(value),
        histograms = histograms(value),
        details = value.details,
        manifest = manifest(value)
    )
end
