"""
$(TYPEDEF)

Store scalar Monte Carlo statistics for one real-valued observation.

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
    isempty(values) && throw(ArgumentError("cannot summarise an empty sample"))
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

Statistics.mean(summary::SampleSummary) = summary.mean
Statistics.std(summary::SampleSummary) = summary.std
Statistics.median(summary::SampleSummary) = summary.median
Base.minimum(summary::SampleSummary) = summary.min
Base.maximum(summary::SampleSummary) = summary.max

"""
$(TYPEDEF)

Store a normalised piecewise-constant probability density. `edges` use the
units of the sampled quantity and `density` uses their reciprocal.

$(TYPEDFIELDS)
"""
struct HistogramDensity{T <: AbstractFloat}
    "Strictly increasing bin edges."
    edges::Vector{T}
    "Normalised density in each bin."
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

Group resistance, inductance, capacitance, and conductance products.

`Basis` is `:per_length` or `:total` and applies to all four fields. `T` is
the common representation used for samples, summaries, or histograms.

$(TYPEDFIELDS)
"""
struct RLCG{T, Basis}
    "Resistance values."
    R::T
    "Inductance values."
    L::T
    "Capacitance values."
    C::T
    "Conductance values."
    G::T
end


function RLCG(R::T, L::T, C::T, G::T; basis::Symbol = :per_length) where {T}
    basis in (:per_length, :total) || throw(
        ArgumentError("basis must be :per_length or :total; got :$basis"),
    )
    return RLCG{T, basis}(R, L, C, G)
end

basis(::RLCG{T, Basis}) where {T, Basis} = Basis
