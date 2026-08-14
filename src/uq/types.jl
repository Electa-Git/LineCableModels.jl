"""
    SampleSummary{T}

Store the fixed summary statistics for one scalar Monte Carlo observable.
All fields use the observable's canonical units.
"""
struct SampleSummary{T <: Real}
    "Sample mean."
    mean::T
    "Corrected sample standard deviation."
    std::T
    "Minimum sampled value."
    min::T
    "Five-percent quantile."
    q05::T
    "Median."
    q50::T
    "Ninety-five-percent quantile."
    q95::T
    "Maximum sampled value."
    max::T
end

function Base.:(==)(left::SampleSummary, right::SampleSummary)
    left.mean == right.mean && left.std == right.std && left.min == right.min &&
        left.q05 == right.q05 && left.q50 == right.q50 && left.q95 == right.q95 &&
        left.max == right.max
end

function SampleSummary(values::AbstractVector{T}) where {T <: Real}
    isempty(values) && throw(ArgumentError("cannot summarize an empty sample"))
    sigma = length(values) == 1 ? zero(float(first(values))) : Statistics.std(values)
    promoted = promote(
        Statistics.mean(values),
        sigma,
        minimum(values),
        StatsBase.quantile(values, 0.05),
        StatsBase.quantile(values, 0.50),
        StatsBase.quantile(values, 0.95),
        maximum(values)
    )
    return SampleSummary(promoted...)
end

"""
    RLCG{T}

Group resistance, inductance, capacitance, and conductance values with one
common storage representation. Units follow the basis of the containing
result: \\[Ω/m, H/m, F/m, S/m\\] for `:per_length`, or
\\[Ω, H, F, S\\] for `:total`.
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

function Base.:(==)(left::RLCG, right::RLCG)
    left.R == right.R && left.L == right.L && left.C == right.C && left.G == right.G
end

"""
    HistogramPDF(edges, density)

Piecewise-constant univariate probability density. `edges` has one more entry
than `density`; density is normalized by the constructor.
"""
struct HistogramPDF{T <: Real} <: ContinuousUnivariateDistribution
    "Strictly increasing histogram bin edges."
    edges::Vector{T}
    "Piecewise-constant probability density for each bin."
    density::Vector{T}

    function HistogramPDF{T}(
            edges::Vector{T},
            density::Vector{T}
    ) where {T <: Real}
        return new{T}(edges, density)
    end
end

function Base.:(==)(left::HistogramPDF, right::HistogramPDF)
    left.edges == right.edges && left.density == right.density
end

"""
    CableConstantsMC

Store Monte Carlo summaries and optional retained data for scalar cable
constants. `surrogate` preserves the covariance of \\[R, L, C\\] with
`Measurement` values. Canonical units are \\[Ω/m, H/m, F/m\\].
"""
struct CableConstantsMC{S, Samples, Distributions, Surrogate, T <: Real}
    "Per-observable summary statistics."
    statistics::CableConstants{S}
    "Joint retained samples, or `nothing`."
    samples::Samples
    "Retained marginal histogram distributions, or `nothing`."
    distributions::Distributions
    "Covariance-preserving cable constants."
    surrogate::Surrogate
    "Number of Monte Carlo trials."
    ntrials::Int
    "Central confidence mass used by the analysis."
    confidence::T

    function CableConstantsMC{S, Samples, Distributions, Surrogate, T}(
            statistics::CableConstants{S},
            samples::Samples,
            distributions::Distributions,
            surrogate::Surrogate,
            ntrials::Int,
            confidence::T
    ) where {S, Samples, Distributions, Surrogate, T <: Real}
        return new{S, Samples, Distributions, Surrogate, T}(
            statistics,
            samples,
            distributions,
            surrogate,
            ntrials,
            confidence
        )
    end
end

function CableConstantsMC(
        statistics::CableConstants{S},
        samples::Samples,
        distributions::Distributions,
        surrogate::Surrogate,
        ntrials::Integer,
        confidence::T
) where {S, Samples, Distributions, Surrogate, T <: Real}
    ntrials > 0 || throw(ArgumentError("ntrials must be positive"))
    0 < confidence < 1 || throw(ArgumentError("confidence must lie between zero and one"))
    all(value -> value isa SampleSummary, (statistics.R, statistics.L, statistics.C)) ||
        throw(ArgumentError("cable-constant statistics must contain SampleSummary values"))
    if samples !== nothing
        samples isa CableConstants || throw(
            ArgumentError("cable-constant samples must be CableConstants or nothing"),
        )
        all(value -> value isa AbstractVector{<:Real}, (samples.R, samples.L, samples.C)) ||
            throw(ArgumentError("retained cable-constant samples must be real vectors"))
        all(length(values) == ntrials for values in (samples.R, samples.L, samples.C)) ||
            throw(DimensionMismatch("retained cable-constant samples must match ntrials"))
    end
    if distributions !== nothing
        distributions isa CableConstants || throw(
            ArgumentError("cable-constant distributions must be CableConstants or nothing"),
        )
        all(value -> value isa HistogramPDF, (
            distributions.R,
            distributions.L,
            distributions.C
        )) ||
            throw(ArgumentError("cable-constant distributions must contain HistogramPDF values"))
    end
    surrogate isa CableConstants || throw(
        ArgumentError("the cable-constant surrogate must be a CableConstants value"),
    )
    return CableConstantsMC{S, Samples, Distributions, Surrogate, T}(
        statistics,
        samples,
        distributions,
        surrogate,
        Int(ntrials),
        confidence
    )
end

"""
    LineParametersMC

Store Monte Carlo summaries and optional retained data for frequency-dependent
line parameters. `surrogate` is a covariance-preserving [`LineParameters`](@ref)
whose canonical units and basis match this result.
"""
struct LineParametersMC{S, Samples, Distributions, Surrogate, T <: Real}
    "RLCG summary tensors with dimensions conductor × conductor × frequency."
    statistics::RLCG{S}
    "Joint RLCG sample tensors, or `nothing`."
    samples::Samples
    "Marginal RLCG histogram distributions, or `nothing`."
    distributions::Distributions
    "Covariance-preserving line parameters."
    surrogate::Surrogate
    "Number of Monte Carlo trials."
    ntrials::Int
    "Central confidence mass used by the analysis."
    confidence::T

    function LineParametersMC{S, Samples, Distributions, Surrogate, T}(
            statistics::RLCG{S},
            samples::Samples,
            distributions::Distributions,
            surrogate::Surrogate,
            ntrials::Int,
            confidence::T
    ) where {S, Samples, Distributions, Surrogate, T <: Real}
        return new{S, Samples, Distributions, Surrogate, T}(
            statistics,
            samples,
            distributions,
            surrogate,
            ntrials,
            confidence
        )
    end
end

function LineParametersMC(
        statistics::RLCG{S},
        samples::Samples,
        distributions::Distributions,
        surrogate::Surrogate,
        ntrials::Integer,
        confidence::T
) where {S, Samples, Distributions, Surrogate, T <: Real}
    ntrials > 0 || throw(ArgumentError("ntrials must be positive"))
    0 < confidence < 1 || throw(ArgumentError("confidence must lie between zero and one"))
    all(value -> value isa AbstractArray{<:SampleSummary, 3},
        (
            statistics.R,
            statistics.L,
            statistics.C,
            statistics.G
        )) ||
        throw(ArgumentError("R, L, C, and G statistics must be three-dimensional SampleSummary arrays"))
    statistic_size = size(statistics.R)
    all(size(values) == statistic_size
    for values in (statistics.L, statistics.C, statistics.G)) ||
        throw(DimensionMismatch("R, L, C, and G statistics must have equal dimensions"))
    if samples !== nothing
        samples isa RLCG || throw(
            ArgumentError("line-parameter samples must be RLCG or nothing"),
        )
        all(value -> value isa AbstractArray{<:Real, 4}, (
            samples.R,
            samples.L,
            samples.C,
            samples.G
        )) ||
            throw(ArgumentError("retained R, L, C, and G samples must be real four-dimensional arrays"))
        sample_size = size(samples.R)
        length(sample_size) == 4 ||
            throw(DimensionMismatch("retained line-parameter samples must be four-dimensional"))
        sample_size[4] == ntrials || throw(
            DimensionMismatch("retained line-parameter samples must match ntrials"),
        )
        all(size(values) == sample_size for values in (samples.L, samples.C, samples.G)) ||
            throw(DimensionMismatch("R, L, C, and G samples must have equal dimensions"))
        sample_size[1:3] == statistic_size || throw(
            DimensionMismatch("sample and statistic dimensions must agree"),
        )
    end
    if distributions !== nothing
        distributions isa RLCG || throw(
            ArgumentError("line-parameter distributions must be RLCG or nothing"),
        )
        all(value -> value isa AbstractArray{<:HistogramPDF, 3},
            (
                distributions.R,
                distributions.L,
                distributions.C,
                distributions.G
            )) ||
            throw(ArgumentError("R, L, C, and G distributions must be three-dimensional HistogramPDF arrays"))
        all(size(values) == statistic_size
        for values in (
            distributions.R,
            distributions.L,
            distributions.C,
            distributions.G
        )) || throw(DimensionMismatch("distribution and statistic dimensions must agree"))
    end
    surrogate isa LineParameters || throw(
        ArgumentError("the line-parameter surrogate must be a LineParameters value"),
    )
    size(surrogate.Z) == statistic_size || throw(
        DimensionMismatch("surrogate and statistic dimensions must agree"),
    )
    return LineParametersMC{S, Samples, Distributions, Surrogate, T}(
        statistics,
        samples,
        distributions,
        surrogate,
        Int(ntrials),
        confidence
    )
end

const _CABLE_CONSTANT_QUANTITIES = (:R, :L, :C)
const _LINE_PARAMETER_QUANTITIES = (:R, :L, :C, :G)

@inline function _quantity(container, quantity::Symbol, supported)
    quantity in supported || throw(
        ArgumentError("unsupported quantity :$quantity; expected one of $(supported)"),
    )
    return getproperty(container, quantity)
end

@inline _selection(value, ::Tuple{}) = value
@inline _selection(value::AbstractArray{<:Any, 3}, indices::Tuple{Int, Int}) = view(
    value, indices[1], indices[2], :)
@inline _selection(value::AbstractArray{<:Any, 3}, indices::Tuple{
    Int, Int, Any}) = value[indices...]
@inline _selection(value::AbstractArray{<:Any, 4}, indices::Tuple{Int, Int}) = view(
    value, indices[1], indices[2], :, :)
@inline _selection(value::AbstractArray{<:Any, 4}, indices::Tuple{Int, Int, Any}) = view(
    value, indices[1], indices[2], indices[3], :)
@inline _selection(value, indices::Tuple) = getindex(value, indices...)

function statistics(result::CableConstantsMC, quantity::Symbol, indices...)
    _selection(_quantity(result.statistics, quantity, _CABLE_CONSTANT_QUANTITIES), indices)
end
function statistics(result::LineParametersMC, quantity::Symbol, indices...)
    _selection(_quantity(result.statistics, quantity, _LINE_PARAMETER_QUANTITIES), indices)
end
statistics(result::Union{CableConstantsMC, LineParametersMC}) = result.statistics

@inline _summary_field(summary::SampleSummary, field::Symbol) = getproperty(summary, field)
function _summary_field(summaries::AbstractArray, field::Symbol)
    map(summary -> getproperty(summary, field), summaries)
end

function Statistics.mean(result::Union{CableConstantsMC, LineParametersMC}, quantity::Symbol, indices...)
    _summary_field(statistics(result, quantity, indices...), :mean)
end
function Statistics.std(result::Union{CableConstantsMC, LineParametersMC}, quantity::Symbol, indices...)
    _summary_field(statistics(result, quantity, indices...), :std)
end

has_samples(result::Union{CableConstantsMC, LineParametersMC}) = result.samples !== nothing
function has_distributions(result::Union{CableConstantsMC, LineParametersMC})
    result.distributions !== nothing
end

function samples(result::CableConstantsMC, quantity::Symbol, indices...)
    has_samples(result) || throw(
        ArgumentError("samples were not retained; rerun mc(...; return_samples=true)"),
    )
    return _selection(_quantity(result.samples, quantity, _CABLE_CONSTANT_QUANTITIES), indices)
end

function samples(result::LineParametersMC, quantity::Symbol, indices...)
    has_samples(result) || throw(
        ArgumentError("samples were not retained; rerun mc(...; return_samples=true)"),
    )
    return _selection(_quantity(result.samples, quantity, _LINE_PARAMETER_QUANTITIES), indices)
end

function distribution(result::CableConstantsMC, quantity::Symbol, indices...)
    has_distributions(result) || throw(
        ArgumentError("distributions were not retained; rerun mc(...; return_pdf=true)"),
    )
    return _selection(
        _quantity(result.distributions, quantity, _CABLE_CONSTANT_QUANTITIES),
        indices
    )
end

function distribution(result::LineParametersMC, quantity::Symbol, indices...)
    has_distributions(result) || throw(
        ArgumentError("distributions were not retained; rerun mc(...; return_pdf=true)"),
    )
    return _selection(
        _quantity(result.distributions, quantity, _LINE_PARAMETER_QUANTITIES),
        indices
    )
end

surrogate(result::Union{CableConstantsMC, LineParametersMC}) = result.surrogate
basis(::CableConstantsMC) = :per_length
ntrials(result::Union{CableConstantsMC, LineParametersMC}) = result.ntrials
confidence(result::Union{CableConstantsMC, LineParametersMC}) = result.confidence
frequencies(result::LineParametersMC) = frequencies(result.surrogate)
function domain(
        ::LineParametersMC{S,
        Samples,
        Distributions,
        Surrogate},
) where {
        S,
        Samples,
        Distributions,
        T,
        U,
        D,
        Basis,
        Surrogate <: LineParameters{T, U, D, Basis}
}
    D
end
function basis(
        ::LineParametersMC{S,
        Samples,
        Distributions,
        Surrogate},
) where {
        S,
        Samples,
        Distributions,
        T,
        U,
        D,
        Basis,
        Surrogate <: LineParameters{T, U, D, Basis}
}
    Basis
end
nconductors(result::LineParametersMC) = nconductors(result.surrogate)
nfrequencies(result::LineParametersMC) = nfrequencies(result.surrogate)

function _fixed_quantile(summary::SampleSummary, probability::Real)
    probability == 0.05 && return summary.q05
    probability == 0.50 && return summary.q50
    probability == 0.95 && return summary.q95
    return nothing
end

function StatsBase.quantile(
        result::Union{CableConstantsMC, LineParametersMC},
        quantity::Symbol,
        probability::Real,
        indices...
)
    0 <= probability <= 1 || throw(ArgumentError("probability must lie in [0, 1]"))
    summary = statistics(result, quantity, indices...)
    if summary isa SampleSummary
        fixed = _fixed_quantile(summary, probability)
        fixed !== nothing && return fixed
    elseif probability in (0.05, 0.50, 0.95)
        field = probability == 0.05 ? :q05 : probability == 0.50 ? :q50 : :q95
        return _summary_field(summary, field)
    end

    if has_samples(result)
        values = samples(result, quantity, indices...)
        values isa AbstractVector || throw(
            ArgumentError("select one scalar observable before requesting an arbitrary quantile"),
        )
        return StatsBase.quantile(values, probability)
    elseif has_distributions(result)
        value = distribution(result, quantity, indices...)
        value isa HistogramPDF || throw(
            ArgumentError("select one scalar distribution before requesting an arbitrary quantile"),
        )
        return Distributions.quantile(value, probability)
    end
    throw(
        ArgumentError(
        "only 0.05, 0.50, and 0.95 are stored; retain samples or distributions for other probabilities",
    ),
    )
end
