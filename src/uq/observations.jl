"""
Return core results in Gridspace traversal order.
"""
Units.label(::Units.Quantity{:sample_count}) = "Count"
Units.symbol(::Units.Quantity{:sample_count}) = "n"
Units.label(::Units.Quantity{:probability}) = "Probability"
Units.symbol(::Units.Quantity{:probability}) = "p"
Units.label(::Units.Quantity{:cumulative_probability}) = "Cumulative probability"
Units.symbol(::Units.Quantity{:cumulative_probability}) = "F"
Units.label(::Units.Quantity{:probability_density}) = "Probability density"
Units.symbol(::Units.Quantity{:probability_density}) = "p"

const _DimensionlessStatisticalQuantity = Union{
    Units.Quantity{:sample_count},
    Units.Quantity{:probability},
    Units.Quantity{:cumulative_probability}
}
Units.native_unit(::_DimensionlessStatisticalQuantity) = Units.units(:base, :dimensionless)
Units.display_unit(::_DimensionlessStatisticalQuantity) = Units.units(:base, :dimensionless)

"""
Return the sample-summary product for each Monte Carlo point.
"""
statistics(value::MonteCarloResult) = value.stats

"""
Return retained sample products, or `nothing` when retention was disabled.
"""
samples(value::MonteCarloResult) = value.sample_values

"""
Return retained histogram products, or `nothing` when retention was disabled.
"""
histograms(value::MonteCarloResult) = value.histogram_values

basis(value::MonteCarloResult) = basis(first(value.values))
basis(value::LinearErrorResult) = basis(first(value.values))

"""
$(TYPEDSIGNATURES)

Return uncertainty-bearing core results, in configuration order, or the core
result at `point`. Monte Carlo reconstruction requires Measurements.jl and
preserves retained marginal means and standard deviations, not joint output
correlations. An indexed Monte Carlo read reconstructs only the selected point.
"""
uncertain(value::LinearErrorResult) = value.values
uncertain(value::AbstractUncertaintyResult, point::Integer) = uncertain(value)[point]

function uncertain(value::MonteCarloResult)
    throw(ArgumentError("uncertain requires a reconstruction for this Monte Carlo result; " *
        "load Measurements.jl for built-in cable and line results"))
end

"""
$(TYPEDSIGNATURES)

Return the resolved root random seed of a Monte Carlo calculation.
"""
root_seed(value::MonteCarloResult) = value.root_seed

"""
$(TYPEDSIGNATURES)

Return the random seed used for the Gridspace point at `index`.
"""
point_seed(value::MonteCarloResult, index::Integer) = value.point_seeds[index]

"""
$(TYPEDSIGNATURES)

Return the number of trials evaluated for the Gridspace point at `index`.
"""
trial_count(value::MonteCarloResult, index::Integer) = value.trial_counts[index]

"""
$(TYPEDSIGNATURES)

Return the simultaneous empirical-CDF confidence of a Monte Carlo calculation.
"""
confidence(value::MonteCarloResult) = value.formulation.options.confidence

"""
$(TYPEDSIGNATURES)

Return the empirical-CDF tolerance used to size a Monte Carlo calculation.
"""
cdf_tolerance(value::MonteCarloResult) = value.formulation.options.cdf_tol

"""
$(TYPEDSIGNATURES)

Return the sampling distribution of a Monte Carlo calculation.
"""
sampling_distribution(value::MonteCarloResult) = value.formulation.options.distribution

"""
$(TYPEDSIGNATURES)

Return the simultaneous DKW bound for all retained real marginal summaries at
one outer point. Both matrix orientations are counted conservatively. This is
an iid-sampling bound, not measured CDF error; retry sampling concerns the
population conditional on success. Fixed trial counts do not certify the
configured target automatically.
"""
function confidence(value::MonteCarloResult, point::Integer)
    count = sum(length, values(value.stats[point]))
    trials = trial_count(value,point)
    bound = sqrt(log(2count / (1-confidence(value))) / (2trials))
    retained=details(value)
    diagnostics = isempty(retained) ? nothing :
        (failures=retained.failures[point],failure_summary=retained.failure_summary[point],
            clearance=haskey(retained,:clearance) ? retained.clearance[point] : nothing)
    mean_standard_error=map(value.stats[point]) do summaries
        map(summary -> summary.n>1 ? summary.std/sqrt(summary.n) : missing,summaries)
    end
    return (point, trials, spread_estimated=trials>1,marginal_count=count, confidence=confidence(value),
        target_cdf=cdf_tolerance(value), cdf_bound=min(1.0,bound),
        target_supported=bound <= cdf_tolerance(value), scope=:point_all_retained_marginals,
        assumption=:iid, distribution=sampling_distribution(value),
        conditioning=value.formulation.options.on_error === :retry ? :successful_realizations : :none,
        diagnostics, mean_standard_error, root_seed=root_seed(value), point_seed=point_seed(value,point),
        samples_retained=samples(value) !== nothing, histograms_retained=histograms(value) !== nothing)
end

const _MonteCarloProductSelector = Union{
    typeof(statistics),
    typeof(samples),
    typeof(histograms)
}

const _MonteCarloScientificSelector = Union{
    typeof(R),
    typeof(L),
    typeof(C),
    typeof(Engine.G), typeof(Engine.X), typeof(Engine.B), typeof(Engine.Z), typeof(Engine.Y)
}

const _StatisticSelector = Union{
    typeof(Statistics.mean),
    typeof(Statistics.std),
    typeof(Statistics.median),
    typeof(minimum),
    typeof(maximum), Base.Fix2{typeof(Statistics.quantile)}
}

function Units.quantity(
        ::_MonteCarloProductSelector,
        selector::_MonteCarloScientificSelector
)
    return Units.quantity(selector)
end

function Units.quantity(::typeof(statistics), selector::_MonteCarloScientificSelector,
        ::_StatisticSelector)
    return Units.quantity(selector)
end

function _monte_carlo_product(value::MonteCarloResult, ::typeof(statistics), point::Integer)
    return value.stats[point]
end

function _monte_carlo_product(value::MonteCarloResult, ::typeof(samples), point::Integer)
    retained = value.sample_values
    retained === nothing && throw(ArgumentError("Monte Carlo samples were not retained"))
    return retained[point]
end

function _monte_carlo_product(value::MonteCarloResult, ::typeof(histograms), point::Integer)
    retained = value.histogram_values
    retained === nothing && throw(ArgumentError(
        "Monte Carlo histograms were not retained or derived",
    ))
    return retained[point]
end

function observe(
        value::MonteCarloResult,
        selector::_MonteCarloScientificSelector,
        point::Integer,
        indices...
)
    return observe(value.values[point], selector, indices...)
end

function observe(
        value::MonteCarloResult{<:Engine.LineParameters},
        ::typeof(frequencies),
        point::Integer,
        indices...
)
    return observe(value.values[point], frequencies, indices...)
end

function observe(
        value::MonteCarloResult,
        product::_MonteCarloProductSelector,
        selector::_MonteCarloScientificSelector,
        point::Integer,
        indices...
)
    stored = if selector in (Engine.X, Engine.B)
        source = selector === Engine.X ? L : C
        values = observe(value, product, source, point)
        angular = reshape(2pi .* frequencies(value[point]), 1, 1, :)
        product === samples && (angular = reshape(angular, 1, 1, :, 1))
        detach.(values, angular)
    elseif selector in (Engine.Z, Engine.Y)
        product === samples || throw(ArgumentError(
            "complex quantities support selected mean/std or retained joint samples, not ordered summaries"))
        real_selector, imaginary_selector = selector === Engine.Z ? (R, Engine.X) : (Engine.G, Engine.B)
        observe(value, samples, real_selector, point) .+
            im .* observe(value, samples, imaginary_selector, point)
    else
        _monte_carlo_field(_monte_carlo_product(value, product, point), selector)
    end
    return _product_value(stored, indices)
end

function _histogram_observation(
        value::MonteCarloResult,
        selector::_MonteCarloScientificSelector,
        point::Integer,
        indices::Tuple,
        bins::Union{Nothing, Integer}
)
    bins === nothing || bins > 0 || throw(ArgumentError("histogram bins must be positive"))
    if selector in (Engine.X,Engine.B)
        base_selector=selector===Engine.X ? L : C
        original=_histogram_observation(value,base_selector,point,indices,bins)
        return detach(original,2pi*frequencies(value[point])[last(indices)])
    end
    selector in (Engine.Z,Engine.Y) && throw(ArgumentError("complex histograms require a real-valued observable"))
    if value.histogram_values !== nothing
        stored = _monte_carlo_field(value.histogram_values[point], selector)
        histogram = _product_value(stored, indices)
        (bins === nothing || length(histogram.density) == bins) && return histogram
        value.sample_values === nothing && throw(ArgumentError(
            "Changing histogram bins requires retained samples; use " *
            "MonteCarlo(...; return_samples=true), or keep the retained bin count",
        ))
    end
    retained = value.sample_values
    retained === nothing && throw(ArgumentError(
        "Monte Carlo histograms were not retained and samples are unavailable for derivation",
    ))
    stored = _monte_carlo_field(retained[point], selector)
    sample = _product_value(stored, (indices..., Colon()))
    return HistogramDensity(collect(sample); bins)
end

"""
$(TYPEDSIGNATURES)

Observe a cable-constant marginal as a [`HistogramDensity`](@ref), using
retained samples when a new binning is requested. No simulation is performed
and retained products are not modified.

# Arguments

- `value`: Monte Carlo result.
- `histograms`: Histogram product selector.
- `selector`: `R`, `L`, `C`, or `G`.
- `point`: Gridspace point index.
- `assembly`: Cable assembly index.
- `bins`: Positive bin count, or `nothing` to reuse the retained model. Without
  a retained model, `nothing` selects the sample-based automatic bin count.

# Returns

- A normalized histogram in the observed quantity's native units. Constant
  samples produce one finite-width bin, irrespective of the requested count.

# Errors

An `ArgumentError` is raised for nonpositive counts or when deriving a new
histogram requires samples that were not retained. A retained model can be
reused without samples when its bin count matches the request.
"""
function observe(
        value::MonteCarloResult{<:Engine.CableConstants},
        ::typeof(histograms),
        selector::_MonteCarloScientificSelector,
        point::Integer,
        assembly::Integer,
        bins::Union{Nothing, Integer}
)
    return _histogram_observation(value, selector, point, (assembly,), bins)
end

"""
$(TYPEDSIGNATURES)

Observe a line-parameter histogram at one matrix element and frequency index.
The `bins` argument follows the same retention and derivation rules as the
cable-constant histogram observation; `row`, `column`, and `frequency` select
the marginal in place of `assembly`.
"""
function observe(
        value::MonteCarloResult{<:Engine.LineParameters},
        ::typeof(histograms),
        selector::_MonteCarloScientificSelector,
        point::Integer,
        row::Integer,
        column::Integer,
        frequency::Integer,
        bins::Union{Nothing, Integer}
)
    return _histogram_observation(
        value,
        selector,
        point,
        (row, column, frequency),
        bins
    )
end

function observe(
        value::MonteCarloResult,
        ::typeof(statistics),
        selector::_MonteCarloScientificSelector,
        transform::_StatisticSelector,
        point::Integer,
        indices...
)
    if selector in (Engine.X,Engine.B)
        base_selector=selector===Engine.X ? L : C
        values=observe(value,statistics,base_selector,transform,point,indices...)
        sample=length(indices)==3 ? last(indices) : Colon()
        angular=2pi .* frequencies(value[point])[sample]
        factor=angular isa AbstractArray ? reshape(angular,ntuple(_ -> 1,ndims(values)-1)...,:) : angular
        return values .* factor
    end
    if selector in (Engine.Z, Engine.Y)
        transform in (Statistics.mean, Statistics.std) || throw(ArgumentError(
            "complex statistics require mean or std; ordered statistics need a real-valued observable"))
        real_selector, imaginary_selector = selector === Engine.Z ? (R, Engine.X) : (Engine.G, Engine.B)
        real_values = observe(value, statistics, real_selector, transform, point, indices...)
        imaginary_values = observe(value, statistics, imaginary_selector, transform, point, indices...)
        return transform === Statistics.mean ? real_values .+ im .* imaginary_values :
            hypot.(real_values, imaginary_values)
    end
    return _statistic(transform, observe(value, statistics, selector, point, indices...))
end

function _monte_carlo_observables(selectors::Tuple)
    product_selectors = (statistics, samples, histograms)
    real_selectors = filter(selector -> !(selector in (Engine.Z, Engine.Y)), selectors)
    products = Tuple((product, selector)
        for product in product_selectors for selector in real_selectors)
    selected = Tuple((statistics, selector, transform)
        for selector in selectors for transform in
            (selector in (Engine.Z, Engine.Y) ? (Statistics.mean, Statistics.std) :
             (Statistics.mean, Statistics.std, minimum, Statistics.median, maximum,
              Base.Fix2(Statistics.quantile, 0.0), Base.Fix2(Statistics.quantile, 0.05),
              Base.Fix2(Statistics.quantile, 0.5), Base.Fix2(Statistics.quantile, 0.95),
              Base.Fix2(Statistics.quantile, 1.0))))
    complex_samples=Tuple((samples,selector) for selector in selectors if selector in (Engine.Z,Engine.Y))
    return (selectors..., products..., selected..., complex_samples...)
end

function observables(
        ::Type{<:MonteCarloResult{T}}
) where {T <: Engine.CableConstants}
    return _monte_carlo_observables((R, L, C, Engine.G))
end

function observables(
        ::Type{<:MonteCarloResult{T}}
) where {T <: Engine.LineParameters}
    return (frequencies, _monte_carlo_observables((R, L, C, Engine.G, Engine.X, Engine.B, Engine.Z, Engine.Y))...)
end

function observables(::Type{<:LinearErrorResult{T}}) where {T}
    selectors = T <: Engine.CableConstants ? (R, L, C, Engine.G) :
        (R, L, C, Engine.G, Engine.X, Engine.B, Engine.Z, Engine.Y)
    selected = Tuple((statistics, selector, transform) for selector in selectors
        for transform in (Statistics.mean, Statistics.std))
    return (selectors..., selected..., (T <: Engine.LineParameters ? (frequencies,) : ())...)
end

function observe(value::LinearErrorResult, selector::_MonteCarloScientificSelector,
        point::Integer, indices...)
    return observe(value[point], selector, indices...)
end

function observe(value::LinearErrorResult{<:Engine.LineParameters}, ::typeof(frequencies),
        point::Integer, indices...)
    return observe(value[point], frequencies, indices...)
end

"""
$(TYPEDSIGNATURES)

Observe first-order nominal values or propagated standard uncertainties in the
quantity's native units. Complex standard deviation is the nonnegative root
sum of component variances; it is not a magnitude-distribution statistic.
No output distribution or independent Measurement values are constructed.
"""
function observe(value::LinearErrorResult, ::typeof(statistics),
        selector::_MonteCarloScientificSelector,
        transform::Union{typeof(Statistics.mean),typeof(Statistics.std)},
        point::Integer, indices...)
    values = observe(value[point], selector, indices...)
    if transform === Statistics.mean
        return nominal.(values)
    end
    return hypot.(uncertainty.(real.(values)), uncertainty.(imag.(values)))
end

function observation_resolution(source::Union{MonteCarloResult,LinearErrorResult}, request;
        atol=nothing, frequencies=nothing)
    identity = request_identity(request)
    if !(identity isa Tuple && length(identity) == 3 && first(identity) === statistics)
        return observation_resolution(nothing, request; atol, frequencies)
    end
    point, indices = _statistics_point(request)
    core = source[point]
    core isa Engine.LineParameters || return observation_resolution(nothing, request; atol, frequencies)
    f = observe(source, LineCableModels.frequencies, point)
    frequencies === nothing || frequencies == f || throw(ArgumentError("supplied frequencies differ from the UQ point"))
    sample = length(indices) == 3 ? last(indices) : Colon()
    values = observe(source, request...)
    return observation_resolution(values, identity[2]; atol, frequencies=f[sample], result_basis=basis(core))
end

@inline _product_value(value, ::Tuple{}) = value
@inline _product_value(value, indices::Tuple) = getindex(value, indices...)

@inline _monte_carlo_field(product, ::typeof(R)) = product.R
@inline _monte_carlo_field(product, ::typeof(L)) = product.L
@inline _monte_carlo_field(product, ::typeof(C)) = product.C
@inline _monte_carlo_field(product, ::typeof(Engine.G)) = product.G

@inline _statistic(transform, value::AbstractArray) = map(transform, value)
@inline _statistic(transform, value) = transform(value)

function detach(summary::SampleSummary, factor)
    return SampleSummary(
        summary.mean * factor,
        summary.std * abs(factor),
        summary.min * factor,
        summary.q05 * factor,
        summary.median * factor,
        summary.q95 * factor,
        summary.max * factor,
        summary.n
    )
end

function detach(summary::SampleSummary, factor, clip::Bool)
    return detach(summary, factor)
end

function detach(
        summaries::AbstractArray{<:SampleSummary},
        factor
)
    return map(summary -> detach(summary, factor), summaries)
end

function detach(
        summaries::AbstractArray{<:SampleSummary},
        factor,
        clip::Bool
)
    return map(summary -> detach(summary, factor, clip), summaries)
end

function detach(histogram::HistogramDensity, factor)
    factor > zero(factor) || throw(ArgumentError("histogram conversion must be positive"))
    return HistogramDensity(
        histogram.edges .* factor,
        histogram.density ./ factor
    )
end
