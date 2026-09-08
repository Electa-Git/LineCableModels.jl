_uq_compact(value) = sprint(show, value; context = :compact => true, sizehint = 80)

TextDisplay.name(::Type{<:LinearError}) = "LinearError"
Base.summary(io::IO, ::LinearError) = print(io, "Linear uncertainty propagation")
function Base.show(io::IO, formulation::LinearError)
    print(io, "LinearError(", _uq_compact(formulation.inner), ")")
end
function Base.show(io::IO, ::MIME"text/plain", formulation::LinearError)
    get(io, :compact, false) && return show(io, formulation)
    return TextDisplay.tree(io, "Linear uncertainty propagation", (
        (label = "inner    $(_uq_compact(formulation.inner))", noun = "fields"),
        (label = "details  $(formulation.options.retain_details ? "retained" : "not retained")", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:MonteCarlo}) = "MonteCarlo"
Base.summary(io::IO, formulation::MonteCarlo) = print(io, "Monte Carlo propagation")
function Base.show(io::IO, formulation::MonteCarlo)
    trials = formulation.trials === nothing ? "DKW-sized" : string(formulation.trials)
    print(io, "MonteCarlo(trials=", trials,
        "; samples=", formulation.return_samples,
        ", histograms=", formulation.return_histograms, ")")
end

TextDisplay.name(::Type{<:PolynomialChaos}) = "PolynomialChaos"
Base.summary(io::IO, ::PolynomialChaos) = print(io, "Polynomial chaos propagation")
function Base.show(io::IO, formulation::PolynomialChaos)
    print(io, "PolynomialChaos(degree=", formulation.degree,
        "; distribution=:", formulation.distribution, ")")
end
function Base.show(io::IO, ::MIME"text/plain", formulation::PolynomialChaos)
    get(io, :compact, false) && return show(io, formulation)
    return TextDisplay.tree(io, "Polynomial chaos propagation", (
        (label = "inner         $(_uq_compact(formulation.inner))", noun = "fields"),
        (label = "degree        $(formulation.degree)", noun = "fields"),
        (label = "quadrature    $(formulation.quadrature_order)", noun = "fields"),
        (label = "distribution  :$(formulation.distribution)", noun = "fields"),
        (label = "validation    $(formulation.validation_points) points", noun = "fields"),
        (label = "tolerance     $(TextDisplay.value(formulation.validation_rtol))", noun = "fields"),
    ))
end
function Base.show(io::IO, ::MIME"text/plain", formulation::MonteCarlo)
    get(io, :compact, false) && return show(io, formulation)
    trials = formulation.trials === nothing ? "DKW-sized" : string(formulation.trials)
    distribution = formulation.distribution isa Symbol ?
                   ":$(formulation.distribution)" : _uq_compact(formulation.distribution)
    children = Any[
        (label = "inner        $(_uq_compact(formulation.inner))", noun = "fields"),
        (label = "trials       $trials", noun = "fields"),
        (label = "confidence   $(TextDisplay.value(formulation.confidence))", noun = "fields"),
        (label = "CDF tol      $(TextDisplay.value(formulation.cdf_tol))", noun = "fields"),
        (label = "distribution $distribution", noun = "fields"),
    ]
    formulation.seed === nothing || push!(children,
        (label = "seed         $(formulation.seed)", noun = "fields"))
    formulation.return_samples && push!(children,
        (label = "samples      retained", noun = "fields"))
    formulation.return_histograms && push!(children,
        (label = "histograms   retained", noun = "fields"))
    return TextDisplay.tree(io, "Monte Carlo propagation", Tuple(children))
end

TextDisplay.@showfields SampleSummary "SampleSummary" summary -> (
    mean = TextDisplay.value(summary.mean),
    std = TextDisplay.value(summary.std),
    min = TextDisplay.value(summary.min),
    q05 = TextDisplay.value(summary.q05),
    median = TextDisplay.value(summary.median),
    q95 = TextDisplay.value(summary.q95),
    max = TextDisplay.value(summary.max),
    n = summary.n
)

TextDisplay.name(::Type{<:HistogramDensity}) = "HistogramDensity"
function Base.summary(io::IO, histogram::HistogramDensity)
    print(io, "Histogram density with $(length(histogram.density)) bins")
end
function Base.show(io::IO, histogram::HistogramDensity)
    print(io, "HistogramDensity(", length(histogram.density), " bins; range=",
        TextDisplay.value(first(histogram.edges)), "…",
        TextDisplay.value(last(histogram.edges)), ")")
end
function Base.show(io::IO, ::MIME"text/plain", histogram::HistogramDensity)
    get(io, :compact, false) && return show(io, histogram)
    return TextDisplay.tree(io, "Histogram density · $(length(histogram.density)) bins", (
        (label = "range    $(TextDisplay.value(first(histogram.edges))) … $(TextDisplay.value(last(histogram.edges)))", noun = "fields"),
        (label = "density  normalized", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:LinearErrorResult}) = "LinearErrorResult"
function Base.summary(io::IO, result::LinearErrorResult)
    print(io, "Linear-error result with $(length(result)) values")
end
Base.show(io::IO, result::LinearErrorResult) =
    print(io, "LinearErrorResult($(length(result)) values)")
function Base.show(io::IO, ::MIME"text/plain", result::LinearErrorResult)
    get(io, :compact, false) && return show(io, result)
    kind = isempty(result.values) ? "none" : String(nameof(eltype(result.values)))
    return TextDisplay.tree(io, "Linear-error result · $(length(result)) values", (
        (label = "result type  $kind", noun = "fields"),
        (label = "details      $(length(result.details)) entries", noun = "fields"),
    ))
end

TextDisplay.name(::Type{<:MonteCarloResult}) = "MonteCarloResult"
function Base.summary(io::IO, result::MonteCarloResult)
    print(io, "Monte Carlo result with $(length(result)) points")
end

TextDisplay.name(::Type{<:PolynomialChaosResult}) = "PolynomialChaosResult"
function Base.summary(io::IO, result::PolynomialChaosResult)
    print(io, "Polynomial chaos result with $(length(result)) points")
end
function Base.show(io::IO, result::PolynomialChaosResult)
    print(io, "PolynomialChaosResult($(length(result)) points; degree=",
        result.formulation.degree, ")")
end
function Base.show(io::IO, ::MIME"text/plain", result::PolynomialChaosResult)
    get(io, :compact, false) && return show(io, result)
    result_type = String(nameof(eltype(result.values)))
    children = Any[
        (label = "result type  $result_type", noun = "fields"),
        (label = "degree       $(result.formulation.degree)", noun = "fields"),
        (label = "distribution :$(result.formulation.distribution)", noun = "fields"),
        (label = "expansions   retained", noun = "fields"),
        (label = "statistics   retained", noun = "fields"),
        (label = "validation   passed", noun = "fields"),
    ]
    isempty(result.details) || push!(children,
        (label = "details      $(length(result.details.points)) points", noun = "fields"))
    return TextDisplay.tree(
        io, "Polynomial chaos result · $(length(result)) points", Tuple(children))
end
function Base.show(io::IO, result::MonteCarloResult)
    trials = all(==(first(result.trial_counts)), result.trial_counts) ?
             string(first(result.trial_counts)) : "varying"
    print(io, "MonteCarloResult($(length(result)) points; trials=", trials, ")")
end
function Base.show(io::IO, ::MIME"text/plain", result::MonteCarloResult)
    get(io, :compact, false) && return show(io, result)
    trial_span = extrema(result.trial_counts)
    trials = first(trial_span) == last(trial_span) ? string(first(trial_span)) :
             "$(first(trial_span))…$(last(trial_span))"
    result_type = String(nameof(eltype(result.values)))
    children = Any[
        (label = "result type  $result_type", noun = "fields"),
        (label = "trials       $trials per point", noun = "fields"),
        (label = "statistics   retained", noun = "fields"),
        (label = "root seed    $(result.root_seed)", noun = "fields"),
    ]
    result.sample_values === nothing || push!(children,
        (label = "samples      retained", noun = "fields"))
    result.histogram_values === nothing || push!(children,
        (label = "histograms   retained", noun = "fields"))
    isempty(result.details) || push!(children,
        (label = "details      $(length(result.details)) entries", noun = "fields"))
    return TextDisplay.tree(io, "Monte Carlo result · $(length(result)) points", Tuple(children))
end
