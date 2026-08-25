"""
$(TYPEDEF)

PlotBuilder recipe for marginal Monte Carlo histograms, probability densities,
cumulative distributions, and Q-Q plots.
"""
struct MCDistributionPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

function _mc_plot_exponent(series, field::Symbol)
    maximum_value = 0.0
    for item in series
        values = field === :x ? item.xdata : item.ydata
        values === nothing && continue
        for sample in values
            sample isa Real || continue
            nominal_value = abs(Float64(sample))
            isfinite(nominal_value) &&
                (maximum_value = max(maximum_value, nominal_value))
        end
    end
    iszero(maximum_value) && return 0
    exponent = floor(Int, log10(maximum_value))
    return abs(exponent) < 3 ? 0 : exponent
end

function _mc_selection(::DataModel.CableConstants, statistics, selector::Function, ijk)
    selector in (R, L, C) || throw(
        ArgumentError("cable-constant selectors are R, L, and C"),
    )
    ijk === nothing || throw(
        ArgumentError("cable-constant Monte Carlo results do not use matrix indices"),
    )
    return nothing
end

function _mc_selection(
        ::Engine.LineParameters,
        statistics,
        selector::Function,
        ijk
)
    selector in (R, L, C, Engine.G) || throw(
        ArgumentError("line-parameter selectors are R, L, C, and G"),
    )
    selection = ijk === nothing ? (1, 1, 1) : ijk
    selection isa NTuple{3, Int} || throw(
        ArgumentError("ijk must be a tuple (i, j, k)"),
    )
    observable = observe(statistics, selector)
    checkbounds(observable, selection...)
    return selection
end

function _mc_request(selector, selection, retained_samples::Bool)
    selection === nothing && return selector
    return retained_samples ? (selector, selection..., Colon()) : (selector, selection...)
end

function _mc_quantity_prefix(quantity_units, fallback::Symbol)
    quantity_units === nothing && return fallback
    quantity_units isa Symbol && return quantity_units
    quantity_units isa Units.UnitExpr && return quantity_units
    throw(ArgumentError("quantity_units must be a prefix, UnitExpr, or nothing"))
end

function _mc_target_unit(product, selector::Function, length_unit, quantity_units)
    tag = Units.quantity(selector)
    default = Units.display_unit(tag, basis(product); length_prefix = length_unit)
    prefix = _mc_quantity_prefix(
        quantity_units,
        first(default.numerator).prefix
    )
    prefix isa Units.UnitExpr && return tag, prefix
    prefix isa Symbol || throw(
        ArgumentError("quantity-unit overrides must be prefixes or UnitExpr values"),
    )
    target = Units.display_unit(
        tag,
        basis(product);
        length_prefix = length_unit,
        prefix
    )
    return tag, target
end

function _mc_histogram_cdf(histogram::HistogramDensity, x::Real)
    x <= first(histogram.edges) && return 0.0
    x >= last(histogram.edges) && return 1.0
    bin = clamp(searchsortedlast(histogram.edges, x), 1, length(histogram.density))
    cumulative = bin == 1 ? 0.0 :
                 sum(
        histogram.density[1:(bin - 1)] .* diff(histogram.edges[1:bin]),
    )
    return clamp(
        cumulative + histogram.density[bin] * (x - histogram.edges[bin]),
        0.0,
        1.0
    )
end

function _mc_histogram_quantile(histogram::HistogramDensity, probability::Real)
    0 <= probability <= 1 || throw(
        DomainError(probability, "probability must lie between zero and one"),
    )
    probability == 0 && return first(histogram.edges)
    probability == 1 && return last(histogram.edges)
    masses = histogram.density .* diff(histogram.edges)
    cumulative = cumsum(masses)
    bin = something(findfirst(>=(probability), cumulative), length(masses))
    previous = bin == 1 ? 0.0 : cumulative[bin - 1]
    density = histogram.density[bin]
    iszero(density) && return histogram.edges[bin]
    return histogram.edges[bin] + (probability - previous) / density
end

function _mc_empirical_cdf(sorted_values, x)
    searchsortedlast(sorted_values, x) /
    length(sorted_values)
end

function _mc_input_defaults()
    return (;
        selector = R,
        ijk = nothing,
        mode = :hist,
        data = :samples,
        length_unit = :kilo,
        quantity_units = nothing,
        nbins = nothing,
        normalization = :none
    )
end

PlotBuilder.dispatch_on(::Type{MCDistributionPlotDefinition}) = MonteCarloResult
function PlotBuilder.input_kwargs(::Type{MCDistributionPlotDefinition})
    return (
        :selector,
        :ijk,
        :mode,
        :data,
        :length_unit,
        :quantity_units,
        :nbins,
        :normalization
    )
end
PlotBuilder.renderer_kwargs(::Type{MCDistributionPlotDefinition}) = (:fig_size,)
function PlotBuilder.input_defaults(::Type{MCDistributionPlotDefinition}, ::MonteCarloResult)
    _mc_input_defaults()
end
function PlotBuilder.renderer_defaults(
        ::Type{MCDistributionPlotDefinition},
        ::MonteCarloResult
)
    return (; fig_size = (800, 400))
end

function _mc_required_data(mode::Symbol, data::Symbol)
    needs_samples = mode === :qq ||
                    (mode in (:hist, :ecdf) && data in (:samples, :both))
    needs_histogram = mode in (:pdf, :qq) ||
                      (mode in (:hist, :ecdf) && data in (:pdf, :both))
    return needs_samples, needs_histogram
end

function PlotBuilder.resolve(
        ::Type{MCDistributionPlotDefinition},
        recipe::PlotBuilder.PlotRecipe
)
    input = recipe.input
    input.selector isa Function || throw(ArgumentError("selector must be a function"))
    input.mode in (:hist, :pdf, :ecdf, :qq) || throw(
        ArgumentError("mode must be :hist, :pdf, :ecdf, or :qq"),
    )
    input.data in (:samples, :pdf, :both) || throw(
        ArgumentError("data must be :samples, :pdf, or :both"),
    )
    input.nbins === nothing || input.nbins isa Int ||
        throw(
            ArgumentError("nbins must be an integer or nothing"),
        )
    input.nbins === nothing || input.nbins > 0 ||
        throw(
            ArgumentError("nbins must be positive"),
        )
    input.normalization in (:none, :pdf, :density, :probability) || throw(
        ArgumentError("normalization must be :none, :pdf, :density, or :probability"),
    )
    recipe.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    return PlotBuilder.PlotRecipe(
        MCDistributionPlotDefinition,
        recipe.object,
        input,
        recipe.renderer
    )
end

function PlotBuilder.fetch(
        ::Type{MCDistributionPlotDefinition},
        recipe::PlotBuilder.PlotRecipe
)
    length(recipe.object) == 1 || throw(ArgumentError(
        "Monte Carlo distribution plots require one outer Gridspace point",
    ))
    input = recipe.input
    representation = only(UQ.result(recipe.object))
    statistic_product = only(statistics(recipe.object))
    selection = _mc_selection(
        representation,
        statistic_product,
        input.selector,
        input.ijk
    )
    tag, target = _mc_target_unit(
        statistic_product,
        input.selector,
        input.length_unit,
        input.quantity_units
    )
    sample_products = samples(recipe.object)
    histogram_products = histograms(recipe.object)
    sample_payload = if sample_products === nothing
        nothing
    else
        request = _mc_request(input.selector, selection, true)
        observables(
            only(sample_products),
            (selected = request,);
            units = (selected = target,)
        ).selected
    end
    histogram_payload = if histogram_products === nothing
        nothing
    else
        request = _mc_request(input.selector, selection, false)
        observables(
            only(histogram_products),
            (selected = request,);
            units = (selected = target,)
        ).selected
    end
    values = sample_payload === nothing ? nothing : sample_payload.values
    histogram = histogram_payload === nothing ? nothing : histogram_payload.values
    needs_samples, needs_histogram = _mc_required_data(input.mode, input.data)
    needs_samples && values === nothing &&
        throw(
            ArgumentError("this plot requires retained Monte Carlo samples"),
        )
    if needs_histogram && histogram === nothing
        values === nothing && throw(
            ArgumentError("this plot requires retained samples or histograms"),
        )
        histogram = _histogram(values, input.nbins)
    end

    sample_histogram = values === nothing ? nothing : _histogram(values, input.nbins)
    bins = if input.nbins !== nothing && sample_histogram !== nothing
        sample_histogram.edges
    elseif histogram !== nothing
        histogram.edges
    elseif sample_histogram !== nothing
        sample_histogram.edges
    else
        Float64[]
    end
    effective_normalization = input.data in (:pdf, :both) ? :pdf :
                              input.normalization
    quantity_payload = sample_payload === nothing ?
                       (; values = nothing, quantity = tag, unit = target) : sample_payload
    resolved = (;
        values,
        histogram,
        bins,
        effective_normalization,
        quantity_payload,
        selection
    )
    return PlotBuilder.PlotRecipe(
        MCDistributionPlotDefinition,
        recipe.object,
        merge(input, resolved),
        recipe.renderer
    )
end

function PlotBuilder._recipe_variant(
        ::Type{MCDistributionPlotDefinition},
        recipe::PlotBuilder.PlotRecipe
)
    return Val(recipe.input.mode)
end

function PlotBuilder._composition(
        ::Type{MCDistributionPlotDefinition},
        ::Val,
        ::PlotBuilder.PlotRecipe
)
    return Val(:overlay)
end

_mc_data_facets(::Val{:samples}, sample, density) = (sample,)
_mc_data_facets(::Val{:pdf}, sample, density) = (density,)
_mc_data_facets(::Val{:both}, sample, density) = (sample, density)

function PlotBuilder._series_items(
        ::Type{MCDistributionPlotDefinition},
        ::Val{:hist},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return _mc_data_facets(
        Val(recipe.input.data),
        Val(:samples),
        Val(:histogram_pdf)
    )
end

function PlotBuilder._series_items(
        ::Type{MCDistributionPlotDefinition},
        ::Val{:pdf},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return (Val(:histogram_pdf),)
end

function PlotBuilder._series_items(
        ::Type{MCDistributionPlotDefinition},
        ::Val{:ecdf},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return _mc_data_facets(
        Val(recipe.input.data),
        Val(:empirical_cdf),
        Val(:histogram_cdf)
    )
end

function PlotBuilder._series_items(
        ::Type{MCDistributionPlotDefinition},
        ::Val{:qq},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return (Val(:quantiles), Val(:reference))
end

function _mc_values(recipe::PlotBuilder.PlotRecipe)
    recipe.input.values === nothing && throw(
        ArgumentError("Monte Carlo samples were not retained"),
    )
    return recipe.input.values
end

function _mc_histogram_model(recipe::PlotBuilder.PlotRecipe)
    recipe.input.histogram === nothing && throw(
        ArgumentError("Monte Carlo histograms were not retained or derived"),
    )
    return recipe.input.histogram
end

function _mc_cdf_grid(recipe::PlotBuilder.PlotRecipe)
    limits = if recipe.input.histogram !== nothing
        extrema(recipe.input.histogram.edges)
    else
        extrema(_mc_values(recipe))
    end
    lower, upper = limits
    padding = iszero(upper - lower) ? max(abs(lower), 1.0) * 0.05 :
              0.05 * (upper - lower)
    return collect(range(lower - padding, upper + padding; length = 500))
end

function _mc_qq_values(recipe::PlotBuilder.PlotRecipe)
    values = sort(_mc_values(recipe))
    histogram = _mc_histogram_model(recipe)
    probabilities = ((1:length(values)) .- 0.5) ./ length(values)
    histogram_values = _mc_histogram_quantile.(Ref(histogram), probabilities)
    return values, histogram_values
end

function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:samples}
)
    :histogram
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:histogram_pdf}
)
    :stairs
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:histogram_cdf}
)
    :line
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:empirical_cdf}
)
    :line
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:quantiles}
)
    :scatter
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:reference}
)
    :line
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:samples}
)
    return _mc_values(recipe)
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:histogram_pdf}
)
    return _mc_histogram_model(recipe).edges
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:histogram_pdf}
)
    density = _mc_histogram_model(recipe).density
    return [density; last(density)]
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Union{Val{:histogram_cdf}, Val{:empirical_cdf}}
)
    return _mc_cdf_grid(recipe)
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:histogram_cdf}
)
    grid = _mc_cdf_grid(recipe)
    histogram = _mc_histogram_model(recipe)
    return _mc_histogram_cdf.(Ref(histogram), grid)
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:empirical_cdf}
)
    grid = _mc_cdf_grid(recipe)
    values = sort(_mc_values(recipe))
    return _mc_empirical_cdf.(Ref(values), grid)
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:quantiles}
)
    values, _ = _mc_qq_values(recipe)
    return values
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Val{:quantiles}
)
    _, values = _mc_qq_values(recipe)
    return values
end

function PlotBuilder.series_values(
        ::Type{MCDistributionPlotDefinition}, ::Val,
        ::Union{Val{:x}, Val{:y}}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:reference}
)
    sample_values, histogram_values = _mc_qq_values(recipe)
    return collect(extrema(vcat(sample_values, histogram_values)))
end

function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:samples}
)
    "samples"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:histogram_pdf}
)
    "model PDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:histogram_cdf}
)
    "model CDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:empirical_cdf}
)
    "empirical CDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:quantiles}
)
    "quantiles"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:reference}
)
    "perfect fit"
end

function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotDefinition}, ::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:samples}
)
    return (;
        bins = recipe.input.bins,
        normalization = recipe.input.effective_normalization
    )
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:histogram_pdf}
)
    (; step = :post, color = :red, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:histogram_cdf}
)
    (; color = :red, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:empirical_cdf}
)
    (; color = :blue, linestyle = :dash, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:quantiles}
)
    (; color = :steelblue, markersize = 6)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::Val{:reference}
)
    (; color = :black, linestyle = :dash, linewidth = 2)
end

function _mc_title(recipe::PlotBuilder.PlotRecipe, suffix::AbstractString)
    symbol = Units.symbol(recipe.input.quantity_payload.quantity)
    selection = recipe.input.selection
    indices = selection === nothing ? "" : "[$(join(selection, ','))]"
    return "$symbol$indices $suffix"
end

function PlotBuilder.default_title(
        ::Type{MCDistributionPlotDefinition}, ::Val{:hist}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "histogram")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotDefinition}, ::Val{:pdf}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "probability density")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotDefinition}, ::Val{:ecdf}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "cumulative distribution")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotDefinition}, ::Val{:qq}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "Q-Q plot")
end

function PlotBuilder.axis_payload(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return recipe.input.quantity_payload
end
function PlotBuilder.axis_payload(
        ::Type{MCDistributionPlotDefinition}, ::Val{:qq}, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return recipe.input.quantity_payload
end

function _statistical_payload(tag::Symbol, unit::Units.UnitExpr = Units.UnitExpr())
    return (; values = nothing, quantity = Units.QuantityTag{tag}(), unit)
end

function PlotBuilder.axis_payload(
        ::Type{MCDistributionPlotDefinition}, ::Val{:hist}, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    normalization = recipe.input.effective_normalization
    normalization === :none && return _statistical_payload(:sample_count)
    normalization === :probability && return _statistical_payload(:probability)
    return _statistical_payload(
        :probability_density,
        inv(recipe.input.quantity_payload.unit)
    )
end

function PlotBuilder.axis_payload(
        ::Type{MCDistributionPlotDefinition}, ::Val{:pdf}, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return _statistical_payload(
        :probability_density,
        inv(recipe.input.quantity_payload.unit)
    )
end

function PlotBuilder.axis_payload(
        ::Type{MCDistributionPlotDefinition}, ::Val{:ecdf}, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return _statistical_payload(:cumulative_probability)
end

function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotDefinition}, ::Val{:qq}, ::Val{:x},
        quantity::Units.QuantityTag, unit::Units.UnitExpr,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return "Sample quantiles [$(Units.label(unit))]"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotDefinition}, ::Val{:qq}, ::Val{:y},
        quantity::Units.QuantityTag, unit::Units.UnitExpr,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return "Model quantiles [$(Units.label(unit))]"
end

function PlotBuilder.view_key(
        ::Type{MCDistributionPlotDefinition}, ::Val,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return (;
        selector = recipe.input.selector,
        selection = recipe.input.selection,
        mode = recipe.input.mode
    )
end

function PlotBuilder.default_figsize(
        ::Type{MCDistributionPlotDefinition}, ::Val,
        recipe::PlotBuilder.PlotRecipe, page_key
)
    return recipe.renderer.fig_size
end

function PlotBuilder.axis_exponent(
        ::Type{MCDistributionPlotDefinition}, ::Val, ::Val{dimension},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        series::Vector{PlotBuilder.SeriesSpec}
) where {dimension}
    return _mc_plot_exponent(series, dimension)
end

function PlotBuilder.page_identity(
        ::Type{MCDistributionPlotDefinition}, ::Val,
        recipe::PlotBuilder.PlotRecipe, page_key
)
    return (;
        selector = recipe.input.selector,
        selection = recipe.input.selection,
        mode = recipe.input.mode,
        data = recipe.input.data
    )
end
