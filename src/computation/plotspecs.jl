"""
$(TYPEDEF)

PlotBuilder recipe for marginal Monte Carlo histograms, probability densities,
cumulative distributions, and Q-Q plots.
"""
struct MCDistributionPlotSpec <: PlotBuilder.AbstractPlotSpec end

struct MCSeriesKey{K} end

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

function _mc_selection(
        result::MonteCarloResult,
        ::DataModel.CableConstants,
        quantity::Symbol,
        ijk
)
    quantity in (:R, :L, :C) || throw(
        ArgumentError("cable-constant quantities are :R, :L, and :C"),
    )
    ijk === nothing || throw(
        ArgumentError("cable-constant Monte Carlo results do not use matrix indices"),
    )
    sample_values = result.samples === nothing ? nothing :
                    getproperty(result.samples, quantity)
    histogram_value = result.histograms === nothing ? nothing :
                      getproperty(result.histograms, quantity)
    return sample_values, histogram_value, nothing
end

function _mc_selection(
        result::MonteCarloResult,
        ::Engine.LineParameters,
        quantity::Symbol,
        ijk
)
    quantity in (:R, :L, :C, :G) || throw(
        ArgumentError("line-parameter quantities are :R, :L, :C, and :G"),
    )
    selection = ijk === nothing ? (1, 1, 1) : ijk
    selection isa NTuple{3, Int} || throw(
        ArgumentError("ijk must be a tuple (i, j, k)"),
    )
    observable = getproperty(result.statistics, quantity)
    checkbounds(observable, selection...)
    sample_values = result.samples === nothing ? nothing :
                    collect(
        view(getproperty(result.samples, quantity), selection..., :),
    )
    histogram_value = result.histograms === nothing ? nothing :
                      getproperty(result.histograms, quantity)[selection...]
    return sample_values, histogram_value, selection
end

function _mc_selection(result::MonteCarloResult, quantity::Symbol, ijk)
    return _mc_selection(result, result.representation, quantity, ijk)
end

function _mc_target_unit(result::MonteCarloResult, quantity, length_unit, quantity_units)
    result_basis = result.representation isa DataModel.CableConstants ?
                   :per_length : basis(result.representation)
    resolved = _mc_unit(
        quantity,
        result_basis,
        length_unit,
        quantity_units
    )
    return resolved.quantity, resolved.units, resolved.scale
end

function _scaled_histogram(histogram::HistogramPDF, conversion)
    conversion > zero(conversion) || throw(
        ArgumentError("unit conversion must be positive"),
    )
    return HistogramPDF(
        histogram.edges .* conversion,
        histogram.density ./ conversion
    )
end

function _mc_histogram_cdf(histogram::HistogramPDF, x::Real)
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

function _mc_histogram_quantile(histogram::HistogramPDF, probability::Real)
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
        quantity = :R,
        ijk = nothing,
        mode = :hist,
        data = :samples,
        length_unit = :kilo,
        quantity_units = nothing,
        nbins = nothing,
        normalization = :none
    )
end

PlotBuilder.dispatch_on(::Type{MCDistributionPlotSpec}) = MonteCarloResult
function PlotBuilder.input_kwargs(::Type{MCDistributionPlotSpec})
    return (
        :quantity,
        :ijk,
        :mode,
        :data,
        :length_unit,
        :quantity_units,
        :nbins,
        :normalization
    )
end
PlotBuilder.renderer_kwargs(::Type{MCDistributionPlotSpec}) = (:fig_size,)
function PlotBuilder.input_defaults(::Type{MCDistributionPlotSpec}, ::MonteCarloResult)
    _mc_input_defaults()
end
function PlotBuilder.renderer_defaults(
        ::Type{MCDistributionPlotSpec},
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

function PlotBuilder.resolve_input(
        ::Type{MCDistributionPlotSpec},
        recipe::PlotBuilder.PlotRecipe
)
    input = recipe.input
    input.quantity isa Symbol || throw(ArgumentError("quantity must be a Symbol"))
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

    sample_values, histogram_value, selection = _mc_selection(
        recipe.object,
        input.quantity,
        input.ijk
    )
    tag, target, conversion = _mc_target_unit(
        recipe.object,
        input.quantity,
        input.length_unit,
        input.quantity_units
    )
    values = sample_values === nothing ? nothing :
             collect(sample_values) .* conversion
    histogram = histogram_value === nothing ? nothing :
                _scaled_histogram(histogram_value, conversion)

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
    resolved = (;
        values,
        histogram,
        bins,
        effective_normalization,
        tag,
        target,
        selection
    )
    return PlotBuilder.PlotRecipe(
        recipe.object,
        merge(input, resolved),
        recipe.renderer
    )
end

function PlotBuilder.recipe_mode(
        ::Type{MCDistributionPlotSpec},
        recipe::PlotBuilder.PlotRecipe
)
    return Val(recipe.input.mode)
end

function PlotBuilder.grouping_mode(
        ::Type{MCDistributionPlotSpec},
        ::Val,
        ::PlotBuilder.PlotRecipe
)
    return Val(:overlay)
end

_mc_data_facets(::Val{:samples}, sample, density) = (sample,)
_mc_data_facets(::Val{:pdf}, sample, density) = (density,)
_mc_data_facets(::Val{:both}, sample, density) = (sample, density)

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:hist},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return _mc_data_facets(
        Val(recipe.input.data),
        MCSeriesKey{:samples}(),
        MCSeriesKey{:histogram_pdf}()
    )
end

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:pdf},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return (MCSeriesKey{:histogram_pdf}(),)
end

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:ecdf},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return _mc_data_facets(
        Val(recipe.input.data),
        MCSeriesKey{:empirical_cdf}(),
        MCSeriesKey{:histogram_cdf}()
    )
end

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:qq},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return (MCSeriesKey{:quantiles}(), MCSeriesKey{:reference}())
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
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples}
)
    :histogram
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:histogram_pdf}
)
    :stairs
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:histogram_cdf}
)
    :line
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf}
)
    :line
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles}
)
    :scatter
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference}
)
    :line
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:samples}
)
    return _mc_values(recipe)
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:histogram_pdf}
)
    return _mc_histogram_model(recipe).edges
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:histogram_pdf}
)
    density = _mc_histogram_model(recipe).density
    return [density; last(density)]
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::Union{MCSeriesKey{:histogram_cdf}, MCSeriesKey{:empirical_cdf}}
)
    return _mc_cdf_grid(recipe)
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:histogram_cdf}
)
    grid = _mc_cdf_grid(recipe)
    histogram = _mc_histogram_model(recipe)
    return _mc_histogram_cdf.(Ref(histogram), grid)
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:empirical_cdf}
)
    grid = _mc_cdf_grid(recipe)
    values = sort(_mc_values(recipe))
    return _mc_empirical_cdf.(Ref(values), grid)
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:quantiles}
)
    values, _ = _mc_qq_values(recipe)
    return values
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        ::MCSeriesKey{:quantiles}
)
    _, values = _mc_qq_values(recipe)
    return values
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, ::Val,
        ::Union{Val{:x}, Val{:y}}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference}
)
    sample_values, histogram_values = _mc_qq_values(recipe)
    return collect(extrema(vcat(sample_values, histogram_values)))
end

function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples}
)
    "samples"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:histogram_pdf}
)
    "model PDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:histogram_cdf}
)
    "model CDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf}
)
    "empirical CDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles}
)
    "quantiles"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference}
)
    "perfect fit"
end

function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, ::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples}
)
    return (;
        bins = recipe.input.bins,
        normalization = recipe.input.effective_normalization
    )
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:histogram_pdf}
)
    (; step = :post, color = :red, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:histogram_cdf}
)
    (; color = :red, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf}
)
    (; color = :blue, linestyle = :dash, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles}
)
    (; color = :steelblue, markersize = 6)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, ::Val, ::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference}
)
    (; color = :black, linestyle = :dash, linewidth = 2)
end

function _mc_title(recipe::PlotBuilder.PlotRecipe, suffix::AbstractString)
    symbol = UnitHandler.get_symbol(recipe.input.tag)
    selection = recipe.input.selection
    indices = selection === nothing ? "" : "[$(join(selection, ','))]"
    return "$symbol$indices $suffix"
end

function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:hist}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "histogram")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:pdf}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "probability density")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:ecdf}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "cumulative distribution")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    _mc_title(recipe, "Q-Q plot")
end

function PlotBuilder.axis_quantity(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return recipe.input.tag
end
function PlotBuilder.axis_quantity(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return recipe.input.tag
end
function PlotBuilder.axis_quantity(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:y},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return UnitHandler.QuantityTag{:dimensionless}()
end

function PlotBuilder.axis_unit(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        quantity::UnitHandler.QuantityTag, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    return recipe.input.target
end
function PlotBuilder.axis_unit(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:y},
        quantity::UnitHandler.QuantityTag, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    return recipe.input.target
end
function PlotBuilder.axis_unit(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:y},
        quantity::UnitHandler.QuantityTag, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key
)
    return UnitHandler.Units()
end

function _mc_quantity_label(quantity, unit)
    unit_label = UnitHandler.get_label(unit)
    return isempty(unit_label) ? UnitHandler.get_label(quantity) :
           "$(UnitHandler.get_label(quantity)) [$unit_label]"
end

function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:x},
        quantity::UnitHandler.QuantityTag, unit::UnitHandler.Units,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return "sample quantiles [$(UnitHandler.get_label(unit))]"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{:x},
        quantity::UnitHandler.QuantityTag, unit::UnitHandler.Units,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return _mc_quantity_label(quantity, unit)
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:hist}, ::Val{:y},
        quantity::UnitHandler.QuantityTag, unit::UnitHandler.Units,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    normalization = recipe.input.effective_normalization
    return normalization === :none ? "count" : String(normalization)
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:pdf}, ::Val{:y},
        quantity::UnitHandler.QuantityTag, unit::UnitHandler.Units,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    "probability density"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:ecdf}, ::Val{:y},
        quantity::UnitHandler.QuantityTag, unit::UnitHandler.Units,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    "cumulative probability"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:y},
        quantity::UnitHandler.QuantityTag, unit::UnitHandler.Units,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return "model quantiles [$(UnitHandler.get_label(unit))]"
end

function PlotBuilder.view_key(
        ::Type{MCDistributionPlotSpec}, ::Val,
        recipe::PlotBuilder.PlotRecipe, page_key, view_key
)
    return (;
        quantity = recipe.input.quantity,
        selection = recipe.input.selection,
        mode = recipe.input.mode
    )
end

function PlotBuilder.default_figsize(
        ::Type{MCDistributionPlotSpec}, ::Val,
        recipe::PlotBuilder.PlotRecipe, page_key
)
    return recipe.renderer.fig_size
end

function PlotBuilder.axis_exponent(
        ::Type{MCDistributionPlotSpec}, ::Val, ::Val{dimension},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        series::Vector{PlotBuilder.SeriesSpec}
) where {dimension}
    return _mc_plot_exponent(series, dimension)
end

function PlotBuilder.page_identity(
        ::Type{MCDistributionPlotSpec}, ::Val,
        recipe::PlotBuilder.PlotRecipe, page_key
)
    return (;
        quantity = recipe.input.quantity,
        selection = recipe.input.selection,
        mode = recipe.input.mode,
        data = recipe.input.data
    )
end
