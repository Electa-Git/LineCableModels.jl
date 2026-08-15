struct MCDistributionPlotSpec <: PlotBuilder.AbstractPlotSpec end

function _plot_exponent(series, field::Symbol)
    maximum_value = 0.0
    for item in series
        samples = field === :x ? item.xdata : item.ydata
        samples === nothing && continue
        for sample in samples
            nominal = abs(value(sample))
            nominal isa Real && isfinite(nominal) &&
                (maximum_value = max(maximum_value, Float64(nominal)))
        end
    end
    iszero(maximum_value) && return 0
    exponent = floor(Int, log10(maximum_value))
    return abs(exponent) < 3 ? 0 : exponent
end

function _mc_selection(result::CableConstantsMC, quantity::Symbol, ijk)
    ijk === nothing || throw(ArgumentError("CableConstantsMC does not use matrix indices"))
    sample_values = has_samples(result) ? samples(result, quantity) : nothing
    distribution_value = has_distributions(result) ? distribution(result, quantity) :
                         nothing
    return sample_values, distribution_value, nothing
end

function _mc_selection(result::LineParametersMC, quantity::Symbol, ijk)
    selection = isnothing(ijk) ? (1, 1, 1) : ijk
    selection isa NTuple{3, Int} || throw(ArgumentError("ijk must be a tuple (i, j, k)"))
    i, j, k = selection
    sample_values = has_samples(result) ? samples(result, quantity, i, j, k) : nothing
    distribution_value = has_distributions(result) ?
                         distribution(result, quantity, i, j, k) : nothing
    return sample_values, distribution_value, selection
end

function _mc_target_unit(result, quantity::Symbol, length_unit, quantity_units)
    quantity in (:R, :L, :C, :G) ||
        throw(ArgumentError("quantity must be :R, :L, :C, or :G"))
    result isa CableConstantsMC && quantity === :G &&
        throw(
            ArgumentError("CableConstantsMC does not contain conductance"),
        )
    result_basis = result isa CableConstantsMC ? :per_length : basis(result)
    resolved = UnitHandler.line_component_unit(
        quantity,
        result_basis;
        length_unit,
        quantity_units
    )
    return resolved.quantity, resolved.units, resolved.scale
end

function _scaled_distribution(distribution_value::HistogramPDF, conversion)
    return HistogramPDF(
        distribution_value.edges .* conversion,
        distribution_value.density ./ conversion
    )
end

struct MCSeriesKey{K} end

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

function PlotBuilder.dispatch_on(::Type{MCDistributionPlotSpec})
    Union{CableConstantsMC, LineParametersMC}
end
function PlotBuilder.input_kwargs(::Type{MCDistributionPlotSpec})
    (
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
function PlotBuilder.input_defaults(
        ::Type{MCDistributionPlotSpec},
        ::Union{CableConstantsMC, LineParametersMC}
)
    _mc_input_defaults()
end
function PlotBuilder.renderer_defaults(
        ::Type{MCDistributionPlotSpec},
        ::Union{CableConstantsMC, LineParametersMC}
)
    (; fig_size = (800, 400))
end

function PlotBuilder.resolve_input(::Type{MCDistributionPlotSpec}, recipe::PlotBuilder.PlotRecipe)
    input = recipe.input
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
        throw(ArgumentError("nbins must be positive"))
    recipe.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )

    sample_values, distribution_value, selection = _mc_selection(
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
    values = sample_values === nothing ? nothing : collect(sample_values) .* conversion
    model = distribution_value === nothing ? nothing :
            _scaled_distribution(distribution_value, conversion)
    bin_count = input.nbins === nothing && values !== nothing ? _auto_nbins(values) :
                input.nbins
    bins = if model !== nothing
        model.edges
    elseif values !== nothing
        collect(fit(Histogram, values; nbins = bin_count, closed = :left).edges[1])
    else
        Float64[]
    end
    effective_normalization = input.data in (:pdf, :both) ? :pdf : input.normalization
    resolved = (;
        values,
        model,
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

function PlotBuilder.recipe_mode(::Type{MCDistributionPlotSpec}, recipe::PlotBuilder.PlotRecipe)
    return Val(recipe.input.mode)
end

function PlotBuilder.grouping_mode(
        ::Type{MCDistributionPlotSpec},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:overlay)
end

_mc_hist_facets(::Val{:samples}) = (MCSeriesKey{:samples}(),)
_mc_hist_facets(::Val{:pdf}) = (MCSeriesKey{:model_pdf}(),)
_mc_hist_facets(::Val{:both}) = (MCSeriesKey{:samples}(), MCSeriesKey{:model_pdf}())

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:hist},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return _mc_hist_facets(Val(recipe.input.data))
end

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:pdf},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return (MCSeriesKey{:model_pdf}(),)
end

function PlotBuilder.group_facets(
        ::Type{MCDistributionPlotSpec},
        ::Val{:ecdf},
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return (MCSeriesKey{:model_cdf}(), MCSeriesKey{:empirical_cdf}())
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
    recipe.input.values === nothing && throw(ArgumentError("samples were not retained"))
    return recipe.input.values
end

function _mc_model(recipe::PlotBuilder.PlotRecipe)
    recipe.input.model === nothing &&
        throw(ArgumentError("distributions were not retained"))
    return recipe.input.model
end

function _mc_cdf_grid(recipe::PlotBuilder.PlotRecipe)
    model = _mc_model(recipe)
    lower, upper = extrema(model.edges)
    padding = iszero(upper - lower) ? one(lower) : 0.05 * (upper - lower)
    return collect(range(lower - padding, upper + padding; length = 500))
end

function _mc_qq_values(recipe::PlotBuilder.PlotRecipe)
    values = sort(_mc_values(recipe))
    model = _mc_model(recipe)
    probabilities = ((1:length(values)) .- 0.5) ./ length(values)
    model_values = Distributions.quantile.(Ref(model), probabilities)
    return values, model_values
end

function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples})
    :histogram
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_pdf})
    :stairs
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_cdf})
    :line
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf})
    :line
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles})
    :scatter
end
function PlotBuilder.plot_kind(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference})
    :line
end

function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples})
    return _mc_values(recipe)
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_pdf})
    return _mc_model(recipe).edges
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:y}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_pdf})
    model = _mc_model(recipe)
    return [model.density; last(model.density)]
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_cdf})
    return _mc_cdf_grid(recipe)
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:y}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_cdf})
    x = _mc_cdf_grid(recipe)
    return Distributions.cdf.(Ref(_mc_model(recipe)), x)
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf})
    return _mc_cdf_grid(recipe)
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:y}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf})
    x = _mc_cdf_grid(recipe)
    return StatsBase.ecdf(_mc_values(recipe)).(x)
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles})
    values, _ = _mc_qq_values(recipe)
    return values
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:y}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles})
    _, values = _mc_qq_values(recipe)
    return values
end
function PlotBuilder.series_data(
        ::Type{MCDistributionPlotSpec}, mode::Val, dim::Union{Val{:x}, Val{:y}},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key, ::MCSeriesKey{:reference})
    sample_values, model_values = _mc_qq_values(recipe)
    return collect(extrema(vcat(sample_values, model_values)))
end

function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples})
    "samples"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_pdf})
    "model PDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_cdf})
    "model CDF"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf})
    "empirical"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles})
    "quantiles"
end
function PlotBuilder.legend_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference})
    "perfect fit"
end

function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:samples})
    return (;
        bins = recipe.input.bins, normalization = recipe.input.effective_normalization)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_pdf})
    return (; step = :post, color = :red, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:model_cdf})
    return (; color = :red, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:empirical_cdf})
    return (; color = :blue, linestyle = :dash, linewidth = 2)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:quantiles})
    return (; color = :steelblue, markersize = 6)
end
function PlotBuilder.series_attributes(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key, ::MCSeriesKey{:reference})
    return (; color = :black, linestyle = :dash, linewidth = 2)
end

function _mc_title(recipe::PlotBuilder.PlotRecipe, suffix::AbstractString)
    symbol = UnitHandler.get_symbol(recipe.input.tag)
    selection = recipe.input.selection
    indices = selection === nothing ? "" : "[$(join(selection, ','))]"
    return "$symbol$indices $suffix"
end

function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:hist}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    _mc_title(recipe, "histogram")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:pdf}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    _mc_title(recipe, "probability density")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:ecdf}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    _mc_title(recipe, "cumulative distribution")
end
function PlotBuilder.default_title(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    _mc_title(recipe, "Q-Q plot")
end

function PlotBuilder.axis_quantity(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    return recipe.input.tag
end
function PlotBuilder.axis_quantity(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:y}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    return recipe.input.tag
end
function PlotBuilder.axis_quantity(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:y}, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    return UnitHandler.QuantityTag{:dimensionless}()
end
function PlotBuilder.axis_unit(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x},
        quantity::UnitHandler.QuantityTag, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return recipe.input.target
end
function PlotBuilder.axis_unit(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:y},
        quantity::UnitHandler.QuantityTag, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return recipe.input.target
end
function PlotBuilder.axis_unit(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:y},
        quantity::UnitHandler.QuantityTag, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return UnitHandler.Units()
end

function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:x},
        quantity::UnitHandler.QuantityTag,
        unit::UnitHandler.Units, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return "sample quantiles [$(UnitHandler.get_label(unit))]"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{:x},
        quantity::UnitHandler.QuantityTag,
        unit::UnitHandler.Units, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return "$(UnitHandler.get_label(quantity)) [$(UnitHandler.get_label(unit))]"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:hist}, ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        unit::UnitHandler.Units, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return String(recipe.input.effective_normalization)
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:pdf}, ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        unit::UnitHandler.Units, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return "density"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:ecdf}, ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        unit::UnitHandler.Units, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return "cumulative probability"
end
function PlotBuilder.axis_label(
        ::Type{MCDistributionPlotSpec}, ::Val{:qq}, ::Val{:y},
        quantity::UnitHandler.QuantityTag,
        unit::UnitHandler.Units, recipe::PlotBuilder.PlotRecipe, page_key, view_key)
    return "model quantiles [$(UnitHandler.get_label(unit))]"
end

function PlotBuilder.view_key(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe,
        page_key, view_key)
    return (;
        quantity = recipe.input.quantity,
        selection = recipe.input.selection,
        mode = recipe.input.mode
    )
end

function PlotBuilder.default_figsize(
        ::Type{MCDistributionPlotSpec}, mode::Val, recipe::PlotBuilder.PlotRecipe, page_key)
    return recipe.renderer.fig_size
end

function PlotBuilder.axis_exponent(
        ::Type{MCDistributionPlotSpec}, mode::Val, ::Val{dim},
        recipe::PlotBuilder.PlotRecipe, page_key, view_key,
        series::Vector{PlotBuilder.SeriesSpec}) where {dim}
    return _plot_exponent(series, dim)
end

function PlotBuilder.page_identity(
        ::Type{MCDistributionPlotSpec}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return (;
        quantity = recipe.input.quantity,
        selection = recipe.input.selection,
        mode = recipe.input.mode,
        data = recipe.input.data
    )
end
