"""
$(TYPEDEF)

Prepare marginal Monte Carlo distributions for a loaded plotting extension.
"""
struct MCDistributionPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
$(TYPEDEF)

Store one detached Monte Carlo distribution page for the standard plotting
shell.

$(TYPEDFIELDS)
"""
struct MCDistributionPagePayload{X, Y, L, S}
    "Published abscissa observation."
    x_observation::X
    "Published ordinate observation."
    y_observation::Y
    "Optional abscissa-label override."
    xlabel::Union{Nothing, String}
    "Optional ordinate-label override."
    ylabel::Union{Nothing, String}
    "Detached distribution layers in draw order."
    layers::L
    "Captured runtime state used for current-state SVG replay."
    runtime::S
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
    return searchsortedlast(sorted_values, x) / length(sorted_values)
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
function PlotBuilder.input_defaults(::Type{MCDistributionPlotDefinition}, ::MonteCarloResult)
    return _mc_input_defaults()
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
        ::MonteCarloResult,
        request::NamedTuple
)
    input = request.input
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
    request.renderer.fig_size isa Tuple{Int, Int} || throw(
        ArgumentError("fig_size must be a tuple of two integers"),
    )
    return request
end

function _mc_values(values)
    values === nothing && throw(
        ArgumentError("Monte Carlo samples were not retained"),
    )
    return values
end

function _mc_histogram_model(histogram)
    histogram === nothing && throw(
        ArgumentError("Monte Carlo histograms were not retained or derived"),
    )
    return histogram
end

function _mc_cdf_grid(values, histogram)
    limits = histogram === nothing ? extrema(_mc_values(values)) : extrema(histogram.edges)
    lower, upper = limits
    padding = iszero(upper - lower) ? max(abs(lower), 1.0) * 0.05 :
              0.05 * (upper - lower)
    return collect(range(lower - padding, upper + padding; length = 500))
end

function _mc_qq_values(values, histogram)
    sorted = sort(_mc_values(values))
    model = _mc_histogram_model(histogram)
    probabilities = ((1:length(sorted)) .- 0.5) ./ length(sorted)
    model_values = _mc_histogram_quantile.(Ref(model), probabilities)
    return sorted, model_values
end

function _mc_histogram_layer(state)
    return (;
        kind = Val(:histogram),
        x = _mc_values(state.values),
        y = nothing,
        label = "samples",
        group = :samples,
        style = (;
            bins = state.bins,
            normalization = state.effective_normalization
        )
    )
end

function _mc_pdf_layer(state)
    histogram = _mc_histogram_model(state.histogram)
    return (;
        kind = Val(:stairs),
        x = histogram.edges,
        y = [histogram.density; last(histogram.density)],
        label = "model PDF",
        group = :histogram_pdf,
        style = (; step = :post, color = :red, linewidth = 2)
    )
end

function _mc_model_cdf_layer(state)
    grid = _mc_cdf_grid(state.values, state.histogram)
    histogram = _mc_histogram_model(state.histogram)
    return (;
        kind = Val(:line),
        x = grid,
        y = _mc_histogram_cdf.(Ref(histogram), grid),
        label = "model CDF",
        group = :histogram_cdf,
        style = (; color = :red, linewidth = 2)
    )
end

function _mc_empirical_cdf_layer(state)
    grid = _mc_cdf_grid(state.values, state.histogram)
    sorted = sort(_mc_values(state.values))
    return (;
        kind = Val(:line),
        x = grid,
        y = _mc_empirical_cdf.(Ref(sorted), grid),
        label = "empirical CDF",
        group = :empirical_cdf,
        style = (; color = :blue, linestyle = :dash, linewidth = 2)
    )
end

function _mc_qq_layers(state)
    sample_values, model_values = _mc_qq_values(state.values, state.histogram)
    reference = collect(extrema(vcat(sample_values, model_values)))
    return (
        (;
            kind = Val(:scatter),
            x = sample_values,
            y = model_values,
            label = "quantiles",
            group = :quantiles,
            style = (; color = :steelblue, markersize = 6)
        ),
        (;
            kind = Val(:line),
            x = reference,
            y = reference,
            label = "perfect fit",
            group = :reference,
            style = (; color = :black, linestyle = :dash, linewidth = 2)
        )
    )
end

_mc_layers(::Val{:hist}, ::Val{:samples}, state) = (_mc_histogram_layer(state),)
_mc_layers(::Val{:hist}, ::Val{:pdf}, state) = (_mc_pdf_layer(state),)
function _mc_layers(::Val{:hist}, ::Val{:both}, state)
    return (_mc_histogram_layer(state), _mc_pdf_layer(state))
end
_mc_layers(::Val{:pdf}, ::Val, state) = (_mc_pdf_layer(state),)
_mc_layers(::Val{:ecdf}, ::Val{:samples}, state) = (_mc_empirical_cdf_layer(state),)
_mc_layers(::Val{:ecdf}, ::Val{:pdf}, state) = (_mc_model_cdf_layer(state),)
function _mc_layers(::Val{:ecdf}, ::Val{:both}, state)
    return (_mc_empirical_cdf_layer(state), _mc_model_cdf_layer(state))
end
_mc_layers(::Val{:qq}, ::Val, state) = _mc_qq_layers(state)

function _mc_axis_values(layers, field::Symbol)
    values = Float64[]
    for layer in layers
        data = getproperty(layer, field)
        data === nothing && continue
        append!(values, Float64.(data))
    end
    return values
end

function _statistical_observation(
        values,
        tag::Symbol,
        unit::Units.UnitExpr = Units.UnitExpr()
)
    return (; values, quantity = Units.Quantity{tag}(), unit)
end

function _mc_y_observation(::Val{:hist}, state, quantity_observation)
    normalization = state.effective_normalization
    normalization === :none && return _statistical_observation(
        _mc_axis_values(state.layers, :y),
        :sample_count
    )
    normalization === :probability && return _statistical_observation(
        _mc_axis_values(state.layers, :y),
        :probability
    )
    return _statistical_observation(
        _mc_axis_values(state.layers, :y),
        :probability_density,
        inv(quantity_observation.unit)
    )
end

function _mc_y_observation(::Val{:pdf}, state, quantity_observation)
    return _statistical_observation(
        _mc_axis_values(state.layers, :y),
        :probability_density,
        inv(quantity_observation.unit)
    )
end

function _mc_y_observation(::Val{:ecdf}, state, quantity_observation)
    return _statistical_observation(
        _mc_axis_values(state.layers, :y),
        :cumulative_probability
    )
end

function _mc_y_observation(::Val{:qq}, state, quantity_observation)
    return (;
        values = _mc_axis_values(state.layers, :y),
        quantity = quantity_observation.quantity,
        unit = quantity_observation.unit
    )
end

_mc_title_suffix(::Val{:hist}) = "histogram"
_mc_title_suffix(::Val{:pdf}) = "probability density"
_mc_title_suffix(::Val{:ecdf}) = "cumulative distribution"
_mc_title_suffix(::Val{:qq}) = "Q-Q plot"

function _mc_title(quantity, selection, mode::Val)
    scientific_symbol = Units.symbol(quantity)
    indices = selection === nothing ? "" : "[$(join(selection, ','))]"
    return "$scientific_symbol$indices $(_mc_title_suffix(mode))"
end

function _mc_axis_labels(::Val{:qq}, quantity_observation)
    displayed_unit = Units.label(quantity_observation.unit)
    return (
        "Sample quantiles [$displayed_unit]",
        "Model quantiles [$displayed_unit]"
    )
end
_mc_axis_labels(::Val, quantity_observation) = (nothing, nothing)

function PlotBuilder.fetch(
        ::Type{MCDistributionPlotDefinition},
        result::MonteCarloResult,
        request::NamedTuple
)
    length(result) == 1 || throw(ArgumentError(
        "Monte Carlo distribution plots require one outer Gridspace point",
    ))
    input = request.input
    representation = only(result)
    statistic_product = only(statistics(result))
    selection = _mc_selection(
        representation,
        statistic_product,
        input.selector,
        input.ijk
    )
    quantity = Units.quantity(input.selector)
    target = only(unit_targets(
        (input.selector,),
        basis(statistic_product);
        length_prefix = input.length_unit,
        overrides = input.quantity_units
    ))
    sample_products = samples(result)
    histogram_products = histograms(result)
    sample_payload = if sample_products === nothing
        nothing
    else
        observable_request = _mc_request(input.selector, selection, true)
        observables(
            only(sample_products),
            (selected = observable_request,);
            units = (selected = target,)
        ).selected
    end
    histogram_payload = if histogram_products === nothing
        nothing
    else
        observable_request = _mc_request(input.selector, selection, false)
        observables(
            only(histogram_products),
            (selected = observable_request,);
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
    effective_normalization = input.data in (:pdf, :both) ? :pdf : input.normalization
    layer_state = (; values, histogram, bins, effective_normalization)
    mode = Val(input.mode)
    layers = _mc_layers(mode, Val(input.data), layer_state)
    quantity_observation = (;
        values = _mc_axis_values(layers, :x),
        quantity,
        unit = target
    )
    state = merge(layer_state, (; layers))
    y_observation = _mc_y_observation(mode, state, quantity_observation)
    xlabel, ylabel = _mc_axis_labels(mode, quantity_observation)
    title = _mc_title(quantity, selection, mode)
    key = (;
        selector = input.selector,
        selection,
        mode = input.mode,
        data = input.data
    )
    page = MCDistributionPagePayload(
        quantity_observation,
        y_observation,
        xlabel,
        ylabel,
        layers,
        nothing
    )
    return PlotBuilder.PlotPage[
        PlotBuilder.PlotPage(
        title,
        request.renderer.fig_size,
        merge((; page = 1), key),
        page;
        legend = PlotBuilder.LegendDefinition(),
        export_definition = PlotBuilder.ExportDefinition(
            theme = request.renderer.export_theme,
            name = title,
            open_file = request.renderer.open_export
        )
    ),
    ]
end
