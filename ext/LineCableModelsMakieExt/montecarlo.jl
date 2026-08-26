const _MonteCarloScientificSelector = Union{
    typeof(LineCableModels.R),
    typeof(LineCableModels.L),
    typeof(LineCableModels.C),
    typeof(LineCableModels.G)
}

function _monte_carlo_marginal(result::LineCableModels.MonteCarloResult, request)
    length(result) == 1 || throw(ArgumentError(
        "Monte Carlo plots require exactly one outer Gridspace point",
    ))
    represented = only(LineCableModels.result(result))
    if represented isa LineCableModels.CableConstants
        request isa _MonteCarloScientificSelector || throw(ArgumentError(
            "cable-constant Monte Carlo plots require a bare R, L, or C selector",
        ))
        request in (LineCableModels.R, LineCableModels.L, LineCableModels.C) ||
            throw(ArgumentError("cable-constant selectors are R, L, and C"))
        return (; scientific_selector = request, indices = ())
    end
    request isa Tuple && length(request) == 4 || throw(ArgumentError(
        "line-parameter Monte Carlo plots require `@observe quantity[i, j, k]`",
    ))
    scientific_selector, i, j, k = request
    scientific_selector isa _MonteCarloScientificSelector || throw(ArgumentError(
        "line-parameter selectors are R, L, C, and G",
    ))
    all(index -> index isa Integer, (i, j, k)) || throw(ArgumentError(
        "Monte Carlo plot matrix and frequency indices must be integers",
    ))
    checkbounds(LineCableModels.observe(represented, scientific_selector), i, j, k)
    return (; scientific_selector, indices = (Int(i), Int(j), Int(k)))
end

function _monte_carlo_request(product, marginal, retained_samples::Bool)
    point = 1
    isempty(marginal.indices) && return retained_samples ?
        (product, marginal.scientific_selector, point, Colon()) :
        (product, marginal.scientific_selector, point)
    return retained_samples ?
           (product, marginal.scientific_selector, point, marginal.indices..., Colon()) :
           (product, marginal.scientific_selector, point, marginal.indices...)
end

function _monte_carlo_publication(
        result,
        request;
        need_samples::Bool,
        need_model::Bool,
        bins,
        length_unit::Symbol,
        quantity_units
)
    marginal = _monte_carlo_marginal(result, request)
    retained_model = LineCableModels.histograms(result) !== nothing
    publish_samples = need_samples || (need_model && !retained_model)
    publish_model = need_model && retained_model
    requests = if publish_samples && publish_model
        (
            sample = _monte_carlo_request(
                LineCableModels.samples, marginal, true),
            model = _monte_carlo_request(
                LineCableModels.histograms, marginal, false)
        )
    elseif publish_samples
        (sample = _monte_carlo_request(LineCableModels.samples, marginal, true),)
    elseif publish_model
        (model = _monte_carlo_request(LineCableModels.histograms, marginal, false),)
    else
        error("Monte Carlo plotting requested no published product")
    end
    target = only(LineCableModels.Grammar.unit_targets(
        (marginal.scientific_selector,),
        LineCableModels.basis(result);
        length_prefix = length_unit,
        overrides = quantity_units
    ))
    units = NamedTuple{keys(requests)}(ntuple(_ -> target, length(requests)))
    published = LineCableModels.observables(result, requests; units)
    sample = haskey(published, :sample) ? published.sample : nothing
    model = if haskey(published, :model)
        published.model
    elseif need_model
        bins === nothing || bins isa Int || throw(ArgumentError(
            "derived HistogramDensity bins must be an integer or nothing",
        ))
        (;
            values = LineCableModels.UQ._histogram(sample.values, bins),
            quantity = sample.quantity,
            unit = sample.unit
        )
    else
        nothing
    end
    return (; marginal, sample, model)
end

function _distribution_observation(values, tag::Symbol, unit = LineCableModels.Units.UnitExpr())
    return (;
        values,
        quantity = LineCableModels.Units.Quantity{tag}(),
        unit
    )
end

function _monte_carlo_title(publication, suffix::AbstractString)
    payload = publication.sample === nothing ? publication.model : publication.sample
    symbol = LineCableModels.Units.symbol(payload.quantity)
    indices = isempty(publication.marginal.indices) ? "" :
              "[$(join(publication.marginal.indices, ','))]"
    return "$symbol$indices $suffix"
end

function _register_distribution!(context, axis, groups, labels, data)
    registration = only(panel for panel in context.panels if panel.axis === axis)
    PlotBuilder.register!(context, axis;
        xmetadata = registration.metadata.xaxis,
        ymetadata = registration.metadata.yaxis,
        groups,
        labels,
        data)
    return axis
end

function _distribution_window(
        callback;
        title,
        fig_size,
        backend,
        display_plot,
        controls,
        export_theme,
        open_export
)
    return PlotBuilder.plotwindow(callback;
        title,
        size = fig_size,
        backend,
        display_plot,
        controls,
        export_theme,
        open_export)
end

function Makie.hist(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        bins = nothing,
        normalization = :none,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        title = nothing,
        fig_size = (800, 400),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        kwargs...
)
    publication = _monte_carlo_publication(result, request;
        need_samples = true,
        need_model = false,
        bins,
        length_unit,
        quantity_units)
    heading = title === nothing ? _monte_carlo_title(publication, "histogram") :
              String(title)
    sample = publication.sample
    y_observation = normalization === :none ?
                    _distribution_observation(Float64[], :sample_count) :
                    normalization === :probability ?
                    _distribution_observation(Float64[], :probability) :
                    _distribution_observation(
        Float64[], :probability_density, inv(sample.unit))
    return _distribution_window(;
        title = heading,
        fig_size,
        backend,
        display_plot,
        controls,
        export_theme,
        open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], sample, y_observation;
            title = heading)
        plot = bins === nothing ?
               hist!(axis, sample.values; normalization, label = "samples", kwargs...) :
               hist!(axis, sample.values; bins, normalization, label = "samples", kwargs...)
        _register_distribution!(context, axis,
            (samples = (plot,),),
            (samples = "samples",),
            ((; xdata = sample.values, ydata = nothing,
                group = :samples, label = "samples"),))
    end
end

function Makie.stairs(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        bins = nothing,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        title = nothing,
        fig_size = (800, 400),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        kwargs...
)
    publication = _monte_carlo_publication(result, request;
        need_samples = false,
        need_model = true,
        bins,
        length_unit,
        quantity_units)
    heading = title === nothing ? _monte_carlo_title(publication, "probability density") :
              String(title)
    model = publication.model
    abscissa = (;
        values = model.values.edges,
        quantity = model.quantity,
        unit = model.unit
    )
    ordinate = _distribution_observation(
        model.values.density, :probability_density, inv(model.unit))
    return _distribution_window(;
        title = heading, fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], abscissa, ordinate;
            title = heading)
        y = [model.values.density; last(model.values.density)]
        plot = stairs!(axis, model.values.edges, y;
            step = :post, label = "model PDF", kwargs...)
        _register_distribution!(context, axis,
            (model = (plot,),),
            (model = "model PDF",),
            ((; xdata = model.values.edges, ydata = y,
                group = :model, label = "model PDF"),))
    end
end

function Makie.ecdfplot(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        title = nothing,
        fig_size = (800, 400),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        kwargs...
)
    publication = _monte_carlo_publication(result, request;
        need_samples = true,
        need_model = false,
        bins = nothing,
        length_unit,
        quantity_units)
    heading = title === nothing ?
              _monte_carlo_title(publication, "empirical cumulative distribution") :
              String(title)
    sample = publication.sample
    ordinate = _distribution_observation(Float64[], :cumulative_probability)
    return _distribution_window(;
        title = heading, fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], sample, ordinate;
            title = heading)
        plot = ecdfplot!(axis, sample.values; label = "empirical CDF", kwargs...)
        _register_distribution!(context, axis,
            (empirical = (plot,),),
            (empirical = "empirical CDF",),
            ((; xdata = sample.values, ydata = nothing,
                group = :empirical, label = "empirical CDF"),))
    end
end

function _model_cdf_grid(histogram)
    lower, upper = extrema(histogram.edges)
    padding = iszero(upper - lower) ? max(abs(lower), 1.0) * 0.05 :
              0.05 * (upper - lower)
    return collect(range(lower - padding, upper + padding; length = 500))
end

function Makie.lines(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        bins = nothing,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        title = nothing,
        fig_size = (800, 400),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        kwargs...
)
    publication = _monte_carlo_publication(result, request;
        need_samples = false,
        need_model = true,
        bins,
        length_unit,
        quantity_units)
    heading = title === nothing ?
              _monte_carlo_title(publication, "model cumulative distribution") :
              String(title)
    model = publication.model
    grid = _model_cdf_grid(model.values)
    probability = LineCableModels.UQ.cumulative_probability.(Ref(model.values), grid)
    abscissa = (; values = grid, quantity = model.quantity, unit = model.unit)
    ordinate = _distribution_observation(probability, :cumulative_probability)
    return _distribution_window(;
        title = heading, fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], abscissa, ordinate;
            title = heading)
        plot = lines!(axis, grid, probability; label = "model CDF", kwargs...)
        _register_distribution!(context, axis,
            (model = (plot,),),
            (model = "model CDF",),
            ((; xdata = grid, ydata = probability,
                group = :model, label = "model CDF"),))
    end
end

function Makie.qqplot(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        qqline::Symbol = :identity,
        bins = nothing,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        title = nothing,
        fig_size = (800, 400),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        kwargs...
)
    qqline in (:identity, :none) || throw(ArgumentError(
        "qqline must be :identity or :none",
    ))
    publication = _monte_carlo_publication(result, request;
        need_samples = true,
        need_model = true,
        bins,
        length_unit,
        quantity_units)
    heading = title === nothing ? _monte_carlo_title(publication, "Q-Q plot") :
              String(title)
    pairs = LineCableModels.UQ.quantile_pairs(
        publication.model.values,
        publication.sample.values
    )
    sample = (;
        values = pairs.sample,
        quantity = publication.sample.quantity,
        unit = publication.sample.unit
    )
    model = (;
        values = pairs.model,
        quantity = publication.model.quantity,
        unit = publication.model.unit
    )
    displayed_unit = LineCableModels.Units.label(sample.unit)
    return _distribution_window(;
        title = heading, fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], sample, model;
            title = heading,
            xlabel = "Sample quantiles [$displayed_unit]",
            ylabel = "Model quantiles [$displayed_unit]")
        points = scatter!(axis, pairs.sample, pairs.model;
            label = "quantiles", kwargs...)
        if qqline === :identity
            reference = collect(pairs.reference)
            line = lines!(axis, reference, reference;
                color = :black, linestyle = :dash,
                linewidth = 2, label = "perfect fit")
            _register_distribution!(context, axis,
                (quantiles = (points,), reference = (line,)),
                (quantiles = "quantiles", reference = "perfect fit"),
                (
                    (; xdata = pairs.sample, ydata = pairs.model,
                        group = :quantiles, label = "quantiles"),
                    (; xdata = reference, ydata = reference,
                        group = :reference, label = "perfect fit")
                ))
        else
            _register_distribution!(context, axis,
                (quantiles = (points,),),
                (quantiles = "quantiles",),
                ((; xdata = pairs.sample, ydata = pairs.model,
                    group = :quantiles, label = "quantiles"),))
        end
    end
end
