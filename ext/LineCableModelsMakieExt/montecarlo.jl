function _monte_carlo_publication(
        result::LineCableModels.MonteCarloResult,
        scientific_selector::Function,
        indices::Tuple;
        need_samples::Bool,
        need_model::Bool,
        bins,
        length_unit::Symbol,
        quantity_units,
        clip::Bool
)
    length(result) == 1 || throw(ArgumentError(
        "Monte Carlo plots require exactly one outer Gridspace point",
    ))
    bins === nothing || bins isa Integer ||
        throw(ArgumentError(
            "Monte Carlo histogram bins must be an integer or nothing",
        ))
    point = 1
    sample_request = isempty(indices) ?
                     (LineCableModels.samples, scientific_selector, point, Colon()) :
                     (
        LineCableModels.samples, scientific_selector, point, indices..., Colon())
    model_request = isempty(indices) ?
                    (LineCableModels.histograms, scientific_selector, point, bins) :
                    (
        LineCableModels.histograms, scientific_selector, point, indices..., bins)
    requests = if need_samples && need_model
        (sample_request, model_request)
    elseif need_samples
        (sample_request,)
    elseif need_model
        (model_request,)
    else
        error("Monte Carlo plotting requested no published product")
    end
    target = only(LineCableModels.Grammar.unit_targets(
        (scientific_selector,),
        LineCableModels.basis(result);
        length_prefix = length_unit,
        overrides = quantity_units
    ))
    units = ntuple(_ -> target, length(requests))
    published = LineCableModels.observables(result, requests; units, clip)
    sample = need_samples ? first(published) : nothing
    model = need_model ? published[need_samples ? 2 : 1] : nothing
    return (; scientific_selector, indices, sample, model)
end

function _monte_carlo_publication(
        result::LineCableModels.MonteCarloResult{<:LineCableModels.CableConstants},
        request;
        kwargs...
)
    if request isa Function
        all(value -> length(value) == 1, result.values) || throw(ArgumentError(
            "a multi-assembly cable-constant plot requires `(selector, assembly)`",
        ))
        publication = _monte_carlo_publication(result, request, (1,); kwargs...)
        return merge(publication, (indices = (),))
    end
    request isa Tuple && length(request) == 2 || throw(ArgumentError(
        "cable-constant Monte Carlo plots require R, L, C, G, or `(selector, assembly)`",
    ))
    scientific_selector, assembly = request
    scientific_selector isa Function && assembly isa Integer || throw(ArgumentError(
        "a cable-constant plot assembly request requires a selector and integer index",
    ))
    return _monte_carlo_publication(
        result, scientific_selector, (Int(assembly),); kwargs...
    )
end

function _monte_carlo_publication(
        result::LineCableModels.MonteCarloResult{<:LineCableModels.LineParameters},
        request;
        kwargs...
)
    request isa Tuple && length(request) == 4 || throw(ArgumentError(
        "line-parameter Monte Carlo plots require `@observe quantity[i, j, k]`",
    ))
    scientific_selector, row, column, frequency = request
    scientific_selector isa Function || throw(ArgumentError(
        "line-parameter Monte Carlo requests must begin with a scientific selector",
    ))
    all(index -> index isa Integer, (row, column, frequency)) || throw(ArgumentError(
        "Monte Carlo plot matrix and frequency indices must be integers",
    ))
    return _monte_carlo_publication(
        result,
        scientific_selector,
        (Int(row), Int(column), Int(frequency));
        kwargs...
    )
end

function _monte_carlo_title(publication, suffix::AbstractString)
    payload = publication.sample === nothing ? publication.model : publication.sample
    symbol = LineCableModels.Units.symbol(payload.quantity)
    indices = isempty(publication.indices) ? "" :
              "[$(join(publication.indices, ','))]"
    return "$symbol$indices $suffix"
end

function Makie.hist(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        bins = nothing,
        normalization = :none,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        clip::Bool = true,
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
        quantity_units,
        clip)
    heading = title === nothing ? _monte_carlo_title(publication, "histogram") :
              String(title)
    sample = publication.sample
    y_observation = normalization === :none ?
                    (;
        values = Float64[],
        quantity = LineCableModels.Units.Quantity{:sample_count}(),
        unit = LineCableModels.Units.UnitExpr()
    ) :
                    normalization === :probability ?
                    (;
        values = Float64[],
        quantity = LineCableModels.Units.Quantity{:probability}(),
        unit = LineCableModels.Units.UnitExpr()
    ) :
                    (;
        values = Float64[],
        quantity = LineCableModels.Units.Quantity{:probability_density}(),
        unit = inv(sample.unit)
    )
    return PlotBuilder.plotwindow(;
        title = heading,
        size = fig_size,
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
        PlotBuilder.register!(context, axis;
            groups = (samples = (plot,),),
            labels = (samples = "samples",),
            data = ((; xdata = sample.values, ydata = nothing,
                group = :samples, label = "samples"),))
    end
end

function Makie.stairs(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        bins = nothing,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        clip::Bool = true,
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
        quantity_units,
        clip)
    heading = title === nothing ? _monte_carlo_title(publication, "probability density") :
              String(title)
    model = publication.model
    abscissa = (;
        values = model.values.edges,
        quantity = model.quantity,
        unit = model.unit
    )
    ordinate = (;
        values = model.values.density,
        quantity = LineCableModels.Units.Quantity{:probability_density}(),
        unit = inv(model.unit)
    )
    return PlotBuilder.plotwindow(;
        title = heading, size = fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], abscissa, ordinate;
            title = heading)
        y = [model.values.density; last(model.values.density)]
        plot = stairs!(axis, model.values.edges, y;
            step = :post, label = "model PDF", kwargs...)
        PlotBuilder.register!(context, axis;
            groups = (model = (plot,),),
            labels = (model = "model PDF",),
            data = ((; xdata = model.values.edges, ydata = y,
                group = :model, label = "model PDF"),))
    end
end

function Makie.ecdfplot(
        result::LineCableModels.MonteCarloResult,
        request = LineCableModels.R;
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        clip::Bool = true,
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
        quantity_units,
        clip)
    heading = title === nothing ?
              _monte_carlo_title(publication, "empirical cumulative distribution") :
              String(title)
    sample = publication.sample
    ordinate = (;
        values = Float64[],
        quantity = LineCableModels.Units.Quantity{:cumulative_probability}(),
        unit = LineCableModels.Units.UnitExpr()
    )
    return PlotBuilder.plotwindow(;
        title = heading, size = fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], sample, ordinate;
            title = heading)
        plot = ecdfplot!(axis, sample.values; label = "empirical CDF", kwargs...)
        PlotBuilder.register!(context, axis;
            groups = (empirical = (plot,),),
            labels = (empirical = "empirical CDF",),
            data = ((; xdata = sample.values, ydata = nothing,
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
        clip::Bool = true,
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
        quantity_units,
        clip)
    heading = title === nothing ?
              _monte_carlo_title(publication, "model cumulative distribution") :
              String(title)
    model = publication.model
    grid = _model_cdf_grid(model.values)
    probability = LineCableModels.UQ.cumulative_probability.(Ref(model.values), grid)
    abscissa = (; values = grid, quantity = model.quantity, unit = model.unit)
    ordinate = (;
        values = probability,
        quantity = LineCableModels.Units.Quantity{:cumulative_probability}(),
        unit = LineCableModels.Units.UnitExpr()
    )
    return PlotBuilder.plotwindow(;
        title = heading, size = fig_size, backend, display_plot, controls,
        export_theme, open_export) do context
        axis = PlotBuilder.axis!(context, context.canvas[1, 1], abscissa, ordinate;
            title = heading)
        plot = lines!(axis, grid, probability; label = "model CDF", kwargs...)
        PlotBuilder.register!(context, axis;
            groups = (model = (plot,),),
            labels = (model = "model CDF",),
            data = ((; xdata = grid, ydata = probability,
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
        clip::Bool = true,
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
        quantity_units,
        clip)
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
    return PlotBuilder.plotwindow(;
        title = heading, size = fig_size, backend, display_plot, controls,
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
            PlotBuilder.register!(context, axis;
                groups = (quantiles = (points,), reference = (line,)),
                labels = (quantiles = "quantiles", reference = "perfect fit"),
                data = (
                    (; xdata = pairs.sample, ydata = pairs.model,
                        group = :quantiles, label = "quantiles"),
                    (; xdata = reference, ydata = reference,
                        group = :reference, label = "perfect fit")
                ))
        else
            PlotBuilder.register!(context, axis;
                groups = (quantiles = (points,),),
                labels = (quantiles = "quantiles",),
                data = ((; xdata = pairs.sample, ydata = pairs.model,
                    group = :quantiles, label = "quantiles"),))
        end
    end
end
