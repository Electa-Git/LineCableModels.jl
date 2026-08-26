const MCDistributionDefinition = LineCableModels.UQ.MCDistributionPlotDefinition

function _distribution_state(payload)
    return payload.runtime === nothing ?
           (;
        xscale = :linear,
        yscale = :linear,
        current_limits = nothing,
        hidden_groups = ()
    ) : payload.runtime
end

function _draw_distribution_layer!(axis, ::Val{:histogram}, layer, visible)
    values, _ = _numeric_values(layer.x)
    return Any[hist!(
        axis,
        values;
        label = layer.label,
        visible,
        layer.style...
    )]
end

function _draw_distribution_layer!(axis, ::Val{:stairs}, layer, visible)
    x, _ = _numeric_values(layer.x)
    y, _ = _numeric_values(layer.y)
    return Any[stairs!(
        axis,
        x,
        y;
        label = layer.label,
        visible,
        layer.style...
    )]
end

function _draw_distribution_layer!(axis, ::Val{:line}, layer, visible)
    return _draw_line!(
        axis,
        layer.x,
        layer.y;
        label = layer.label,
        visible,
        attributes = layer.style
    )
end

function _draw_distribution_layer!(axis, ::Val{:scatter}, layer, visible)
    x, _ = _numeric_values(layer.x)
    y, _ = _numeric_values(layer.y)
    return Any[scatter!(
        axis,
        x,
        y;
        label = layer.label,
        visible,
        layer.style...
    )]
end

function _draw_distribution_page!(context::UIContext, page::PlotPage)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = page.payload
    state = _distribution_state(payload)
    axis = PlotBuilder.axis!(
        context,
        context.canvas[1, 1],
        payload.x_observation,
        payload.y_observation;
        title = page.title,
        xlabel = payload.xlabel,
        ylabel = payload.ylabel,
        xscale = state.xscale,
        yscale = state.yscale,
        xscales = (:linear,),
        yscales = (:linear,)
    )
    registration = only(
        candidate for candidate in context.panels if candidate.axis === axis
    )
    plot_objects = Tuple(map(payload.layers) do layer
        _draw_distribution_layer!(
            axis,
            layer.kind,
            layer,
            layer.group ∉ state.hidden_groups
        )
    end)
    group_names = Tuple(layer.group for layer in payload.layers)
    groups = NamedTuple{group_names}(plot_objects)
    labels = NamedTuple{group_names}(Tuple(layer.label for layer in payload.layers))
    data = Tuple((;
                     xdata = layer.x,
                     ydata = layer.y,
                     group = layer.group,
                     label = layer.label
                 ) for layer in payload.layers)
    PlotBuilder.register!(
        context,
        axis;
        xmetadata = registration.metadata.xaxis,
        ymetadata = registration.metadata.yaxis,
        groups,
        labels,
        data,
        limits = state.current_limits
    )
    registered = only(candidate for candidate in context.panels if candidate.axis === axis)
    if state.current_limits === nothing
        _reset_panel_limits!(registered)
    else
        xlimits, ylimits = state.current_limits
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
    end
    return context
end

function draw!(
        context::UIContext,
        ::Type{MCDistributionDefinition},
        page::PlotPage
)
    return _draw_distribution_page!(context, page)
end

function _replay_page(plot::UIPlot, ::Type{MCDistributionDefinition})
    original = plot.page.payload
    panel = only(plot.panels)
    runtime = _current_panel_state(panel)
    payload = LineCableModels.UQ.MCDistributionPagePayload(
        original.x_observation,
        original.y_observation,
        original.xlabel,
        original.ylabel,
        original.layers,
        runtime
    )
    return PlotPage(
        plot.page.title,
        plot.page.size,
        plot.page.key,
        payload;
        legend = plot.page.legend,
        colorbars = plot.page.colorbars,
        widgets = plot.page.widgets,
        export_definition = plot.page.export_definition
    )
end
