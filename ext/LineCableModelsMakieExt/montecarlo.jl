const MCDistributionDefinition = LineCableModels.UQ.MCDistributionPlotDefinition

function _distribution_payload(recipe::PlotRecipe, page::PageSpec)
    return recipe.input.pages[page.key.page]
end

function _distribution_state(recipe::PlotRecipe)
    return get(recipe.input,
        :runtime,
        (;
            xscale = :linear,
            yscale = :linear,
            current_limits = nothing,
            hidden_groups = ()
        ))
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

function _draw_distribution_page!(
        context::UIContext,
        recipe::PlotRecipe,
        page::PageSpec
)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = _distribution_payload(recipe, page)
    state = _distribution_state(recipe)
    axis = PlotBuilder.axis!(
        context,
        context.canvas[1, 1],
        payload.x_observation,
        payload.y_observation;
        title = payload.title,
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
        xmetadata = registration.view.xaxis,
        ymetadata = registration.view.yaxis,
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
        recipe::PlotRecipe,
        page::PageSpec
)
    return _draw_distribution_page!(context, recipe, page)
end

function _current_input(plot::UIPlot, ::Type{MCDistributionDefinition})
    panel = only(plot.panels)
    hidden_groups = Tuple(
        group for group in panel.group_order
        if any(plot_object -> !plot_object.visible[], panel.groups[group])
    )
    runtime = (;
        xscale = _current_scale(panel.axis.xscale[]),
        yscale = _current_scale(panel.axis.yscale[]),
        current_limits = _current_limits(panel.axis),
        hidden_groups
    )
    return merge(plot.render.input, (; runtime))
end
