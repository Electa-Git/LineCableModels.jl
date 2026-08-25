const LineParameterDefinition = LineCableModels.Engine.LineParameterPlotDefinition
const LineComparisonDefinition =
    LineCableModels.Engine.LineParametersBenchmarkPlotDefinition

function _line_page(recipe::PlotRecipe, page::PageSpec)
    return recipe.input.pages[page.key.page]
end

function _line_panel_state(recipe::PlotRecipe, page::PageSpec, panel_index::Int, panel)
    runtime = get(recipe.input, :runtime, nothing)
    if runtime !== nothing && runtime.page == page.key.page
        return runtime.panels[panel_index]
    end
    return (;
        xscale = panel.xscale,
        yscale = panel.yscale,
        current_limits = nothing,
        hidden_groups = ()
    )
end

function _line_registry(entries)
    order = Symbol[]
    grouped = Dict{Symbol, Vector{Any}}()
    labels = Dict{Symbol, String}()
    for (curve, plot_objects) in entries
        haskey(grouped, curve.group) || push!(order, curve.group)
        append!(get!(grouped, curve.group, Any[]), plot_objects)
        labels[curve.group] = curve.label
    end
    groups = NamedTuple{Tuple(order)}(Tuple(grouped[group] for group in order))
    group_labels = NamedTuple{Tuple(order)}(Tuple(labels[group] for group in order))
    return groups, group_labels
end

function _draw_line_panel!(
        context::UIContext,
        panel,
        state::NamedTuple
)
    row, column = panel.position
    axis = PlotBuilder.axis!(
        context,
        context.canvas[row, column],
        panel.x_observation,
        panel.y_observation;
        title = panel.title,
        xscale = state.xscale,
        yscale = state.yscale,
        xscales = panel.xscales,
        yscales = panel.yscales,
        panel.attributes...
    )
    registration = only(candidate for candidate in context.panels if candidate.axis === axis)
    entries = Pair{Any, Vector{Any}}[]
    for curve in panel.curves
        visible = curve.group ∉ state.hidden_groups
        plot_objects = _draw_line!(
            axis,
            panel.x_observation.values,
            curve.values;
            label = curve.label,
            visible,
            attributes = curve.style
        )
        push!(entries, curve => plot_objects)
    end
    groups, labels = _line_registry(entries)
    data = Tuple((;
        xdata = panel.x_observation.values,
        ydata = curve.values,
        group = curve.group,
        label = curve.label
    ) for curve in panel.curves)
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
    return registered
end

function _draw_line_page!(
        context::UIContext,
        recipe::PlotRecipe,
        page::PageSpec
)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = _line_page(recipe, page)
    for (panel_index, panel) in enumerate(payload.panels)
        state = _line_panel_state(recipe, page, panel_index, panel)
        _draw_line_panel!(context, panel, state)
    end
    return context
end

function draw!(
        context::UIContext,
        ::Type{LineParameterDefinition},
        recipe::PlotRecipe,
        page::PageSpec
)
    return _draw_line_page!(context, recipe, page)
end

function draw!(
        context::UIContext,
        ::Type{LineComparisonDefinition},
        recipe::PlotRecipe,
        page::PageSpec
)
    return _draw_line_page!(context, recipe, page)
end

function _current_line_input(plot::UIPlot)
    panel_states = Tuple(map(plot.panels) do panel
        hidden_groups = Tuple(
            group for group in panel.group_order
            if any(plot_object -> !plot_object.visible[], panel.groups[group])
        )
        return (;
            xscale = _current_scale(panel.axis.xscale[]),
            yscale = _current_scale(panel.axis.yscale[]),
            current_limits = _current_limits(panel.axis),
            hidden_groups
        )
    end)
    runtime = (; page = plot.page.key.page, panels = panel_states)
    return merge(plot.render.input, (; runtime))
end

function _current_input(plot::UIPlot, ::Type{LineParameterDefinition})
    return _current_line_input(plot)
end

function _current_input(plot::UIPlot, ::Type{LineComparisonDefinition})
    return _current_line_input(plot)
end
