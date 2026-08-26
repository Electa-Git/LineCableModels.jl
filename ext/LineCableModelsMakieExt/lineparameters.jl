const LineParameterDefinition = LineCableModels.Engine.LineParameterPlotDefinition
const LineComparisonDefinition = LineCableModels.Engine.LineParametersBenchmarkPlotDefinition

function _line_panel_state(payload, panel_index::Int, panel)
    payload.runtime === nothing && return (;
        xscale = panel.xscale,
        yscale = panel.yscale,
        current_limits = nothing,
        hidden_groups = ()
    )
    return payload.runtime.panels[panel_index]
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
    registration = only(candidate
    for candidate in context.panels if candidate.axis === axis)
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
    return registered
end

function _draw_line_page!(context::UIContext, page::PlotPage)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = page.payload
    for (panel_index, panel) in enumerate(payload.panels)
        state = _line_panel_state(payload, panel_index, panel)
        _draw_line_panel!(context, panel, state)
    end
    return context
end

function draw!(
        context::UIContext,
        ::Type{LineParameterDefinition},
        page::PlotPage
)
    return _draw_line_page!(context, page)
end

function draw!(
        context::UIContext,
        ::Type{LineComparisonDefinition},
        page::PlotPage
)
    return _draw_line_page!(context, page)
end

function _current_line_page(plot::UIPlot)
    original = plot.page.payload
    runtime = (; panels = Tuple(_current_panel_state.(plot.panels)))
    payload = LineCableModels.Engine.LinePagePayload(
        original.panels,
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

function _replay_page(plot::UIPlot, ::Type{LineParameterDefinition})
    return _current_line_page(plot)
end

function _replay_page(plot::UIPlot, ::Type{LineComparisonDefinition})
    return _current_line_page(plot)
end
