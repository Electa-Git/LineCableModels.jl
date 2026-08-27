const LineParameterDefinition = LineCableModels.Engine.LineParameterPlotDefinition
const LineComparisonDefinition = LineCableModels.Engine.LineParametersBenchmarkPlotDefinition

function _line_panel_state(payload, panel_index::Int)
    payload.runtime.panels === nothing && return (;
        xscale = payload.runtime.xscale,
        yscale = payload.runtime.yscale,
        current_limits = nothing,
        hidden_groups = ()
    )
    return payload.runtime.panels[panel_index]
end

function _line_registration(groups, objects, labels)
    names = Tuple(groups)
    return (
        NamedTuple{names}(Tuple(objects)),
        NamedTuple{names}(Tuple(labels))
    )
end

function _line_axis!(context, payload, panel_index, y_observation, state)
    row, column = payload.positions[panel_index]
    axis = PlotBuilder.axis!(
        context,
        context.canvas[row, column],
        payload.frequency,
        y_observation;
        title = payload.titles[panel_index],
        xscale = state.xscale,
        yscale = state.yscale,
        xscales = payload.xscales[panel_index],
        yscales = payload.yscales[panel_index],
        payload.attributes[panel_index]...
    )
    return axis
end

function _finish_line_axis!(context, axis, groups, objects, labels, data, state)
    registration = only(candidate
    for candidate in context.panels if candidate.axis === axis)
    registered_groups, registered_labels = _line_registration(groups, objects, labels)
    PlotBuilder.register!(context, axis;
        xmetadata = registration.metadata.xaxis,
        ymetadata = registration.metadata.yaxis,
        groups = registered_groups,
        labels = registered_labels,
        data = Tuple(data),
        limits = state.current_limits)
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

function _draw_single_line_page!(context::UIContext, page::PlotPage)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = page.payload
    for panel_index in eachindex(payload.requests)
        state = _line_panel_state(payload, panel_index)
        observation = payload.observations[panel_index]
        rows, columns, _ = payload.coordinates[panel_index]
        y_observation = (;
            values = collect(vec(observation.values)),
            quantity = observation.quantity,
            unit = observation.unit
        )
        axis = _line_axis!(context, payload, panel_index, y_observation, state)
        groups = Symbol[]
        objects = Vector{Any}[]
        labels = String[]
        data = NamedTuple[]
        family = LineCableModels.Units.family(observation.quantity) === Val(:series) ?
                 "series" : "shunt"
        for (local_row, row) in enumerate(rows),
            (local_column, column) in enumerate(columns)

            curve = collect(view(observation.values, local_row, local_column, :))
            label_index = (local_row - 1) * length(columns) + local_column
            label = if payload.legend_labels === nothing
                join(
                    (
                        "$(LineCableModels.Units.symbol(item.quantity))[$row,$column]"
                    for item in payload.observations),
                    ", "
                )
            else
                String(payload.legend_labels[label_index])
            end
            group = Symbol("$(family)_$(row)_$(column)")
            plots = _draw_line!(axis, payload.frequency.values, curve;
                label,
                visible = group ∉ state.hidden_groups,
                attributes = (; linewidth = 2))
            push!(groups, group)
            push!(objects, plots)
            push!(labels, label)
            push!(data, (;
                xdata = payload.frequency.values,
                ydata = curve,
                group,
                label
            ))
        end
        _finish_line_axis!(context, axis, groups, objects, labels, data, state)
    end
    return context
end

function _draw_comparison_page!(context::UIContext, page::PlotPage)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = page.payload
    for panel_index in eachindex(payload.positions)
        state = _line_panel_state(payload, panel_index)
        local_row, local_column = payload.positions[panel_index]
        first_observation = first(payload.observations)
        values = collect(Iterators.flatten(
            view(observation.values, local_row, local_column, :)
        for observation in payload.observations))
        y_observation = (;
            values,
            quantity = first_observation.quantity,
            unit = first_observation.unit
        )
        axis = _line_axis!(context, payload, panel_index, y_observation, state)
        groups = Symbol[]
        objects = Vector{Any}[]
        labels = String[]
        data = NamedTuple[]
        for source_index in eachindex(payload.observations)
            curve = collect(view(
                payload.observations[source_index].values,
                local_row,
                local_column,
                :
            ))
            label = String(payload.legend_labels[source_index])
            group = Symbol("line_parameters_$source_index")
            plots = _draw_line!(axis, payload.frequency.values, curve;
                label,
                visible = group ∉ state.hidden_groups,
                attributes = (;
                    color = payload.colors[source_index],
                    linestyle = :solid,
                    linewidth = 2
                ))
            push!(groups, group)
            push!(objects, plots)
            push!(labels, label)
            push!(data, (;
                xdata = payload.frequency.values,
                ydata = curve,
                group,
                label
            ))
        end
        _finish_line_axis!(context, axis, groups, objects, labels, data, state)
    end
    return context
end

function draw!(
        context::UIContext,
        ::Type{LineParameterDefinition},
        page::PlotPage
)
    return _draw_single_line_page!(context, page)
end

function draw!(
        context::UIContext,
        ::Type{LineComparisonDefinition},
        page::PlotPage
)
    return _draw_comparison_page!(context, page)
end

function _current_line_page(plot::UIPlot)
    original = plot.page.payload
    runtime = (; panels = Tuple(_current_panel_state.(plot.panels)))
    payload = LineCableModels.Engine.LineDashboardPayload(
        original.frequency,
        original.requests,
        original.observations,
        original.coordinates,
        original.positions,
        original.titles,
        original.xscales,
        original.yscales,
        original.attributes,
        original.legend_labels,
        original.colors,
        (; xscale = original.runtime.xscale,
            yscale = original.runtime.yscale,
            panels = runtime.panels)
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
