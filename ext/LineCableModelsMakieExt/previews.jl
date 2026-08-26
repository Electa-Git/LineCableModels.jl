const CablePreviewDefinition = LineCableModels.DataModel.CablePreviewPlotDefinition
const SystemPreviewDefinition = LineCableModels.DataModel.SystemPreviewPlotDefinition
const MaterialScaleDefinition = LineCableModels.DataModel.MaterialScalePlotDefinition

_shell_kind(::Type{MaterialScaleDefinition}) = Val(:colorbar)

function _preview_axis!(context::UIContext, title::AbstractString)
    unit = LineCableModels.Units.units(:base, :meter)
    unit_label = label(unit)
    return PlotBuilder.axis!(
        context,
        context.canvas[1, 1];
        title,
        xlabel = "y [$unit_label]",
        ylabel = "z [$unit_label]",
        aspect = :data
    )
end

function _preview_registry(entries)
    order = Symbol[]
    grouped = Dict{Symbol, Vector{Any}}()
    labels = Dict{Symbol, String}()
    for (entry, plot_object) in entries
        haskey(grouped, entry.group) || push!(order, entry.group)
        push!(get!(grouped, entry.group, Any[]), plot_object)
        if hasproperty(entry, :label) && entry.label !== nothing
            labels[entry.group] = entry.label
        end
    end
    groups = NamedTuple{Tuple(order)}(Tuple(grouped[group] for group in order))
    label_order = Tuple(group for group in order if haskey(labels, group))
    group_labels = NamedTuple{label_order}(Tuple(labels[group] for group in label_order))
    return groups, group_labels
end

function _draw_preview!(context::UIContext, page::PlotPage)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = page.payload
    axis = _preview_axis!(context, page.title)
    state = payload.runtime === nothing ?
            (;
        hidden_groups = (),
        current_limits = payload.limits
    ) : payload.runtime
    entries = Pair{Any, Any}[]
    for reference in payload.references
        plot_object = hlines!(
            axis,
            reference.values;
            color = reference.color,
            linewidth = reference.width,
            visible = reference.group ∉ state.hidden_groups
        )
        push!(entries, reference => plot_object)
    end
    for polygon in payload.polygons
        plot_object = poly!(
            axis,
            polygon.geometry;
            label = polygon.label,
            color = polygon.color,
            strokecolor = polygon.stroke,
            strokewidth = polygon.width,
            visible = polygon.group ∉ state.hidden_groups
        )
        push!(entries, polygon => plot_object)
    end
    limits = state.current_limits
    if limits === nothing
        autolimits!(axis)
    else
        xlimits, ylimits = limits
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
    end
    groups, labels = _preview_registry(entries)
    PlotBuilder.register!(
        context,
        axis;
        groups,
        labels,
        limits,
        aspect = :data
    )
    return context
end

function draw!(
        context::UIContext,
        ::Type{CablePreviewDefinition},
        page::PlotPage
)
    return _draw_preview!(context, page)
end

function draw!(
        context::UIContext,
        ::Type{SystemPreviewDefinition},
        page::PlotPage
)
    return _draw_preview!(context, page)
end

function _current_preview_page(plot::UIPlot)
    original = plot.page.payload
    panel = only(plot.panels)
    hidden_groups = Tuple(
        group
    for group in panel.group_order
    if any(plot_object -> !plot_object.visible[], panel.groups[group])
    )
    runtime = (;
        hidden_groups,
        current_limits = _current_limits(panel.axis)
    )
    payload = LineCableModels.DataModel.PreviewPayload(
        original.polygons,
        original.references,
        original.limits,
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

function _replay_page(plot::UIPlot, ::Type{CablePreviewDefinition})
    return _current_preview_page(plot)
end

function _replay_page(plot::UIPlot, ::Type{SystemPreviewDefinition})
    return _current_preview_page(plot)
end

function draw!(
        context::UIContext,
        ::Type{MaterialScaleDefinition},
        ::PlotPage
)
    return context
end
