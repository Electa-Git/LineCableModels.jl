const CablePreviewDefinition = LineCableModels.DataModel.CablePreviewPlotDefinition
const CableCollectionPreviewDefinition =
    LineCableModels.DataModel.CableCollectionPreviewPlotDefinition
const SystemPreviewDefinition = LineCableModels.DataModel.SystemPreviewPlotDefinition
const MaterialScaleDefinition = LineCableModels.DataModel.MaterialScalePlotDefinition

function build_shell(
        context::UIContext,
        ::Type{MaterialScaleDefinition},
        page::PlotPage
)
    return _colorbar_shell!(context, page)
end

function _preview_axis!(context::UIContext, position, title::AbstractString)
    unit = LineCableModels.Units.units(:base, :meter)
    unit_label = label(unit)
    return PlotBuilder.axis!(
        context,
        position;
        title,
        xlabel = "y [$unit_label]",
        ylabel = "z [$unit_label]",
        aspect = :data
    )
end

_preview_label(entry::LineCableModels.DataModel.PreviewPolygon) = entry.label
_preview_label(::LineCableModels.DataModel.PreviewReferenceLine) = nothing

function _preview_registry(entries)
    order = Symbol[]
    grouped = Dict{Symbol, Vector{Any}}()
    labels = Dict{Symbol, String}()
    for (entry, plot_object) in entries
        haskey(grouped, entry.group) || push!(order, entry.group)
        push!(get!(grouped, entry.group, Any[]), plot_object)
        entry_label = _preview_label(entry)
        if entry_label !== nothing
            labels[entry.group] = entry_label
        end
    end
    groups = NamedTuple{Tuple(order)}(Tuple(grouped[group] for group in order))
    label_order = Tuple(group for group in order if haskey(labels, group))
    group_labels = NamedTuple{label_order}(Tuple(labels[group] for group in label_order))
    return groups, group_labels
end

function _draw_preview!(context::UIContext, position, title, payload)
    context.canvas === nothing && error("the standard shell has no canvas")
    axis = _preview_axis!(context, position, title)
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

function _draw_preview!(context::UIContext, page::PlotPage)
    return _draw_preview!(context, context.canvas[1, 1], page.title, page.payload)
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

function draw!(
        context::UIContext,
        ::Type{CableCollectionPreviewDefinition},
        page::PlotPage
)
    context.canvas === nothing && error("the standard shell has no canvas")
    rows, columns = page.payload.layout

    # The detached recipe has already selected the grid. Equal relative tracks
    # make occupied and deliberately empty slots share the same cell geometry.
    # A sized child grid is necessary because the shell canvas itself begins
    # without tracks and cannot assign a size to an unoccupied trailing slot.
    grid = GridLayout(rows, columns)
    grid.default_rowgap = Fixed(GRID_ROW_GAP)
    grid.default_colgap = Fixed(GRID_COLUMN_GAP)
    context.canvas[1, 1] = grid
    for row in 1:rows
        rowsize!(grid, row, Relative(1 / rows))
    end
    for column in 1:columns
        colsize!(grid, column, Relative(1 / columns))
    end

    # Each detached scalar payload follows the same drawing path as
    # `preview(design)`. This specialization only chooses the Makie grid slot
    # and uses the cable identifier already fixed as the panel title.
    for panel in page.payload.panels
        row, column = panel.position
        _draw_preview!(
            context,
            grid[row, column],
            panel.title,
            panel.payload
        )
    end
    return context
end

function place_colorbars!(
        context::UIContext,
        ::Type{CableCollectionPreviewDefinition},
        page::PlotPage
)
    # This recipe deliberately publishes no legend, so the generic
    # definition-dispatched legend stage collapses that dock row. If a later
    # collection recipe publishes one, its renderer can specialize
    # `place_legend!` on this same definition type and either call the standard
    # two-argument method or build its own Makie legend.

    # `PlotBuilder.fetch` owns which material scales belong to this page. This
    # renderer method owns only their Makie placement, so the detached recipe
    # remains independent of shell coordinates and backend objects.
    isempty(page.colorbars) && return place_colorbars!(context, page)

    # With the legend row collapsed, the standard shell's auto-sized side grid
    # would otherwise remain vertically centered. This recipe aligns that grid
    # to the top before reusing the standard horizontal colorbar implementation.
    # Another definition can overload the same method and select another slot,
    # orientation, or arrangement.
    if context.shell.kind === :standard
        context.shell.side.valign[] = :top
    end
    return place_colorbars!(context, page)
end

function _current_preview_payload(original, panel)
    hidden_groups = Tuple(
        group
    for group in panel.group_order
    if any(plot_object -> !plot_object.visible[], panel.groups[group])
    )
    runtime = (;
        hidden_groups,
        current_limits = _current_limits(panel.axis)
    )
    return LineCableModels.DataModel.PreviewPayload(
        original.polygons,
        original.references,
        original.limits,
        runtime
    )
end

function _current_preview_page(plot::UIPlot)
    original = plot.page.payload
    panel = only(plot.panels)
    payload = _current_preview_payload(original, panel)
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

function _replay_page(plot::UIPlot, ::Type{CableCollectionPreviewDefinition})
    original = plot.page.payload
    length(original.panels) == length(plot.panels) || error(
        "rendered cable-preview panels no longer match the declarative page",
    )

    # SVG replay preserves the current limits of every subplot independently.
    # The collection payload itself remains the same NamedTuple grammar used by
    # the detached recipe; only each existing `PreviewPayload.runtime` changes.
    panels = map(enumerate(original.panels)) do (index, source_panel)
        payload = _current_preview_payload(source_panel.payload, plot.panels[index])
        return merge(source_panel, (; payload))
    end
    payload = (; panels, layout = original.layout)
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
