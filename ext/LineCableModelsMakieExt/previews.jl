const CablePreviewDefinition = LineCableModels.DataModel.CablePreviewPlotDefinition
const SystemPreviewDefinition = LineCableModels.DataModel.SystemPreviewPlotDefinition
const MaterialScaleDefinition = LineCableModels.DataModel.MaterialScalePlotDefinition

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

function _draw_preview!(context::UIContext, payload, state::NamedTuple)
    context.canvas === nothing && error("the standard shell has no canvas")
    axis = _preview_axis!(context, payload.title)
    hidden_groups = get(state, :hidden_groups, ())
    entries = Pair{Any, Any}[]
    for reference in payload.references
        plot_object = hlines!(
            axis,
            reference.values;
            color = reference.color,
            linewidth = reference.width,
            visible = reference.group ∉ hidden_groups
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
            visible = polygon.group ∉ hidden_groups
        )
        push!(entries, polygon => plot_object)
    end
    limits = get(state, :current_limits, payload.limits)
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
        recipe::PlotRecipe,
        ::PageSpec
)
    return _draw_preview!(context, recipe.input.payload, recipe.input)
end

function draw!(
        context::UIContext,
        ::Type{SystemPreviewDefinition},
        recipe::PlotRecipe,
        ::PageSpec
)
    return _draw_preview!(context, recipe.input.payload, recipe.input)
end

function _current_preview_input(plot::UIPlot)
    panel = only(plot.panels)
    hidden_groups = Tuple(
        group for group in panel.group_order
        if any(plot_object -> !plot_object.visible[], panel.groups[group])
    )
    current_limits = _current_limits(panel.axis)
    return merge(plot.render.input, (; hidden_groups, current_limits))
end

function _current_input(plot::UIPlot, ::Type{CablePreviewDefinition})
    return _current_preview_input(plot)
end

function _current_input(plot::UIPlot, ::Type{SystemPreviewDefinition})
    return _current_preview_input(plot)
end

function draw!(
        context::UIContext,
        ::Type{MaterialScaleDefinition},
        ::PlotRecipe,
        ::PageSpec
)
    return context
end
