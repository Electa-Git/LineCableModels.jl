struct NativeCanvasPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
Store a native-canvas callback and the state needed to replay it.
"""
struct NativeCanvasPayload{F, S}
    "Concrete callable invoked with the Makie UI context."
    callback::F
    "Captured runtime state used for current-state SVG replay."
    runtime::S
end

function _native_observation(observation, name::Symbol)
    observation === nothing && return nothing
    keys(observation) == (:values, :quantity, :unit) || throw(ArgumentError(
        "$name observation must contain only values, quantity, and unit",
    ))
    return observation
end

function _common_exponent(observation)
    observation === nothing && return 0
    finite = Float64[]
    values = observation.values isa Number ? (observation.values,) : observation.values
    for value in values
        numeric = nominal(value)
        numeric isa Real && isfinite(numeric) && push!(finite, abs(Float64(numeric)))
    end
    filter!(!iszero, finite)
    isempty(finite) && return 0
    return floor(Int, log10(maximum(finite)))
end

function _native_label(observation, override)
    override === nothing || return String(override)
    observation === nothing && return ""
    return label(observation.quantity, observation.unit)
end

function _native_axis_metadata(label, scale::Symbol, allowed_scales, exponent::Int)
    scale in (:linear, :log10) || throw(ArgumentError(
        "axis scale must be :linear or :log10",
    ))
    allowed_scales isa Tuple && !isempty(allowed_scales) || throw(ArgumentError(
        "allowed axis scales must be a nonempty tuple",
    ))
    all(item -> item in (:linear, :log10), allowed_scales) ||
        throw(ArgumentError(
            "allowed axis scales must contain only :linear and :log10",
        ))
    scale in allowed_scales || throw(ArgumentError(
        "current axis scale must occur in its allowed scales",
    ))
    return (; label, scale, allowed_scales, exponent)
end

function _native_groups(groups::NamedTuple, labels::NamedTuple)
    unknown_labels = Tuple(key for key in keys(labels) if key ∉ keys(groups))
    isempty(unknown_labels) || throw(ArgumentError(
        "legend labels do not match registered groups: $(join(unknown_labels, ", "))",
    ))
    registered = Dict{Symbol, Vector{Any}}()
    for (name, objects) in pairs(groups)
        values = objects isa Union{Tuple, AbstractVector} ? Any[objects...] : Any[objects]
        registered[name] = values
    end
    return registered,
    Dict{Symbol, String}(
        name => String(value) for (name, value) in pairs(labels)
    ), Symbol[keys(groups)...]
end

function _native_series(data::Tuple, groups::Dict{Symbol, Vector{Any}}, order)
    series = Any[]
    for (index, entry) in enumerate(data)
        entry isa NamedTuple || throw(ArgumentError(
            "registered axis data entries must be named tuples",
        ))
        required = (:xdata, :ydata)
        all(name -> haskey(entry, name), required) || throw(ArgumentError(
            "registered axis data entries must define xdata and ydata",
        ))
        group = get(entry, :group, Symbol("series_$index"))
        group isa Symbol || throw(ArgumentError("registered plot groups must be symbols"))
        haskey(groups, group) || begin
            groups[group] = Any[]
            push!(order, group)
        end
        push!(series, (;
            xdata = entry.xdata,
            ydata = entry.ydata,
            group,
            label = get(entry, :label, nothing)
        ))
    end
    return series
end

function PlotBuilder.register!(
        ui::UIContext,
        axis;
        xmetadata = nothing,
        ymetadata = nothing,
        groups::NamedTuple = (;),
        labels::NamedTuple = (;),
        data::Tuple = (),
        limits = nothing,
        aspect = nothing
)
    index = findfirst(candidate -> candidate.axis === axis, ui.panels)
    previous = index === nothing ? nothing : ui.panels[index]
    if previous !== nothing
        xmetadata === nothing && (xmetadata = previous.metadata.xaxis)
        ymetadata === nothing && (ymetadata = previous.metadata.yaxis)
        limits === nothing && (limits = previous.metadata.limits)
        aspect === nothing && (aspect = previous.metadata.aspect)
    end
    registered_groups, group_labels, group_order = _native_groups(groups, labels)
    series = _native_series(data, registered_groups, group_order)
    plots = Any[]
    for objects in values(registered_groups), object in objects

        object in plots || push!(plots, object)
    end
    metadata = (;
        xaxis = xmetadata,
        yaxis = ymetadata,
        series,
        limits,
        aspect
    )
    panel = UIPanel(
        metadata,
        axis,
        plots,
        registered_groups,
        group_labels,
        group_order
    )
    if index === nothing
        push!(ui.panels, panel)
    else
        ui.panels[index] = panel
    end
    xmetadata === nothing || _set_axis_scale!(
        axis, xmetadata, :x, xmetadata.exponent, xmetadata.scale)
    ymetadata === nothing || _set_axis_scale!(
        axis, ymetadata, :y, ymetadata.exponent, ymetadata.scale)
    return axis
end

function PlotBuilder.axis!(
        ui::UIContext,
        position,
        x_observation = nothing,
        y_observation = nothing;
        title::AbstractString = "",
        xlabel = nothing,
        ylabel = nothing,
        xscale::Symbol = :linear,
        yscale::Symbol = :linear,
        xscales::Tuple = (xscale,),
        yscales::Tuple = (yscale,),
        aspect = nothing,
        kwargs...
)
    x_observation = _native_observation(x_observation, :x)
    y_observation = _native_observation(y_observation, :y)
    xmetadata = _native_axis_metadata(
        _native_label(x_observation, xlabel),
        xscale,
        xscales,
        _common_exponent(x_observation)
    )
    ymetadata = _native_axis_metadata(
        _native_label(y_observation, ylabel),
        yscale,
        yscales,
        _common_exponent(y_observation)
    )
    axis = Axis(
        position;
        xlabel = _axis_label(xmetadata, xmetadata.exponent, xscale),
        ylabel = _axis_label(ymetadata, ymetadata.exponent, yscale),
        title = String(title),
        xscale = _scale(xscale),
        yscale = _scale(yscale),
        xticks = _ticks(xscale),
        yticks = _ticks(yscale),
        xtickformat = _tickformat(xmetadata.exponent, xscale),
        ytickformat = _tickformat(ymetadata.exponent, yscale),
        aspect = aspect === :data ? DataAspect() : aspect,
        tellwidth = false,
        tellheight = false,
        kwargs...
    )
    data = x_observation === nothing && y_observation === nothing ? () :
           ((;
        xdata = x_observation === nothing ? nothing : x_observation.values,
        ydata = y_observation === nothing ? nothing : y_observation.values,
        group = :axis_values,
        label = nothing
    ),)
    PlotBuilder.register!(
        ui,
        axis;
        xmetadata,
        ymetadata,
        data,
        aspect
    )
    return axis
end

function _refresh_native_panels!(ui::UIContext)
    for panel in ui.panels
        isempty(panel.plots) || continue
        append!(panel.plots, Any[panel.axis.scene.plots...])
    end
    return ui
end

function _apply_panel_state!(panel::UIPanel, state)
    xmetadata = panel.metadata.xaxis
    ymetadata = panel.metadata.yaxis
    xmetadata === nothing || _set_axis_scale!(
        panel.axis, xmetadata, :x, xmetadata.exponent, state.xscale)
    ymetadata === nothing || _set_axis_scale!(
        panel.axis, ymetadata, :y, ymetadata.exponent, state.yscale)
    for group in panel.group_order, object in panel.groups[group]

        object.visible[] = group ∉ state.hidden_groups
    end
    if state.current_limits === nothing
        _reset_panel_limits!(panel)
    else
        xlimits, ylimits = state.current_limits
        xlims!(panel.axis, xlimits...)
        ylims!(panel.axis, ylimits...)
    end
    return panel
end

function _apply_native_runtime!(context::UIContext, runtime)
    runtime === nothing && return context
    length(runtime.panels) == length(context.panels) || throw(DimensionMismatch(
        "native-canvas replay produced a different number of registered axes",
    ))
    foreach(_apply_panel_state!, context.panels, runtime.panels)
    return context
end

function draw!(
        context::UIContext,
        ::Type{NativeCanvasPlotDefinition},
        page::PlotPage
)
    context.canvas === nothing && error("the standard shell has no canvas")
    payload = page.payload
    payload.callback(context)
    _refresh_native_panels!(context)
    _apply_native_runtime!(context, payload.runtime)
    return context
end

function _native_colorbar(definition::NamedTuple)
    keys(definition) == (:label, :colormap, :limits, :ticks) || throw(ArgumentError(
        "native colorbars must define label, colormap, limits, and ticks",
    ))
    return ColorbarDefinition(
        definition.label,
        definition.colormap,
        definition.limits,
        definition.ticks
    )
end

function _current_panel_state(panel::UIPanel)
    hidden_groups = Tuple(
        group
    for group in panel.group_order
    if any(plot_object -> !plot_object.visible[], panel.groups[group])
    )
    return (;
        xscale = _current_scale(panel.axis.xscale[]),
        yscale = _current_scale(panel.axis.yscale[]),
        current_limits = _current_limits(panel.axis),
        hidden_groups
    )
end

function _replay_page(plot::UIPlot, ::Type{NativeCanvasPlotDefinition})
    original = plot.page.payload
    runtime = (; panels = Tuple(_current_panel_state.(plot.panels)))
    payload = NativeCanvasPayload(
        original.callback,
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

function PlotBuilder.plotwindow(
        callback::F;
        title::AbstractString,
        size::Tuple{<:Integer, <:Integer} = (800, 400),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        legend::Bool = true,
        colorbars::Tuple = (),
        export_theme::Symbol = :default,
        open_export::Bool = true,
        export_name::AbstractString = title,
        export_mode::Bool = false
) where {F}
    payload = NativeCanvasPayload(
        callback,
        nothing
    )
    page = PlotPage(
        title,
        size,
        (; kind = :native_canvas),
        payload;
        legend = LegendDefinition(enabled = legend),
        colorbars = Tuple(_native_colorbar(definition) for definition in colorbars),
        export_definition = ExportDefinition(
            theme = export_theme,
            name = export_name,
            open_file = open_export
        )
    )
    recipe = PlotRecipe(NativeCanvasPlotDefinition, (page,))
    return only(build(
        recipe;
        backend,
        display = display_plot,
        controls,
        export_mode,
        export_theme
    ))
end
