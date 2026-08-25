struct NativeCanvasPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

function _native_observation(observation, name::Symbol)
    observation === nothing && return nothing
    keys(observation) == (:values, :quantity, :unit) || throw(ArgumentError(
        "$name observation must contain only values, quantity, and unit",
    ))
    return observation
end

function _native_exponent(observation)
    observation === nothing && return 0
    finite = Float64[]
    values = observation.values isa Number ? (observation.values,) : observation.values
    for value in values
        numeric = nominal(value)
        numeric isa Real && isfinite(numeric) && push!(finite, abs(Float64(numeric)))
    end
    filter!(!iszero, finite)
    isempty(finite) && return 0
    exponent = floor(Int, log10(maximum(finite)))
    return abs(exponent) >= 3 ? exponent : 0
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
    return registered, Dict{Symbol, String}(
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
    registered_groups, group_labels, group_order = _native_groups(groups, labels)
    series = _native_series(data, registered_groups, group_order)
    plots = Any[]
    for objects in values(registered_groups), object in objects
        object in plots || push!(plots, object)
    end
    view = (;
        xaxis = xmetadata,
        yaxis = ymetadata,
        zaxis = nothing,
        series,
        limits,
        aspect
    )
    panel = UIPanel(
        view,
        axis,
        plots,
        registered_groups,
        group_labels,
        group_order
    )
    index = findfirst(candidate -> candidate.axis === axis, ui.panels)
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
        _native_exponent(x_observation)
    )
    ymetadata = _native_axis_metadata(
        _native_label(y_observation, ylabel),
        yscale,
        yscales,
        _native_exponent(y_observation)
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
    data = x_observation === nothing && y_observation === nothing ? () : ((;
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

function draw!(
        context::UIContext,
        ::Type{NativeCanvasPlotDefinition},
        recipe::PlotRecipe,
        ::PageSpec
)
    context.canvas === nothing && error("the standard shell has no canvas")
    recipe.object(context)
    _refresh_native_panels!(context)
    return context
end

function _native_colorbar(definition::NamedTuple)
    keys(definition) == (:label, :colormap, :limits, :ticks) || throw(ArgumentError(
        "native colorbars must define label, colormap, limits, and ticks",
    ))
    return PlotBuilder.ColorbarSpec(
        definition.label,
        definition.colormap,
        definition.limits,
        definition.ticks
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
    page = PlotBuilder.PageSpec(
        title,
        size,
        (; kind = :native_canvas),
        PlotBuilder.layout_preset(:single, 0),
        PlotBuilder.ViewSpec[];
        legend = PlotBuilder.LegendSpec(enabled = legend),
        colorbars = PlotBuilder.ColorbarSpec[
            _native_colorbar(definition) for definition in colorbars
        ],
        export_spec = PlotBuilder.ExportSpec(
            theme = export_theme,
            name = export_name,
            open_file = open_export
        )
    )
    recipe = PlotBuilder.PlotRecipe(
        NativeCanvasPlotDefinition,
        callback,
        (;),
        (; export_theme, open_export, layout = nothing),
        PlotBuilder.PageSpec[page]
    )
    return only(build(
        recipe;
        backend,
        display = display_plot,
        controls,
        export_mode,
        export_theme
    ))
end
