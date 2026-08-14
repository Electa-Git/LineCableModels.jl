module UIComponents

using Makie
using Measurements
using Dates

import LineCableModels.PlotBuilder
using LineCableModels.PlotBuilder:
                                   AxisSpec, SeriesSpec, ViewSpec, PageSpec,
                                   RenderSpec, UIPlot
const BackendHandler = PlotBuilder.BackendHandler

export build, export_svg

const FIGURE_PADDING = (80, 60, 40, 40)
const LEGEND_WIDTH = 220
const COLORBAR_WIDTH = 160
const EXPORT_TIMESTAMP_FORMAT = "yyyymmdd_HHMMSS"

mutable struct UIContext
    backend::Symbol
    interactive::Bool
    window::Any
    status::Observable{String}
end

struct UIPanel
    view::ViewSpec
    axis::Any
    plots::Vector{Any}
    groups::Dict{Symbol, Vector{Any}}
    group_labels::Dict{Symbol, String}
end

function _theme(; export_mode::Bool = false)
    background = export_mode ? :white : :grey95
    return Theme(
        backgroundcolor = background,
        Axis = (
            titlesize = 15,
            xlabelsize = 14,
            ylabelsize = 14,
            xticklabelsize = 14,
            yticklabelsize = 14
        ),
        Legend = (; fontsize = 14, labelsize = 14),
        Colorbar = (; labelsize = 14, ticklabelsize = 14)
    )
end

function _context(backend, display)
    active = BackendHandler.ensure_backend!(backend)
    interactive = display && active in (:gl, :wgl)
    window = interactive && active === :gl ?
             BackendHandler.make_screen("LineCableModels Plot"; backend = :gl) : nothing
    return UIContext(active, interactive, window, Observable("Ready."))
end

_scale(symbol::Symbol) = symbol === :log10 ? Makie.log10 : Makie.identity

function _numeric_values(values)
    values === nothing && return nothing, nothing
    nominal = Measurements.value.(values)
    errors = Measurements.uncertainty.(values)
    return nominal, any(error -> !iszero(error), errors) ? errors : nothing
end

function _without(attributes::NamedTuple, excluded::Tuple)
    return (; (key => value for (key, value) in pairs(attributes) if key ∉ excluded)...)
end

function _draw!(axis, series::SeriesSpec)
    attributes = _without(series.attributes, (:group,))
    plots = Any[]
    if series.kind === :line
        x, xerror = _numeric_values(series.xdata)
        y, yerror = _numeric_values(series.ydata)
        line = lines!(axis, x, y; label = series.label, attributes...)
        push!(plots, line)
        if yerror !== nothing
            push!(
                plots,
                errorbars!(
                    axis,
                    x,
                    y,
                    yerror;
                    color = :black,
                    direction = :y,
                    whiskerwidth = 3,
                    linewidth = 1
                )
            )
        end
        if xerror !== nothing
            push!(
                plots,
                errorbars!(
                    axis,
                    x,
                    y,
                    xerror;
                    color = :black,
                    direction = :x,
                    whiskerwidth = 3,
                    linewidth = 1
                )
            )
        end
    elseif series.kind === :scatter
        x, _ = _numeric_values(series.xdata)
        y, _ = _numeric_values(series.ydata)
        push!(plots, scatter!(axis, x, y; label = series.label, attributes...))
    elseif series.kind === :histogram
        values, _ = _numeric_values(series.xdata)
        push!(plots, hist!(axis, values; label = series.label, attributes...))
    elseif series.kind === :stairs
        x, _ = _numeric_values(series.xdata)
        y, _ = _numeric_values(series.ydata)
        push!(plots, stairs!(axis, x, y; label = series.label, attributes...))
    elseif series.kind === :heatmap
        push!(
            plots,
            heatmap!(
                axis,
                series.xdata,
                series.ydata,
                series.zdata;
                attributes...
            )
        )
    elseif series.kind === :polygon
        push!(plots, poly!(axis, series.zdata; label = series.label, attributes...))
    elseif series.kind === :hline
        push!(plots, hlines!(axis, series.ydata; label = series.label, attributes...))
    else
        throw(ArgumentError("unsupported PlotBuilder primitive :$(series.kind)"))
    end
    return plots
end

function _axis(parent, view::ViewSpec)
    xaxis = view.xaxis
    yaxis = view.yaxis
    attributes = _without(view.attributes, (:aspect, :limits))
    aspect = get(view.attributes, :aspect, nothing)
    axis = Axis(
        parent;
        xlabel = xaxis === nothing ? "" : xaxis.label,
        ylabel = yaxis === nothing ? "" : yaxis.label,
        title = view.title,
        xscale = xaxis === nothing ? Makie.identity : _scale(xaxis.scale),
        yscale = yaxis === nothing ? Makie.identity : _scale(yaxis.scale),
        aspect = aspect === :data ? DataAspect() : aspect,
        attributes...
    )
    plots = Any[]
    groups = Dict{Symbol, Vector{Any}}()
    group_labels = Dict{Symbol, String}()
    for (index, series) in enumerate(view.series)
        drawn = _draw!(axis, series)
        append!(plots, drawn)
        group = get(series.attributes, :group, Symbol("series_$index"))
        append!(get!(groups, group, Any[]), drawn)
        if series.label !== nothing && !isempty(series.label)
            group_labels[group] = series.label
        end
    end
    if haskey(view.attributes, :limits)
        xlimits, ylimits = view.attributes.limits
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
    end
    return UIPanel(view, axis, plots, groups, group_labels)
end

function _sanitize_filename(value::AbstractString)
    sanitized = lowercase(strip(value))
    sanitized = replace(sanitized, r"[^0-9a-z]+" => "_")
    sanitized = strip(sanitized, '_')
    return isempty(sanitized) ? "linecablemodels_plot" : sanitized
end

function _available_path(page::PageSpec)
    base = _sanitize_filename(get(page.kwargs, :export_name, page.title))
    stamp = Dates.format(Dates.now(), EXPORT_TIMESTAMP_FORMAT)
    candidate = joinpath(pwd(), "$(base)_$(stamp).svg")
    index = 2
    while ispath(candidate)
        candidate = joinpath(pwd(), "$(base)_$(stamp)_$(index).svg")
        index += 1
    end
    return candidate
end

function _visibility_groups(panels)
    groups = Dict{Symbol, Vector{Any}}()
    labels = Dict{Symbol, String}()
    for panel in panels
        for (group, plots) in panel.groups
            append!(get!(groups, group, Any[]), plots)
        end
        merge!(labels, panel.group_labels)
    end
    return groups, labels
end

function _colorbars!(slot, descriptors)
    isempty(descriptors) && return nothing
    grid = GridLayout()
    slot[] = grid
    colsize!(grid, 1, Relative(1))
    for (row, descriptor) in enumerate(descriptors)
        Colorbar(
            grid[row, 1];
            colormap = descriptor.colormap,
            limits = descriptor.limits,
            ticks = descriptor.ticks,
            label = descriptor.label,
            vertical = false,
            width = COLORBAR_WIDTH,
            height = 14
        )
    end
    return grid
end

function _legend!(slot, panels)
    groups, group_labels = _visibility_groups(panels)
    entries = Any[]
    labels = String[]
    for group in sort!(collect(keys(group_labels)); by = string)
        push!(entries, groups[group])
        push!(labels, group_labels[group])
    end
    return isempty(entries) ? nothing : Legend(slot, entries, labels; valign = :top)
end

function _build_page(
        render_spec::RenderSpec,
        page::PageSpec,
        context::UIContext;
        controls::Bool,
        export_mode::Bool
)
    figure = Figure(size = page.size, figure_padding = FIGURE_PADDING)
    toolbar_row = controls ? 1 : 0
    canvas_row = toolbar_row + 1
    status_row = canvas_row + 1
    canvas = GridLayout()
    figure[canvas_row, 1] = canvas
    panel_count = length(page.views)
    columns = max(1, ceil(Int, sqrt(panel_count)))
    panels = UIPanel[]
    for (index, view) in enumerate(page.views)
        row = (index - 1) ÷ columns + 1
        column = (index - 1) % columns + 1
        push!(panels, _axis(canvas[row, column], view))
    end
    side = GridLayout()
    side_column = isempty(page.views) ? 1 : 2
    figure[canvas_row, side_column] = side
    if isempty(page.views)
        colsize!(figure.layout, 1, Relative(1))
        colsize!(side, 1, Relative(1))
    end
    side_row = 1
    legend = nothing
    if get(page.kwargs, :display_legend, true)
        legend = _legend!(side[side_row, 1], panels)
        side_row += 1
    end
    colorbars = get(page.kwargs, :colorbars, ())
    if !isempty(colorbars)
        _colorbars!(side[side_row, 1], colorbars)
        side_row += 1
    end
    isempty(page.views) || colsize!(figure.layout, 2, Fixed(LEGEND_WIDTH))

    widgets = Dict{Symbol, Any}()
    plot_reference = Ref{Any}(nothing)
    if controls
        definitions = get(
            page.kwargs,
            :controls,
            PlotBuilder.control_definitions()
        )
        toolbar = GridLayout()
        figure[1, 1:2] = toolbar
        column = 1
        if definitions.reset
            reset = Button(toolbar[1, column], label = "Reset")
            column += 1
            widgets[:reset] = reset
            on(reset.clicks) do _
                foreach(panel -> autolimits!(panel.axis), panels)
                context.status[] = "Axis limits reset"
            end
        end
        if definitions.export_svg
            save_button = Button(toolbar[1, column], label = "Export SVG")
            column += 1
            widgets[:export_svg] = save_button
            on(save_button.clicks) do _
                try
                    PlotBuilder.export_svg(plot_reference[])
                catch error
                    context.status[] = sprint(showerror, error)
                end
            end
        end
        if definitions.xlog
            active = !isempty(panels) && all(
                panel -> panel.axis.xscale[] === Makie.log10,
                panels
            )
            xlog = Toggle(toolbar[1, column], active = active)
            column += 1
            Label(toolbar[1, column], "log x")
            column += 1
            widgets[:xlog] = xlog
            on(xlog.active) do enabled
                scale = enabled ? Makie.log10 : Makie.identity
                foreach(panel -> panel.axis.xscale[] = scale, panels)
                foreach(panel -> autolimits!(panel.axis), panels)
                context.status[] = enabled ?
                                   "x-axis scale set to log" :
                                   "x-axis scale set to linear"
            end
        end
        if definitions.ylog
            active = !isempty(panels) && all(
                panel -> panel.axis.yscale[] === Makie.log10,
                panels
            )
            ylog = Toggle(toolbar[1, column], active = active)
            column += 1
            Label(toolbar[1, column], "log y")
            widgets[:ylog] = ylog
            on(ylog.active) do enabled
                scale = enabled ? Makie.log10 : Makie.identity
                foreach(panel -> panel.axis.yscale[] = scale, panels)
                foreach(panel -> autolimits!(panel.axis), panels)
                context.status[] = enabled ?
                                   "y-axis scale set to log" :
                                   "y-axis scale set to linear"
            end
        end
        definitions.legend && legend !== nothing && (widgets[:legend] = legend)
        Label(figure[status_row, 1:2], context.status; halign = :left, fontsize = 11)
        rowsize!(figure.layout, 1, Fixed(48))
        rowsize!(figure.layout, canvas_row, Relative(1))
        rowsize!(figure.layout, status_row, Fixed(24))
    end

    built = UIPlot(render_spec, page, figure, panels, widgets, context)
    plot_reference[] = built
    return built
end

function build(
        render_spec::RenderSpec;
        backend = nothing,
        display::Bool = true,
        controls::Bool = true,
        export_mode::Bool = false
)
    context = _context(backend, display)
    built = UIPlot[]
    with_theme(_theme(; export_mode)) do
        for page in render_spec.figures
            plot = _build_page(
                render_spec,
                page,
                context;
                controls,
                export_mode
            )
            push!(built, plot)
            if display
                if context.interactive && context.window !== nothing
                    Base.display(context.window, plot.figure)
                else
                    BackendHandler.renderfig(plot.figure)
                end
            end
        end
    end
    return built
end

function _current_scale(scale)
    scale === Makie.log10 && return :log10
    scale === Makie.identity && return :linear
    throw(ArgumentError("SVG export supports linear and log10 axis scales"))
end

function _axis_with_scale(spec::Union{Nothing, AxisSpec}, scale)
    spec === nothing && return nothing
    return AxisSpec(spec.dim, spec.quantity, spec.units, spec.label, _current_scale(scale))
end

function _current_limits(axis)
    limits = axis.finallimits[]
    xlimits = (limits.origin[1], limits.origin[1] + limits.widths[1])
    ylimits = (limits.origin[2], limits.origin[2] + limits.widths[2])
    return xlimits, ylimits
end

function _current_series(series::SeriesSpec, panel::UIPanel, index::Int)
    group = get(series.attributes, :group, Symbol("series_$index"))
    visible = all(plot_object -> plot_object.visible[], panel.groups[group])
    attributes = merge(series.attributes, (; visible))
    return SeriesSpec(
        series.kind,
        series.xdata,
        series.ydata,
        series.zdata,
        series.label;
        attributes
    )
end

function _current_view(view::ViewSpec, panel::UIPanel)
    series = [_current_series(item, panel, index)
              for (index, item) in enumerate(view.series)]
    attributes = merge(view.attributes, (; limits = _current_limits(panel.axis)))
    return ViewSpec(
        _axis_with_scale(view.xaxis, panel.axis.xscale[]),
        _axis_with_scale(view.yaxis, panel.axis.yscale[]),
        view.zaxis,
        view.title,
        series,
        view.key;
        attributes
    )
end

function _current_page(plot::UIPlot)
    isempty(plot.page.views) && return plot.page
    length(plot.page.views) == length(plot.panels) || throw(
        DimensionMismatch("built panels no longer match the declarative page"),
    )
    views = [_current_view(view, panel)
             for (view, panel) in zip(plot.page.views, plot.panels)]
    return PageSpec(
        plot.page.title, plot.page.size, plot.page.layout, views, plot.page.kwargs)
end

function PlotBuilder.export_svg(plot::UIPlot; path::Union{Nothing, AbstractString} = nothing)
    BackendHandler.backend_available(:cairo) || throw(
        ArgumentError(
        "SVG export requires CairoMakie; load CairoMakie first with `using CairoMakie`",
    ),
    )
    output = path === nothing ? _available_path(plot.page) : abspath(String(path))
    lowercase(splitext(output)[2]) == ".svg" || throw(
        ArgumentError("SVG export paths must use the .svg extension"),
    )
    ispath(output) && throw(ArgumentError("refusing to overwrite existing file: $output"))
    plot.context.status[] = "Exporting SVG..."
    BackendHandler.with_backend(:cairo) do
        one_page = RenderSpec(plot.render.spec, PageSpec[_current_page(plot)])
        exported = build(
            one_page;
            backend = :cairo,
            display = false,
            controls = false,
            export_mode = true
        )
        Makie.save(output, only(exported).figure)
    end
    plot.context.status[] = "Saved SVG to $(basename(output))"
    return output
end

export_svg(args...; kwargs...) = PlotBuilder.export_svg(args...; kwargs...)

end # module UIComponents
