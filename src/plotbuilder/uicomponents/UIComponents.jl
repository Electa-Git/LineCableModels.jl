module UIComponents

using Makie
using Measurements
using Dates
using Printf: @sprintf

import LineCableModels.PlotBuilder
using LineCableModels.PlotBuilder:
                                   AxisSpec, SeriesSpec, ViewSpec, PageSpec,
                                   RenderSpec, UIPlot
const BackendHandler = PlotBuilder.BackendHandler

export build, export_svg

const FIGURE_PADDING = (20, 20, 28, 28)
const COLORBAR_WIDTH = 140
const COLORBAR_TICK_LABEL_SIZE = 12
const COLORBAR_LABEL_SIZE = 14
const COLORBAR_END_PADDING = 28
const GRID_ROW_GAP = 6
const GRID_COLUMN_GAP = 6
const LEGEND_GAP = 4
const TOOLBAR_HEIGHT = 36
const STATUSBAR_HEIGHT = 20
const BUTTON_SIZE = 32
const BUTTON_ICON_SIZE = 18
const BACKGROUND_INTERACTIVE = :grey90
const BACKGROUND_EXPORT = :white
const BUTTON_BACKGROUND = Makie.RGBf(0.94, 0.94, 0.94)
const ICON_COLOR = Makie.RGBAf(0.15, 0.15, 0.15, 1.0)
const MI_REFRESH = "\uE5D5"
const MI_SAVE = "\uE161"
const ICON_FONT = joinpath(
    pkgdir(PlotBuilder),
    "assets",
    "fonts",
    "material-icons",
    "MaterialIcons-Regular.ttf"
)
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
    group_order::Vector{Symbol}
end

function _theme(; export_mode::Bool = false)
    background = export_mode ? BACKGROUND_EXPORT : BACKGROUND_INTERACTIVE
    return Theme(
        backgroundcolor = background,
        fonts = (; icons = ICON_FONT),
        Axis = (
            titlesize = 15,
            xlabelsize = 14,
            ylabelsize = 14,
            xticklabelsize = 14,
            yticklabelsize = 14,
            xminorgridvisible = false,
            yminorgridvisible = false,
            xminorticksvisible = false,
            yminorticksvisible = false
        ),
        Button = (; buttoncolor = BUTTON_BACKGROUND),
        Legend = (; fontsize = 14, labelsize = 14),
        Colorbar = (; labelsize = 14, ticklabelsize = 14)
    )
end

function _context(active::Symbol, display::Bool, title::AbstractString)
    interactive = display && active in (:gl, :wgl)
    window = interactive && active === :gl ?
             BackendHandler.make_screen(
        "Fig. $(BackendHandler.next_fignum()) – $title";
        backend = :gl
    ) : nothing
    return UIContext(active, interactive, window, Observable("Ready."))
end

_scale(symbol::Symbol) = symbol === :log10 ? Makie.log10 : Makie.identity

function _linear_tickformat(exponent::Int)
    scale = 10.0^exponent
    return values -> [@sprintf("%.4g", value / scale) for value in values]
end

function _decade_ticks(vmin, vmax)
    isfinite(vmin) && isfinite(vmax) && 0 < vmin <= vmax || return (Float64[], String[])
    first_exponent = ceil(Int, log10(vmin))
    last_exponent = floor(Int, log10(vmax))
    first_exponent > last_exponent && return (Float64[], String[])
    exponents = first_exponent:last_exponent
    values = 10.0 .^ exponents
    labels = [Makie.rich(
                  "10",
                  Makie.superscript(
                      replace(string(exponent), "-" => "−");
                      offset = Makie.Vec2f(0.1, 0.0)
                  )
              ) for exponent in exponents]
    return values, labels
end

function _axis_label(spec::Union{Nothing, AxisSpec}, exponent::Int, scale::Symbol)
    spec === nothing && return ""
    scale === :log10 && return spec.label
    iszero(exponent) && return spec.label
    formatted_exponent = replace(string(exponent), "-" => "−")
    return Makie.rich(
        spec.label,
        "  × 10",
        Makie.superscript(
            formatted_exponent;
            offset = Makie.Vec2f(0.1, 0.0)
        )
    )
end

function _tickformat(exponent::Int, scale::Symbol)
    return scale === :log10 ? Makie.automatic : _linear_tickformat(exponent)
end

_ticks(scale::Symbol) = scale === :log10 ? _decade_ticks : Makie.automatic

function _set_axis_scale!(
        axis, spec::Union{Nothing, AxisSpec}, dim::Symbol, exponent::Int, scale::Symbol)
    ticks = _ticks(scale)
    formatter = _tickformat(exponent, scale)
    label = _axis_label(spec, exponent, scale)
    if dim === :x
        axis.xticks[] = ticks
        axis.xtickformat[] = formatter
        axis.xlabel[] = label
        axis.xscale[] = _scale(scale)
    elseif dim === :y
        axis.yticks[] = ticks
        axis.ytickformat[] = formatter
        axis.ylabel[] = label
        axis.yscale[] = _scale(scale)
    else
        throw(ArgumentError("axis dimension must be :x or :y"))
    end
    return axis
end

function _axis_values(view::ViewSpec, dim::Symbol; include_uncertainty::Bool = false)
    values = Float64[]
    for series in view.series
        data = dim === :x ? series.xdata : series.ydata
        data === nothing && continue
        for sample in data
            nominal = Measurements.value(sample)
            nominal isa Real || continue
            numeric = Float64(nominal)
            isfinite(numeric) || continue
            uncertainty = abs(Float64(Measurements.uncertainty(sample)))
            if include_uncertainty && isfinite(uncertainty) && !iszero(uncertainty)
                push!(values, numeric - uncertainty, numeric + uncertainty)
            else
                push!(values, numeric)
            end
        end
    end
    return values
end

function _nearly_constant(values)
    isempty(values) && return false
    lower, upper = extrema(values)
    scale = max(abs(lower), abs(upper), floatmin(Float64))
    return upper - lower <= 64eps(Float64) * scale
end

function _linear_constant_limits(values, interval_values)
    center = sum(extrema(values)) / 2
    base_halfspan = iszero(center) ? 1.0 : 0.05abs(center)
    interval_halfspan = maximum(abs(value - center) for value in interval_values)
    halfspan = max(base_halfspan, 2interval_halfspan)
    return center - halfspan, center + halfspan
end

function _log_decade_limits(values)
    all(>(0), values) || throw(
        DomainError(values, "logarithmic axes require strictly positive data"),
    )
    lower, upper = extrema(values)
    lower_exponent = floor(Int, log10(lower))
    upper_exponent = ceil(Int, log10(upper))
    if lower_exponent == upper_exponent
        lower_exponent -= 1
        upper_exponent += 1
    end
    return 10.0^lower_exponent, 10.0^upper_exponent
end

function _reset_panel_limits!(panel::UIPanel)
    axis = panel.axis
    view = panel.view
    autolimits!(axis)
    haskey(view.attributes, :limits) && return axis
    for dim in (:x, :y)
        values = _axis_values(view, dim)
        isempty(values) && continue
        scale = dim === :x ? axis.xscale[] : axis.yscale[]
        _nearly_constant(values) || continue
        interval_values = _axis_values(view, dim; include_uncertainty = true)
        limits = scale === Makie.log10 ?
                 _log_decade_limits(interval_values) :
                 _linear_constant_limits(values, interval_values)
        dim === :x ? xlims!(axis, limits...) : ylims!(axis, limits...)
    end
    return axis
end

function _icon_label(glyph::AbstractString)
    return Makie.rich(
        glyph;
        font = :icons,
        fontsize = BUTTON_ICON_SIZE,
        color = ICON_COLOR,
        offset = (0, -0.18)
    )
end

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

function _axis(parent, view::ViewSpec, page::PageSpec)
    xaxis = view.xaxis
    yaxis = view.yaxis
    x_exponent = get(page.kwargs, :x_exponent, 0)
    y_exponent = get(page.kwargs, :y_exponent, 0)
    xscale = xaxis === nothing ? :linear : xaxis.scale
    yscale = yaxis === nothing ? :linear : yaxis.scale
    attributes = _without(view.attributes, (:aspect, :limits))
    aspect = get(view.attributes, :aspect, nothing)
    axis = Axis(
        parent;
        xlabel = _axis_label(xaxis, x_exponent, xscale),
        ylabel = _axis_label(yaxis, y_exponent, yscale),
        title = view.title,
        xscale = _scale(xscale),
        yscale = _scale(yscale),
        xticks = _ticks(xscale),
        yticks = _ticks(yscale),
        xtickformat = _tickformat(x_exponent, xscale),
        ytickformat = _tickformat(y_exponent, yscale),
        aspect = aspect === :data ? DataAspect() : aspect,
        attributes...
    )
    plots = Any[]
    groups = Dict{Symbol, Vector{Any}}()
    group_labels = Dict{Symbol, String}()
    group_order = Symbol[]
    for (index, series) in enumerate(view.series)
        drawn = _draw!(axis, series)
        append!(plots, drawn)
        group = get(series.attributes, :group, Symbol("series_$index"))
        haskey(groups, group) || push!(group_order, group)
        append!(get!(groups, group, Any[]), drawn)
        if series.label !== nothing && !isempty(series.label)
            group_labels[group] = series.label
        end
    end
    if haskey(view.attributes, :limits)
        xlimits, ylimits = view.attributes.limits
        xlims!(axis, xlimits...)
        ylims!(axis, ylimits...)
    else
        _reset_panel_limits!(UIPanel(view, axis, plots, groups, group_labels, group_order))
    end
    return UIPanel(view, axis, plots, groups, group_labels, group_order)
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
    order = Symbol[]
    for panel in panels
        for group in panel.group_order
            haskey(groups, group) || push!(order, group)
            plots = panel.groups[group]
            append!(get!(groups, group, Any[]), plots)
        end
        merge!(labels, panel.group_labels)
    end
    return groups, labels, order
end

function _rich_scientific_label(label::AbstractString)
    matched = match(r"^([+−-]?[0-9]+(?:\.[0-9]+)?)e([+−-]?[0-9]+)$", label)
    matched === nothing && return label
    coefficient, raw_exponent = matched.captures
    exponent = parse(Int, replace(raw_exponent, "−" => "-"))
    formatted_exponent = replace(string(exponent), "-" => "−")
    prefix = coefficient == "1" ? "" : coefficient == "-1" ? "−" : "$coefficient×"
    return Makie.rich(
        prefix,
        "10",
        Makie.superscript(
            formatted_exponent;
            offset = Makie.Vec2f(0.1, 0.0)
        )
    )
end

function _colorbar_ticks(ticks)
    values, labels = ticks
    return values, map(_rich_scientific_label, labels)
end

function _colorbars!(slot, descriptors)
    isempty(descriptors) && return nothing
    grid = GridLayout()
    grid.default_colgap = Fixed(0)
    slot[] = grid
    colsize!(grid, 1, Fixed(COLORBAR_WIDTH))
    for (row, descriptor) in enumerate(descriptors)
        Colorbar(
            grid[row, 1];
            colormap = descriptor.colormap,
            limits = descriptor.limits,
            ticks = _colorbar_ticks(descriptor.ticks),
            label = descriptor.label,
            vertical = false,
            width = COLORBAR_WIDTH,
            height = 14,
            ticklabelsize = COLORBAR_TICK_LABEL_SIZE,
            labelsize = COLORBAR_LABEL_SIZE
        )
        Label(grid[row, 2], ""; tellheight = false)
    end
    colsize!(grid, 2, Fixed(COLORBAR_END_PADDING))
    return grid
end

function _legend!(slot, panels)
    groups, group_labels, group_order = _visibility_groups(panels)
    entries = Any[]
    labels = String[]
    for group in group_order
        haskey(group_labels, group) || continue
        push!(entries, groups[group])
        push!(labels, group_labels[group])
    end
    return isempty(entries) ? nothing :
           Legend(
        slot,
        entries,
        labels;
        halign = :right,
        valign = :top
    )
end

function _build_page(
        render_spec::RenderSpec,
        page::PageSpec,
        context::UIContext;
        controls::Bool,
        export_mode::Bool
)
    figure = Figure(size = page.size, figure_padding = FIGURE_PADDING)
    figure.layout.default_rowgap = Fixed(GRID_ROW_GAP)
    figure.layout.default_colgap = Fixed(GRID_COLUMN_GAP)
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
        push!(panels, _axis(canvas[row, column], view, page))
    end
    side = GridLayout()
    side.default_rowgap = Fixed(LEGEND_GAP)
    side.halign = :right
    side_column = isempty(page.views) ? 1 : 2
    figure[canvas_row, side_column] = side
    if isempty(page.views)
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
    isempty(page.views) || colsize!(figure.layout, 2, Auto(true))

    widgets = Dict{Symbol, Any}()
    plot_reference = Ref{Any}(nothing)
    if controls
        definitions = get(
            page.kwargs,
            :controls,
            PlotBuilder.control_definitions()
        )
        toolbar = GridLayout()
        toolbar.default_colgap = Fixed(4)
        toolbar.halign = :left
        toolbar.valign = :bottom
        figure[1, 1:2] = toolbar
        column = 1
        if definitions.reset
            reset = Button(
                toolbar[1, column];
                label = _icon_label(MI_REFRESH),
                width = BUTTON_SIZE,
                height = BUTTON_SIZE,
                buttoncolor = BUTTON_BACKGROUND
            )
            column += 1
            widgets[:reset] = reset
            on(reset.clicks) do _
                foreach(_reset_panel_limits!, panels)
                context.status[] = "Axis limits reset"
            end
        end
        if definitions.export_svg
            save_button = Button(
                toolbar[1, column];
                label = _icon_label(MI_SAVE),
                width = BUTTON_SIZE,
                height = BUTTON_SIZE,
                buttoncolor = BUTTON_BACKGROUND
            )
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
                scale = enabled ? :log10 : :linear
                foreach(
                    panel -> _set_axis_scale!(
                        panel.axis,
                        panel.view.xaxis,
                        :x,
                        get(page.kwargs, :x_exponent, 0),
                        scale
                    ),
                    panels
                )
                foreach(_reset_panel_limits!, panels)
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
                scale = enabled ? :log10 : :linear
                foreach(
                    panel -> _set_axis_scale!(
                        panel.axis,
                        panel.view.yaxis,
                        :y,
                        get(page.kwargs, :y_exponent, 0),
                        scale
                    ),
                    panels
                )
                foreach(_reset_panel_limits!, panels)
                context.status[] = enabled ?
                                   "y-axis scale set to log" :
                                   "y-axis scale set to linear"
            end
        end
        definitions.legend && legend !== nothing && (widgets[:legend] = legend)
        Label(figure[status_row, 1:2], context.status; halign = :left, fontsize = 11)
        rowsize!(figure.layout, 1, Fixed(TOOLBAR_HEIGHT))
        rowsize!(figure.layout, canvas_row, Relative(1))
        rowsize!(figure.layout, status_row, Fixed(STATUSBAR_HEIGHT))
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
    active = BackendHandler.ensure_backend!(backend)
    built = UIPlot[]
    with_theme(_theme(; export_mode)) do
        for page in render_spec.figures
            context = _context(active, display, page.title)
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
    message = "Saved SVG to $output"
    plot.context.status[] = message
    @info message
    return output
end

export_svg(args...; kwargs...) = PlotBuilder.export_svg(args...; kwargs...)

end # module UIComponents
