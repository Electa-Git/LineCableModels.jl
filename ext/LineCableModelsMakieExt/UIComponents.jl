"""
    LineCableModelsMakieExt.UIComponents

Draw detached PlotBuilder pages through the standard Makie shell.
"""
module UIComponents

using Makie
using Dates
using Printf: @sprintf

import LineCableModels
import LineCableModels.PlotBuilder
import LineCableModels.PlotBuilder: validate_export_theme
using LineCableModels: nominal, standard_uncertainty
using LineCableModels.Units: label
using LineCableModels.PlotBuilder:
                                   ColorbarDefinition, ExportDefinition,
                                   LegendDefinition, PlotPage, PlotRecipe, UIPlot

export build, export_svg

include("context.jl")

const COLORBAR_WIDTH = 140
const COLORBAR_TICK_LABEL_SIZE = 12
const COLORBAR_LABEL_SIZE = 14
const COLORBAR_LABEL_GAP = 8
const COLORBAR_ROW_GAP = 8
const LEGEND_DOCK_WIDTH = 220
const LEGEND_HEIGHT_TOLERANCE = 1
const LEGEND_MARGIN = (3.0f0, 3.0f0, 3.0f0, 3.0f0)
const OUTER_WINDOW_INSET = 3.0
const GRID_ROW_GAP = 6
const GRID_COLUMN_GAP = 6
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
const EXPORT_FALLBACK_DIRECTORY = "linecablemodels-exports"

mutable struct ResponsiveLegend
    legend::Any
    title::Any
    entries::Any
    ellipsis_entry::Any
    capacity::Int
    heights::Dict{Int, Float64}
    fitting::Bool
end

function _theme(; export_mode::Bool = false, export_theme::Symbol = :default)
    validate_export_theme(export_theme)
    background = export_mode ? BACKGROUND_EXPORT : BACKGROUND_INTERACTIVE
    base = export_mode && export_theme === :publication ? Makie.theme_latexfonts() : Theme()
    custom = Theme(
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
        Legend = (; fontsize = 14, labelsize = 14, margin = LEGEND_MARGIN),
        Colorbar = (; labelsize = 14, ticklabelsize = 14)
    )
    return merge(base, custom)
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

function _axis_label(spec, exponent::Int, scale::Symbol)
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
        axis, spec, dim::Symbol, exponent::Int, scale::Symbol)
    spec === nothing && throw(ArgumentError("cannot set an absent axis scale"))
    scale in spec.allowed_scales || throw(
        ArgumentError("axis :$dim does not allow scale :$scale"),
    )
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

function _series_group(series, index::Int)
    return series.group === nothing ? Symbol("series_$index") : series.group
end

function _series_visible(panel::UIPanel, series, index::Int)
    group = _series_group(series, index)
    return all(plot_object -> plot_object.visible[], panel.groups[group])
end

function _axis_values(panel::UIPanel, dim::Symbol; include_uncertainty::Bool = false)
    values = Float64[]
    for (index, series) in enumerate(panel.metadata.series)
        _series_visible(panel, series, index) || continue
        data = dim === :x ? series.xdata : series.ydata
        data === nothing && continue
        for sample in data
            nominal_value = nominal(sample)
            nominal_value isa Real || continue
            numeric = Float64(nominal_value)
            isfinite(numeric) || continue
            uncertainty = abs(Float64(standard_uncertainty(sample)))
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
    metadata = panel.metadata
    autolimits!(axis)
    all(isempty(_axis_values(panel, dim)) for dim in (:x, :y)) && return axis
    metadata.limits !== nothing && return axis
    for dim in (:x, :y)
        values = _axis_values(panel, dim)
        isempty(values) && continue
        scale = dim === :x ? axis.xscale[] : axis.yscale[]
        _nearly_constant(values) || continue
        interval_values = _axis_values(panel, dim; include_uncertainty = true)
        limits = scale === Makie.log10 ?
                 _log_decade_limits(interval_values) :
                 _linear_constant_limits(values, interval_values)
        dim === :x ? xlims!(axis, limits...) : ylims!(axis, limits...)
    end
    return axis
end

function _observe_visibility_limits!(panels, context::UIContext)
    for panel in panels
        panel.metadata.limits === nothing || continue
        panel.metadata.aspect === :data && continue
        for plots in values(panel.groups), plot_object in plots

            observer = on(plot_object.visible) do _
                _reset_panel_limits!(panel)
                context.status[] = "Axis limits fitted to visible series"
                return nothing
            end
            push!(context.observers, observer)
        end
    end
    return panels
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
    nominal_values = nominal.(values)
    errors = standard_uncertainty.(values)
    return nominal_values, any(error -> !iszero(error), errors) ? errors : nothing
end

function _line_errors!(plots, axis, x, y, xerror, yerror, visible)
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
                linewidth = 1,
                visible
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
                linewidth = 1,
                visible
            )
        )
    end
    return plots
end

function _draw_line!(
        axis,
        xdata,
        ydata;
        label = nothing,
        visible = true,
        attributes::NamedTuple = (;)
)
    plots = Any[]
    x, xerror = _numeric_values(xdata)
    y, yerror = _numeric_values(ydata)
    push!(plots, lines!(
        axis, x, y; label, visible, attributes...))
    return _line_errors!(plots, axis, x, y, xerror, yerror, visible)
end

function _sanitize_filename(value::AbstractString)
    sanitized = lowercase(strip(value))
    sanitized = replace(sanitized, r"[^0-9a-z]+" => "_")
    sanitized = strip(sanitized, '_')
    return isempty(sanitized) ? "linecablemodels_plot" : sanitized
end

function _normalized_path_parts(path::AbstractString)
    parts = collect(splitpath(normpath(realpath(path))))
    return Sys.iswindows() ? lowercase.(parts) : parts
end

function _path_within(path::AbstractString, root::AbstractString)
    path_parts = _normalized_path_parts(path)
    root_parts = _normalized_path_parts(root)
    length(path_parts) >= length(root_parts) || return false
    return path_parts[1:length(root_parts)] == root_parts
end

function _export_directory()
    current = abspath(pwd())
    package = abspath(pkgdir(PlotBuilder))
    _path_within(current, package) || return current
    fallback = joinpath(tempdir(), EXPORT_FALLBACK_DIRECTORY)
    mkpath(fallback)
    return fallback
end

function _available_path(page::PlotPage)
    base = _sanitize_filename(_page_export(page).name)
    stamp = Dates.format(Dates.now(), EXPORT_TIMESTAMP_FORMAT)
    directory = _export_directory()
    candidate = joinpath(directory, "$(base)_$(stamp).svg")
    index = 2
    while ispath(candidate)
        candidate = joinpath(directory, "$(base)_$(stamp)_$(index).svg")
        index += 1
    end
    return candidate
end

function _open_command(path::AbstractString)
    if Sys.iswindows()
        return Cmd(["cmd", "/c", "start", "", path])
    elseif Sys.isapple()
        executable = Sys.which("open")
        return executable === nothing ? nothing : `$executable $path`
    end
    executable = Sys.which("xdg-open")
    executable !== nothing && return `$executable $path`
    executable = Sys.which("gio")
    return executable === nothing ? nothing : `$executable open $path`
end

function _open_export(path::AbstractString)
    command = _open_command(path)
    command === nothing && return false
    try
        process = run(pipeline(ignorestatus(command); stdout = devnull, stderr = devnull))
        return success(process)
    catch error
        @warn "could not open exported SVG with the system application" path exception = (
            error,
            catch_backtrace()
        )
        return false
    end
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

_makie_alignment(value::Symbol) = value === :stretch ? :center : value

function _colorbars!(
        slot,
        descriptors;
        halign::Symbol = :left,
        valign::Symbol = :top
)
    isempty(descriptors) && return nothing
    grid = GridLayout(
        width = LEGEND_DOCK_WIDTH,
        tellwidth = true,
        halign = _makie_alignment(halign),
        valign = _makie_alignment(valign)
    )
    grid.default_colgap = Fixed(COLORBAR_LABEL_GAP)
    grid.default_rowgap = Fixed(COLORBAR_ROW_GAP)
    slot[] = grid
    colorbars = Any[]
    for (row, descriptor) in enumerate(descriptors)
        Label(
            grid[row, 1],
            descriptor.label;
            halign = :right,
            valign = :center,
            fontsize = COLORBAR_LABEL_SIZE
        )
        push!(colorbars,
            Colorbar(
                grid[row, 2];
                colormap = descriptor.colormap,
                limits = descriptor.limits,
                ticks = _colorbar_ticks(descriptor.ticks),
                label = "",
                labelvisible = false,
                vertical = false,
                height = 14,
                ticklabelsize = COLORBAR_TICK_LABEL_SIZE,
                alignmode = Mixed(left = 0, right = 0),
                tellwidth = false
            ))
    end
    colsize!(grid, 1, Auto(true))
    colsize!(grid, 2, Fixed(COLORBAR_WIDTH))
    return (; grid, colorbars)
end

function _set_legend_capacity!(state::ResponsiveLegend, capacity::Int)
    total = length(state.entries)
    0 <= capacity <= total || throw(BoundsError(state.entries, capacity))
    capacity == state.capacity && return state
    displayed = copy(state.entries[1:capacity])
    capacity < total && push!(displayed, state.ellipsis_entry)
    state.legend.entrygroups[] = [(state.title, displayed)]
    state.capacity = capacity
    return state
end

function _legend_height!(state::ResponsiveLegend, capacity::Int)
    return get!(state.heights, capacity) do
        _set_legend_capacity!(state, capacity)
        height = state.legend.layoutobservables.autosize[][2]
        height === nothing && return 0.0
        return Float64(height)
    end
end

function _fit_legend!(state::ResponsiveLegend, available_height::Real)
    state.fitting && return state
    state.fitting = true
    try
        total = length(state.entries)
        available = max(0.0, Float64(available_height) - LEGEND_HEIGHT_TOLERANCE)
        if _legend_height!(state, total) <= available
            return _set_legend_capacity!(state, total)
        end
        lower = 0
        upper = max(0, total - 1)
        best = 0
        while lower <= upper
            middle = (lower + upper) ÷ 2
            if _legend_height!(state, middle) <= available
                best = middle
                lower = middle + 1
            else
                upper = middle - 1
            end
        end
        return _set_legend_capacity!(state, best)
    finally
        state.fitting = false
    end
end

function _observe_legend!(
        state::ResponsiveLegend,
        slot_grid,
        context::UIContext
)
    observer = on(slot_grid.layoutobservables.computedbbox) do bounding_box
        _fit_legend!(state, bounding_box.widths[2])
        return nothing
    end
    push!(context.observers, observer)
    return state
end

function _legend!(
        slot,
        panels;
        width = nothing,
        overflow::Symbol = :ellipsis
)
    groups, group_labels, group_order = _visibility_groups(panels)
    entries = Any[]
    labels = String[]
    for group in group_order
        haskey(group_labels, group) || continue
        push!(entries, groups[group])
        push!(labels, group_labels[group])
    end
    isempty(entries) && return nothing
    dimensions = width === nothing ? (;) : (; width)
    ellipsis = PolyElement(color = :transparent, strokecolor = :transparent)
    legend_entries = [entry => (; polystrokecolor = :transparent, polystrokewidth = 0)
                      for entry in entries]
    contents = overflow === :ellipsis ? Any[legend_entries..., ellipsis] : legend_entries
    legend_labels = overflow === :ellipsis ? [labels; "(...)"] : labels
    legend = Legend(
        slot,
        contents,
        legend_labels;
        dimensions...,
        tellheight = overflow === :show_all,
        halign = :left,
        valign = :top
    )
    overflow === :show_all && return (; legend, responsive = nothing)
    title, legend_entries = only(legend.entrygroups[])
    complete_entries = copy(legend_entries[1:(end - 1)])
    responsive = ResponsiveLegend(
        legend,
        title,
        complete_entries,
        last(legend_entries),
        -1,
        Dict{Int, Float64}(),
        false
    )
    _set_legend_capacity!(responsive, length(complete_entries))
    return (; legend, responsive)
end

function _window_padding(padding::NTuple{4, <:Real})
    return ntuple(index -> max(OUTER_WINDOW_INSET, padding[index]), 4)
end

function _page_supports_log(panels, dim::Symbol)
    isempty(panels) && return false
    return all(panels) do panel
        specification = dim === :x ? panel.metadata.xaxis : panel.metadata.yaxis
        specification !== nothing && :log10 in specification.allowed_scales
    end
end

include("shell.jl")

function _current_scale(scale)
    scale === Makie.log10 && return :log10
    scale === Makie.identity && return :linear
    throw(ArgumentError("SVG export supports linear and log10 axis scales"))
end

function _current_limits(axis)
    limits = axis.finallimits[]
    xlimits = (limits.origin[1], limits.origin[1] + limits.widths[1])
    ylimits = (limits.origin[2], limits.origin[2] + limits.widths[2])
    return xlimits, ylimits
end

_replay_page(plot::UIPlot, ::Type{<:PlotBuilder.AbstractPlotDefinition}) = plot.page

function _current_recipe(plot::UIPlot)
    page = _replay_page(plot, plot.render.definition)
    return PlotRecipe(plot.render.definition, (page,))
end

function _block_vertical_bounds(block)
    layout = block.layoutobservables
    bounding_box = layout.computedbbox[]
    protrusions = layout.protrusions[]
    bottom = bounding_box.origin[2] - protrusions.bottom
    top = bounding_box.origin[2] + bounding_box.widths[2] + protrusions.top
    return Float64(bottom), Float64(top)
end

function _export_dock_growth(figure, page::PlotPage)
    legends = filter(block -> block isa Legend, figure.content)
    isempty(legends) && return 0.0
    legend_bottom = minimum(first(_block_vertical_bounds(legend)) for legend in legends)
    all_colorbars = filter(block -> block isa Colorbar, figure.content)
    length(all_colorbars) == length(_page_colorbars(page)) || error(
        "rendered colorbars no longer match the declarative page",
    )
    required_bottom = 0.0
    if !isempty(all_colorbars)
        scale_top = mapreduce(
            block -> last(_block_vertical_bounds(block)),
            max,
            all_colorbars
        )
        required_bottom = scale_top + COLORBAR_ROW_GAP
    end
    return max(0.0, required_bottom - legend_bottom)
end

function _fit_export_content!(figure, page::PlotPage)
    fitted_size = Makie.resize_to_layout!(figure)
    target_size = Tuple(max.(page.size, ceil.(Int, fitted_size)))
    Makie.resize!(figure, target_size...)
    for _ in 1:4
        Makie.update_state_before_display!(figure)
        growth = _export_dock_growth(figure, page)
        growth <= LEGEND_HEIGHT_TOLERANCE && break
        target_size = (target_size[1], target_size[2] + ceil(Int, growth))
        Makie.resize!(figure, target_size...)
    end
    Makie.update_state_before_display!(figure)
    return target_size
end

function PlotBuilder.export_svg(
        plot::UIPlot;
        path::Union{Nothing, AbstractString} = nothing,
        theme::Union{Nothing, Symbol} = nothing,
        open_file::Union{Nothing, Bool} = nothing
)
    PlotBuilder.backend_available(:cairo) || throw(
        ArgumentError(
        "SVG export requires CairoMakie; load CairoMakie first with `using CairoMakie`",
    ),
    )
    output = path === nothing ? _available_path(plot.page) : abspath(String(path))
    definition = _page_export(plot.page)
    export_theme = theme === nothing ? definition.theme : theme
    should_open = open_file === nothing ? definition.open_file : open_file
    lowercase(splitext(output)[2]) == ".svg" || throw(
        ArgumentError("SVG export paths must use the .svg extension"),
    )
    ispath(output) && throw(ArgumentError("refusing to overwrite existing file: $output"))
    plot.context.status[] = "Exporting SVG..."
    PlotBuilder.with_backend(:cairo) do
        one_page = _current_recipe(plot)
        exported = build(
            one_page;
            backend = :cairo,
            display = false,
            controls = false,
            export_mode = true,
            export_theme
        )
        exported_plot = only(exported)
        _fit_export_content!(exported_plot.figure, exported_plot.page)
        Makie.save(output, exported_plot.figure)
    end
    opened = should_open && _open_export(output)
    message = if opened
        "Saved SVG to $output and opened it with the system application"
    elseif should_open
        "Saved SVG to $output; automatic opening was unavailable"
    else
        "Saved SVG to $output"
    end
    plot.context.status[] = message
    @info message
    return output
end

export_svg(args...; kwargs...) = PlotBuilder.export_svg(args...; kwargs...)

include("native.jl")
include("lineparameters.jl")
include("montecarlo.jl")
include("previews.jl")

end # module UIComponents
