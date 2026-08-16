module UIComponents

using Makie
using Dates
using Printf: @sprintf

import LineCableModels.PlotBuilder
using LineCableModels.Utils: to_nominal, uncertainty_value
using LineCableModels.PlotBuilder:
                                   AbstractTrackSize, FixedTrack, RelativeTrack,
                                   ContentTrack, GridArea, GridSpec, SlotSpec,
                                   LayoutSpec, AxisSpec, SeriesSpec, ViewSpec,
                                   PageSpec, RenderSpec, UIPlot
const BackendHandler = PlotBuilder.BackendHandler

export build, export_svg

const COLORBAR_WIDTH = 140
const COLORBAR_TICK_LABEL_SIZE = 12
const COLORBAR_LABEL_SIZE = 14
const COLORBAR_LABEL_GAP = 8
const COLORBAR_ROW_GAP = 8
const LEGEND_DOCK_WIDTH = 220
const LEGEND_HEIGHT_TOLERANCE = 1
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

mutable struct UIContext
    backend::Symbol
    interactive::Bool
    window::Any
    status::Observable{String}
    observers::Vector{Any}
end

struct UIPanel
    view::ViewSpec
    axis::Any
    plots::Vector{Any}
    groups::Dict{Symbol, Vector{Any}}
    group_labels::Dict{Symbol, String}
    group_order::Vector{Symbol}
end

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
    PlotBuilder._validate_export_theme(export_theme)
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
        Legend = (; fontsize = 14, labelsize = 14),
        Colorbar = (; labelsize = 14, ticklabelsize = 14)
    )
    return merge(base, custom)
end

function _context(
        active::Symbol,
        display::Bool,
        title::AbstractString,
        initial_status::AbstractString
)
    interactive = display && active in (:gl, :wgl)
    window = interactive && active === :gl ?
             BackendHandler.make_screen(
        "Fig. $(BackendHandler.next_fignum()) – $title";
        backend = :gl
    ) : nothing
    return UIContext(
        active,
        interactive,
        window,
        Observable(String(initial_status)),
        Any[]
    )
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

function _series_group(series::SeriesSpec, index::Int)
    return series.group === nothing ? Symbol("series_$index") : series.group
end

function _series_visible(panel::UIPanel, series::SeriesSpec, index::Int)
    group = _series_group(series, index)
    return all(plot_object -> plot_object.visible[], panel.groups[group])
end

function _axis_values(panel::UIPanel, dim::Symbol; include_uncertainty::Bool = false)
    values = Float64[]
    for (index, series) in enumerate(panel.view.series)
        _series_visible(panel, series, index) || continue
        data = dim === :x ? series.xdata : series.ydata
        data === nothing && continue
        for sample in data
            nominal = to_nominal(sample)
            nominal isa Real || continue
            numeric = Float64(nominal)
            isfinite(numeric) || continue
            uncertainty = abs(Float64(uncertainty_value(sample)))
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
    all(isempty(_axis_values(panel, dim)) for dim in (:x, :y)) && return axis
    autolimits!(axis)
    view.limits !== nothing && return axis
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
        panel.view.limits === nothing || continue
        panel.view.aspect === :data && continue
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
    nominal = to_nominal.(values)
    errors = uncertainty_value.(values)
    return nominal, any(error -> !iszero(error), errors) ? errors : nothing
end

function _line_errors!(plots, axis, series::SeriesSpec, x, y, xerror, yerror)
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
                visible = series.visible
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
                visible = series.visible
            )
        )
    end
    return plots
end

function draw!(axis, ::Val{:line}, series::SeriesSpec)
    plots = Any[]
    x, xerror = _numeric_values(series.xdata)
    y, yerror = _numeric_values(series.ydata)
    push!(plots, lines!(
        axis, x, y; label = series.label, visible = series.visible, series.attributes...))
    return _line_errors!(plots, axis, series, x, y, xerror, yerror)
end

function draw!(axis, ::Val{:scatter}, series::SeriesSpec)
    x, _ = _numeric_values(series.xdata)
    y, _ = _numeric_values(series.ydata)
    return Any[scatter!(axis, x, y;
        label = series.label, visible = series.visible, series.attributes...)]
end

function draw!(axis, ::Val{:histogram}, series::SeriesSpec)
    values, _ = _numeric_values(series.xdata)
    return Any[hist!(axis, values;
        label = series.label, visible = series.visible, series.attributes...)]
end

function draw!(axis, ::Val{:stairs}, series::SeriesSpec)
    x, _ = _numeric_values(series.xdata)
    y, _ = _numeric_values(series.ydata)
    return Any[stairs!(axis, x, y;
        label = series.label, visible = series.visible, series.attributes...)]
end

function draw!(axis, ::Val{:heatmap}, series::SeriesSpec)
    return Any[heatmap!(axis, series.xdata, series.ydata, series.zdata;
        visible = series.visible, series.attributes...)]
end

function draw!(axis, ::Val{:polygon}, series::SeriesSpec)
    return Any[poly!(axis, series.zdata;
        label = series.label, visible = series.visible, series.attributes...)]
end

function draw!(axis, ::Val{:hline}, series::SeriesSpec)
    return Any[hlines!(axis, series.ydata;
        label = series.label, visible = series.visible, series.attributes...)]
end

function draw!(axis, ::Val{kind}, series::SeriesSpec) where {kind}
    throw(ArgumentError("unsupported PlotBuilder primitive :$kind"))
end

function _axis_attributes(view::ViewSpec)
    axes = (view.xaxis, view.yaxis, view.zaxis)
    axis_attributes = foldl(
        merge,
        (axis.attributes for axis in axes if axis !== nothing);
        init = (;)
    )
    return merge(axis_attributes, view.attributes)
end

function _axis(parent, view::ViewSpec, page::PageSpec)
    xaxis = view.xaxis
    yaxis = view.yaxis
    x_exponent = xaxis === nothing ? 0 : xaxis.exponent
    y_exponent = yaxis === nothing ? 0 : yaxis.exponent
    xscale = xaxis === nothing ? :linear : xaxis.scale
    yscale = yaxis === nothing ? :linear : yaxis.scale
    attributes = _axis_attributes(view)
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
        aspect = view.aspect === :data ? DataAspect() : view.aspect,
        attributes...
    )
    plots = Any[]
    groups = Dict{Symbol, Vector{Any}}()
    group_labels = Dict{Symbol, String}()
    group_order = Symbol[]
    for (index, series) in enumerate(view.series)
        drawn = draw!(axis, Val(series.kind), series)
        append!(plots, drawn)
        group = series.group === nothing ? Symbol("series_$index") : series.group
        haskey(groups, group) || push!(group_order, group)
        append!(get!(groups, group, Any[]), drawn)
        if series.label !== nothing && !isempty(series.label)
            group_labels[group] = series.label
        end
    end
    if view.limits !== nothing
        xlimits, ylimits = view.limits
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

function _available_path(page::PageSpec)
    base = _sanitize_filename(page.export_spec.name)
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

function _colorbars!(slot, descriptors, specification::SlotSpec)
    isempty(descriptors) && return nothing
    grid = GridLayout(
        width = LEGEND_DOCK_WIDTH,
        tellwidth = true,
        halign = _makie_alignment(specification.halign),
        valign = _makie_alignment(specification.valign)
    )
    grid.default_colgap = Fixed(COLORBAR_LABEL_GAP)
    grid.default_rowgap = Fixed(COLORBAR_ROW_GAP)
    slot[] = grid
    for (row, descriptor) in enumerate(descriptors)
        Label(
            grid[row, 1],
            descriptor.label;
            halign = :right,
            valign = :center,
            fontsize = COLORBAR_LABEL_SIZE
        )
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
        )
    end
    colsize!(grid, 1, Auto(true))
    colsize!(grid, 2, Fixed(COLORBAR_WIDTH))
    return grid
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
        panels,
        specification::SlotSpec;
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
        halign = _makie_alignment(specification.halign),
        valign = _makie_alignment(specification.valign)
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

function _shares_side_dock(page::PageSpec, colorbar_slot_name::Symbol)
    page.legend.enabled || return false
    legend_slot = only(slot for slot in page.layout.slots if slot.name === page.legend.slot)
    colorbar_slot = only(
        slot for slot in page.layout.slots if slot.name === colorbar_slot_name)
    return legend_slot.parent === colorbar_slot.parent &&
           legend_slot.area.columns == colorbar_slot.area.columns &&
           last(legend_slot.area.rows) < first(colorbar_slot.area.rows)
end

function _legend_dock_width(page::PageSpec)
    any(colorbar -> _shares_side_dock(page, colorbar.slot), page.colorbars) ||
        return nothing
    return LEGEND_DOCK_WIDTH
end

_track_size(track::FixedTrack) = Fixed(track.value)
_track_size(track::RelativeTrack) = Relative(track.weight)
_track_size(::ContentTrack) = Auto(true)

function _apply_grid_spec!(grid, specification::GridSpec)
    grid.default_rowgap = Fixed(specification.rowgap)
    grid.default_colgap = Fixed(specification.columngap)
    isempty(grid.addedrowgaps) || rowgap!(grid, Fixed(specification.rowgap))
    isempty(grid.addedcolgaps) || colgap!(grid, Fixed(specification.columngap))
    for (index, track) in enumerate(specification.rows)
        index <= length(grid.rowsizes) && rowsize!(grid, index, _track_size(track))
    end
    for (index, track) in enumerate(specification.columns)
        index <= length(grid.colsizes) && colsize!(grid, index, _track_size(track))
    end
    return grid
end

function _grid_position(parent, area::GridArea)
    return parent[area.rows, area.columns]
end

function _materialize_layout(figure, specification::LayoutSpec)
    PlotBuilder.validate(specification)
    root_specification = only(filter(grid -> grid.parent === nothing, specification.grids))
    grids = Dict{Symbol, Any}(root_specification.name => figure.layout)

    pending = [grid for grid in specification.grids if grid.parent !== nothing]
    while !isempty(pending)
        progressed = false
        for grid_specification in copy(pending)
            haskey(grids, grid_specification.parent) || continue
            padding = grid_specification.padding
            grid = GridLayout(
                length(grid_specification.rows),
                length(grid_specification.columns);
                alignmode = all(iszero, padding) ? Inside() : Outside(padding...)
            )
            _grid_position(grids[grid_specification.parent], grid_specification.area)[] = grid
            grids[grid_specification.name] = grid
            deleteat!(pending, findfirst(==(grid_specification), pending))
            progressed = true
        end
        progressed || error("validated layout could not be materialized")
    end

    slot_specs = Dict{Symbol, SlotSpec}()
    for slot_specification in specification.slots
        slot_specs[slot_specification.name] = slot_specification
    end
    for grid_specification in specification.grids
        _apply_grid_spec!(grids[grid_specification.name], grid_specification)
    end
    return (; grids, slot_specs, collapsed = Set{Tuple{Symbol, Int}}())
end

function _collapse_slot!(layout::LayoutSpec, materialized, name::Symbol)
    index = findfirst(slot -> slot.name === name, layout.slots)
    index === nothing && return nothing
    slot = layout.slots[index]
    parent = materialized.grids[slot.parent]
    sibling_areas = GridArea[]
    append!(
        sibling_areas,
        [grid.area
         for grid in layout.grids
         if grid.parent === slot.parent && grid.area !== nothing]
    )
    append!(
        sibling_areas,
        [other.area
         for other in layout.slots
         if other.parent === slot.parent && other.name !== slot.name]
    )
    rows = filter(slot.area.rows) do row
        all(area -> row ∉ area.rows, sibling_areas)
    end
    foreach(row -> push!(materialized.collapsed, (slot.parent, row)), rows)
    foreach(
        row -> row <= length(parent.rowsizes) && rowsize!(parent, row, Fixed(0)),
        rows
    )
    return nothing
end

function _apply_layout_specs!(layout::LayoutSpec, materialized)
    for specification in layout.grids
        _apply_grid_spec!(materialized.grids[specification.name], specification)
    end
    for (grid_name, row) in materialized.collapsed
        parent = materialized.grids[grid_name]
        row <= length(parent.rowsizes) && rowsize!(parent, row, Fixed(0))
    end
    return nothing
end

function _view_positions(views::AbstractVector{<:ViewSpec})
    isempty(views) && return Tuple{ViewSpec, GridArea}[]
    if all(view -> view.placement.area === nothing, views)
        columns = max(1, ceil(Int, sqrt(length(views))))
        return [(
                    view,
                    GridArea((index - 1) ÷ columns + 1, (index - 1) % columns + 1)
                ) for (index, view) in enumerate(views)]
    end
    return [(view, view.placement.area) for view in views]
end

function _slot_position(materialized, name::Symbol)
    specification = materialized.slot_specs[name]
    return _grid_position(materialized.grids[specification.parent], specification.area)
end

function _slot_grid(
        materialized,
        name::Symbol;
        width = Auto(),
        height = Auto(),
        tellwidth::Bool = true,
        tellheight::Bool = true
)
    specification = materialized.slot_specs[name]
    grid = GridLayout(
        width = width,
        height = height,
        tellwidth = tellwidth,
        tellheight = tellheight,
        halign = _makie_alignment(specification.halign),
        valign = _makie_alignment(specification.valign)
    )
    _slot_position(materialized, name)[] = grid
    return grid
end

function _build_panels(page::PageSpec, materialized)
    panels = UIPanel[]
    for slot_name in unique(view.placement.slot for view in page.views)
        slot = _slot_grid(materialized, slot_name)
        slot.default_rowgap = Fixed(GRID_ROW_GAP)
        slot.default_colgap = Fixed(GRID_COLUMN_GAP)
        views = [view for view in page.views if view.placement.slot === slot_name]
        for (view, area) in _view_positions(views)
            push!(panels, _axis(_grid_position(slot, area), view, page))
        end
    end
    return panels
end

function _page_supports_log(panels, dim::Symbol)
    isempty(panels) && return false
    return all(panels) do panel
        specification = dim === :x ? panel.view.xaxis : panel.view.yaxis
        specification !== nothing && :log10 in specification.allowed_scales
    end
end

function _build_colorbars!(page::PageSpec, materialized)
    for slot_name in unique(colorbar.slot for colorbar in page.colorbars)
        descriptors = [colorbar
                       for colorbar in page.colorbars if colorbar.slot === slot_name]
        _colorbars!(
            _slot_position(materialized, slot_name),
            descriptors,
            materialized.slot_specs[slot_name]
        )
    end
    return nothing
end

function _build_page(
        render_spec::RenderSpec,
        page::PageSpec,
        context::UIContext;
        controls::Bool,
        export_mode::Bool
)
    PlotBuilder.validate(render_spec)
    root = only(filter(grid -> grid.parent === nothing, page.layout.grids))
    figure = Figure(size = page.size, figure_padding = root.padding)
    materialized = _materialize_layout(figure, page.layout)
    panels = _build_panels(page, materialized)
    legend = nothing
    responsive_legend = nothing
    legend_slot_grid = nothing
    if page.legend.enabled
        overflow = export_mode ? :show_all : page.legend.overflow
        dock_width = _legend_dock_width(page)
        responsive = overflow === :ellipsis
        legend_slot_grid = _slot_grid(
            materialized,
            page.legend.slot;
            width = dock_width === nothing ? Auto() : dock_width,
            height = responsive ? nothing : Auto(),
            tellheight = !responsive
        )
        built_legend = _legend!(
            legend_slot_grid[1, 1],
            panels,
            materialized.slot_specs[page.legend.slot];
            width = dock_width,
            overflow
        )
        if built_legend === nothing
            _collapse_slot!(page.layout, materialized, page.legend.slot)
        else
            legend = built_legend.legend
            responsive_legend = built_legend.responsive
        end
    else
        _collapse_slot!(page.layout, materialized, page.legend.slot)
    end
    if isempty(page.colorbars)
        for slot in page.layout.slots
            slot.name === :colorbars &&
                _collapse_slot!(page.layout, materialized, slot.name)
        end
    else
        _build_colorbars!(page, materialized)
    end

    widgets = Dict{Symbol, Any}()
    plot_reference = Ref{Any}(nothing)
    if controls
        definitions = page.controls
        xlog_available = _page_supports_log(panels, :x)
        ylog_available = _page_supports_log(panels, :y)
        toolbar_enabled = definitions.reset || definitions.export_svg ||
                          xlog_available || ylog_available
        if toolbar_enabled
            toolbar = _slot_grid(materialized, definitions.slot)
            toolbar.default_colgap = Fixed(4)
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
            if xlog_available
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
                            panel.view.xaxis.exponent,
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
            if ylog_available
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
                            panel.view.yaxis.exponent,
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
        else
            _collapse_slot!(page.layout, materialized, page.controls.slot)
        end
        page.legend.interactive && legend !== nothing && (widgets[:legend] = legend)
        if page.status.enabled
            Label(
                _slot_position(materialized, page.status.slot),
                context.status;
                halign = :left,
                fontsize = 11
            )
        else
            _collapse_slot!(page.layout, materialized, page.status.slot)
        end
    else
        _collapse_slot!(page.layout, materialized, page.controls.slot)
        _collapse_slot!(page.layout, materialized, page.status.slot)
    end

    _apply_layout_specs!(page.layout, materialized)
    if responsive_legend !== nothing
        _observe_legend!(responsive_legend, legend_slot_grid, context)
        Makie.update_state_before_display!(figure)
        bounding_box = legend_slot_grid.layoutobservables.computedbbox[]
        _fit_legend!(responsive_legend, bounding_box.widths[2])
    end
    page.legend.interactive && legend !== nothing &&
        _observe_visibility_limits!(panels, context)
    built = UIPlot(render_spec, page, figure, panels, widgets, context)
    plot_reference[] = built
    return built
end

function build(
        render_spec::RenderSpec;
        backend = nothing,
        display::Bool = true,
        controls::Bool = true,
        export_mode::Bool = false,
        export_theme::Union{Nothing, Symbol} = nothing
)
    PlotBuilder.validate(render_spec)
    active = BackendHandler.ensure_backend!(backend)
    built = UIPlot[]
    for page in render_spec.figures
        page_export_theme = export_theme === nothing ?
                            page.export_spec.theme : export_theme
        with_theme(_theme(; export_mode, export_theme = page_export_theme)) do
            context = _context(active, display, page.title, page.status.initial)
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
    return AxisSpec(
        spec.dim,
        spec.quantity,
        spec.units,
        spec.label,
        _current_scale(scale);
        allowed_scales = spec.allowed_scales,
        exponent = spec.exponent,
        attributes = spec.attributes
    )
end

function _current_limits(axis)
    limits = axis.finallimits[]
    xlimits = (limits.origin[1], limits.origin[1] + limits.widths[1])
    ylimits = (limits.origin[2], limits.origin[2] + limits.widths[2])
    return xlimits, ylimits
end

function _current_series(series::SeriesSpec, panel::UIPanel, index::Int)
    visible = _series_visible(panel, series, index)
    return SeriesSpec(
        series.kind,
        series.xdata,
        series.ydata,
        series.zdata,
        series.label;
        group = series.group,
        visible,
        attributes = series.attributes
    )
end

function _current_view(view::ViewSpec, panel::UIPanel)
    series = [_current_series(item, panel, index)
              for (index, item) in enumerate(view.series)]
    return ViewSpec(
        _axis_with_scale(view.xaxis, panel.axis.xscale[]),
        _axis_with_scale(view.yaxis, panel.axis.yscale[]),
        view.zaxis,
        view.title,
        series,
        view.key;
        placement = view.placement,
        aspect = view.aspect,
        limits = _current_limits(panel.axis),
        attributes = view.attributes
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
        plot.page.title,
        plot.page.size,
        plot.page.key,
        plot.page.layout,
        views;
        controls = plot.page.controls,
        legend = plot.page.legend,
        colorbars = plot.page.colorbars,
        status = plot.page.status,
        export_spec = plot.page.export_spec
    )
end

function _block_vertical_bounds(block)
    layout = block.layoutobservables
    bounding_box = layout.computedbbox[]
    protrusions = layout.protrusions[]
    bottom = bounding_box.origin[2] - protrusions.bottom
    top = bounding_box.origin[2] + bounding_box.widths[2] + protrusions.top
    return Float64(bottom), Float64(top)
end

function _export_dock_growth(figure, page::PageSpec)
    legends = filter(block -> block isa Legend, figure.content)
    isempty(legends) && return 0.0
    legend_bottom = minimum(first(_block_vertical_bounds(legend)) for legend in legends)
    all_colorbars = filter(block -> block isa Colorbar, figure.content)
    rendered_slots = Symbol[]
    for slot_name in unique(colorbar.slot for colorbar in page.colorbars)
        descriptor_count = count(colorbar -> colorbar.slot === slot_name, page.colorbars)
        append!(rendered_slots, fill(slot_name, descriptor_count))
    end
    length(all_colorbars) == length(rendered_slots) || error(
        "rendered colorbars no longer match the declarative page",
    )
    shared_indices = findall(
        slot_name -> _shares_side_dock(page, slot_name),
        rendered_slots
    )
    colorbar_content = all_colorbars[shared_indices]
    required_bottom = 0.0
    if !isempty(colorbar_content)
        scale_top = mapreduce(
            block -> last(_block_vertical_bounds(block)),
            max,
            colorbar_content
        )
        required_bottom = scale_top + COLORBAR_ROW_GAP
    end
    return max(0.0, required_bottom - legend_bottom)
end

function _fit_export_content!(figure, page::PageSpec)
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
    BackendHandler.backend_available(:cairo) || throw(
        ArgumentError(
        "SVG export requires CairoMakie; load CairoMakie first with `using CairoMakie`",
    ),
    )
    output = path === nothing ? _available_path(plot.page) : abspath(String(path))
    export_theme = theme === nothing ? plot.page.export_spec.theme : theme
    should_open = open_file === nothing ? plot.page.export_spec.open_file : open_file
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

end # module UIComponents
