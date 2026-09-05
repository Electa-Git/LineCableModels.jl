const _ADDON_BUTTON_SIZE = 32
const _ADDON_BUTTON_BACKGROUND = Makie.RGBf(0.94, 0.94, 0.94)
const _ADDON_ICON_COLOR = Makie.RGBAf(0.15, 0.15, 0.15, 1.0)
const _ADDON_COLORBAR_DOCK_LENGTH = 140
const _ADDON_MIN_AXIS_CELL_WIDTH = 240
const _ADDON_MIN_AXIS_CELL_HEIGHT = 220
const _ADDON_MIN_WINDOW_SIZE = (600, 320)
const _ADDON_REFRESH_ICON = "\uE5D5"
const _ADDON_SAVE_ICON = "\uE161"
const _ADDON_ICON_FONT = joinpath(
    pkgdir(LineCableModels),
    "assets",
    "fonts",
    "material-icons",
    "MaterialIcons-Regular.ttf"
)
const _ADDON_ZERO_TOLERANCE = sqrt(eps(Float64))

function _addon_theme(; export_mode::Bool = false, export_theme::Symbol = :default)
    export_theme in (:default, :publication) || throw(ArgumentError(
        "export_theme must be :default or :publication",
    ))
    base = export_mode && export_theme === :publication ? Makie.theme_latexfonts() : Theme()
    custom = Theme(
        backgroundcolor = export_mode ? :white : :grey90,
        fonts = (; icons = _ADDON_ICON_FONT),
        Axis = (;
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
        Button = (; buttoncolor = _ADDON_BUTTON_BACKGROUND),
        Legend = (; fontsize = 14, labelsize = 14),
        Colorbar = (; labelsize = 14, ticklabelsize = 14)
    )
    return merge(base, custom)
end

function _addon_activate_backend(backend)
    backend === nothing && begin
        Makie.current_backend() isa Module || throw(ArgumentError(
            "No Makie backend is active. Load CairoMakie, GLMakie, or WGLMakie first.",
        ))
        return nothing
    end
    backend isa Symbol || throw(ArgumentError(
        "backend must be :cairo, :gl, :wgl, or nothing",
    ))
    extension_name,
    package = if backend === :cairo
        (:LineCableModelsCairoMakieExt, "CairoMakie")
    elseif backend === :gl
        (:LineCableModelsGLMakieExt, "GLMakie")
    elseif backend === :wgl
        (:LineCableModelsWGLMakieExt, "WGLMakie")
    else
        throw(ArgumentError("Unknown backend :$backend. Use :cairo, :gl, or :wgl."))
    end
    extension = Base.get_extension(LineCableModels, extension_name)
    extension === nothing && throw(ArgumentError(
        "Backend :$backend is not loaded. Run `using $package` first.",
    ))
    Base.invokelatest(extension.activate!)
    return nothing
end

function _addon_display!(figure, title::AbstractString)
    if current_backend_symbol() === :gl
        extension = Base.get_extension(LineCableModels, :LineCableModelsGLMakieExt)
        extension === nothing &&
            error("GLMakie is active but its LineCableModels extension is unavailable")
        viewport = figure.scene.viewport[]
        minimum_size = Tuple(
            min(Int(viewport.widths[index]), _ADDON_MIN_WINDOW_SIZE[index])
        for index in 1:2
        )
        screen = Base.invokelatest(
            extension.make_screen,
            String(title);
            minimum_size
        )
        Base.display(screen, figure)
    else
        display(figure)
    end
    return figure
end

function _addon_shell(; size, controls::Bool)
    size isa Tuple{Int, Int} && all(>(0), size) || throw(ArgumentError(
        "fig_size must be a tuple of two positive integers",
    ))
    figure = Figure(size = size, figure_padding = (12, 12, 12, 12))
    root = figure.layout
    root.default_rowgap = Fixed(4)
    rowgap!(root, 4)
    body = GridLayout(3, 3; tellwidth = false, tellheight = false)
    body.default_rowgap = Fixed(0)
    body.default_colgap = Fixed(0)
    rowgap!(body, 0)
    colgap!(body, 0)
    canvas = GridLayout(
        ; width = Relative(1), height = Relative(1),
        tellwidth = false, tellheight = false
    )
    canvas.default_rowgap = Fixed(6)
    canvas.default_colgap = Fixed(6)
    rowgap!(canvas, 6)
    colgap!(canvas, 6)
    body[2, 2] = canvas
    toolbar = GridLayout(; halign = :left, valign = :bottom)
    toolbar.default_colgap = Fixed(4)
    colgap!(toolbar, 4)
    status = Observable("Ready.")
    if controls
        root[2, 1] = body
        root[1, 1] = toolbar
        Label(root[3, 1], status; halign = :left, fontsize = 11)
        rowsize!(root, 1, Fixed(32))
        rowsize!(root, 2, Auto(false, 1))
        rowsize!(root, 3, Fixed(16))
    else
        root[1, 1] = body
        rowsize!(root, 1, Auto(false, 1))
    end
    colsize!(root, 1, Auto(false, 1))
    rowsize!(body, 1, Fixed(0))
    rowsize!(body, 2, Auto(false, 1))
    rowsize!(body, 3, Fixed(0))
    colsize!(body, 1, Fixed(0))
    colsize!(body, 2, Auto(false, 1))
    colsize!(body, 3, Fixed(0))
    return (; figure, root, body, canvas, toolbar, status)
end

function _addon_icon(value)
    return Makie.rich(
        value;
        font = :icons,
        fontsize = 18,
        color = _ADDON_ICON_COLOR,
        offset = (0, -0.18)
    )
end

function _addon_button!(toolbar, column::Int, icon)
    return Button(
        toolbar[1, column];
        label = _addon_icon(icon),
        width = _ADDON_BUTTON_SIZE,
        height = _ADDON_BUTTON_SIZE,
        buttoncolor = _ADDON_BUTTON_BACKGROUND
    )
end

_addon_scale(symbol::Symbol) = symbol === :log10 ? Makie.log10 : Makie.identity

function _addon_scientific_exponent(values)
    magnitudes = Float64[]
    for value in values
        numeric = LineCableModels.nominal(value)
        numeric isa Real || continue
        converted = Float64(numeric)
        isfinite(converted) && !iszero(converted) && push!(magnitudes, abs(converted))
    end
    isempty(magnitudes) && return nothing
    return floor(Int, log10(maximum(magnitudes)))
end

function _addon_linear_tickformat(exponent::Int)
    scale = 10.0^exponent
    return values -> [@sprintf("%.4g", value / scale) for value in values]
end

function _addon_decade_ticks(vmin, vmax)
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

function _addon_axis_label(label, exponent::Int, scale::Symbol)
    scale === :log10 && return label
    iszero(exponent) && return label
    formatted = replace(string(exponent), "-" => "−")
    return Makie.rich(
        label,
        "  × 10",
        Makie.superscript(formatted; offset = Makie.Vec2f(0.1, 0.0))
    )
end

function _addon_set_axis!(axis, dim::Symbol, label, allowed, exponent::Int, scale::Symbol)
    scale in allowed || throw(ArgumentError("axis :$dim does not allow scale :$scale"))
    ticks = scale === :log10 ? _addon_decade_ticks : Makie.automatic
    tickformat = scale === :log10 ? Makie.automatic : _addon_linear_tickformat(exponent)
    axis_label = _addon_axis_label(label, exponent, scale)
    if dim === :x
        axis.xticks[] = ticks
        axis.xtickformat[] = tickformat
        axis.xlabel[] = axis_label
        axis.xscale[] = _addon_scale(scale)
    elseif dim === :y
        axis.yticks[] = ticks
        axis.ytickformat[] = tickformat
        axis.ylabel[] = axis_label
        axis.yscale[] = _addon_scale(scale)
    else
        throw(ArgumentError("axis dimension must be :x or :y"))
    end
    return axis
end

function _addon_numeric_values(values)
    nominal_values = LineCableModels.nominal.(values)
    errors = LineCableModels.uncertainty.(values)
    return nominal_values, any(error -> !iszero(error), errors) ? errors : nothing
end

function _addon_line!(axis, xdata, ydata; label, color = nothing, visible = true)
    x, xerror = _addon_numeric_values(xdata)
    y, yerror = _addon_numeric_values(ydata)
    attributes = color === nothing ? (; linewidth = 2) : (; linewidth = 2, color)
    plots = Any[lines!(axis, x, y; label, visible, attributes...)]
    error_color = color === nothing ? :black : color
    yerror === nothing || push!(plots,
        errorbars!(
            axis,
            x,
            y,
            yerror;
            color = error_color,
            direction = :y,
            whiskerwidth = 3,
            linewidth = 1,
            visible
        ))
    xerror === nothing || push!(plots,
        errorbars!(
            axis,
            x,
            y,
            xerror;
            color = error_color,
            direction = :x,
            whiskerwidth = 3,
            linewidth = 1,
            visible
        ))
    return plots
end

function _addon_visible_values(series, dim::Symbol; include_uncertainty::Bool = false)
    values = Float64[]
    for item in series
        all(plot -> plot.visible[], item.plots) || continue
        data = dim === :x ? item.xdata : item.ydata
        data === nothing && continue
        for sample in data
            nominal_value = LineCableModels.nominal(sample)
            nominal_value isa Real || continue
            numeric = Float64(nominal_value)
            isfinite(numeric) || continue
            interval = abs(Float64(LineCableModels.uncertainty(sample)))
            if include_uncertainty && isfinite(interval) && !iszero(interval)
                push!(values, numeric - interval, numeric + interval)
            else
                push!(values, numeric)
            end
        end
    end
    return values
end

function LineCableModels.plotwindow(
        callback::F;
        title::AbstractString,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        size::Tuple{Int, Int} = (800, 400),
        layout = nothing,
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        export_name::AbstractString = title
) where {F}
    _addon_activate_backend(backend)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls)
        dimensions = layout === nothing ? nothing :
                     _native_preview_layout(1, layout)
        callback(shell.canvas)
        if dimensions !== nothing
            rows, columns = dimensions
            for row in 1:rows
                rowsize!(shell.canvas, row, Relative(1 / rows))
            end
            for column in 1:columns
                colsize!(shell.canvas, column, Relative(1 / columns))
            end
        end
        axes = Any[content for content in shell.figure.content if content isa Axis]
        resets = Function[(() -> autolimits!(axis)) for axis in axes]
        foreach(callback -> callback(), resets)
        _addon_finish!(
            shell,
            axes,
            resets,
            Function[],
            Function[],
            Dict{Symbol, Vector{Any}}(),
            Symbol[],
            Dict{Symbol, String}();
            title = String(title),
            figure_title,
            title_attributes,
            legend_position = nothing,
            legend_attributes = (;),
            controls,
            display_plot,
            export_name = String(export_name),
            export_theme,
            open_export
        )
    end
end

function _addon_statistical_plot(
        draw::F,
        xobservation,
        yobservation;
        title,
        figure_title = nothing,
        title_attributes = (;),
        panel_titles = nothing,
        fig_size,
        backend,
        display_plot,
        controls,
        export_theme,
        open_export,
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_labels = nothing,
        legend_attributes = (;),
        legend_overflow = :ellipsis,
        panel_legends = (),
        xlabel = nothing,
        ylabel = nothing
) where {F}
    _addon_activate_backend(backend)
    resolved_panel_titles = _addon_panel_titles(panel_titles, 1)
    panel_title = resolved_panel_titles === nothing ? title :
                  only(resolved_panel_titles)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size = fig_size, controls)
        panel = _addon_panel!(shell, (1, 1))
        axis,
        automatic_labels = _addon_axis!(
            panel.content,
            xobservation,
            yobservation;
            title = panel_title,
            xscale = :linear,
            yscale = :linear,
            xscales = (:linear,),
            yscales = (:linear,)
        )
        axis_labels = (;
            x = xlabel === nothing ? automatic_labels.x : String(xlabel),
            y = ylabel === nothing ? automatic_labels.y : String(ylabel)
        )
        xlabel === nothing || (axis.xlabel[] = axis_labels.x)
        ylabel === nothing || (axis.ylabel[] = axis_labels.y)
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, String}()
        series = NamedTuple[]
        draw(axis, groups, order, labels, series)
        _addon_relabel_legend!(labels, groups, order, legend_labels)
        reset! = () -> _addon_reset!(axis, series, axis_labels)
        reset!()
        _addon_finish!(
            shell,
            Any[axis],
            Function[reset!],
            Function[],
            Function[],
            groups,
            order,
            labels;
            title,
            figure_title,
            title_attributes,
            legend_position,
            legend_anchor,
            legend_title,
            legend_attributes,
            legend_overflow,
            panels = (panel,),
            panel_legends,
            controls,
            display_plot,
            export_name = title,
            export_theme,
            open_export
        )
    end
end

function _addon_nearly_constant(values)
    isempty(values) && return false
    lower, upper = extrema(values)
    scale = max(abs(lower), abs(upper), floatmin(Float64))
    return upper - lower <= max(_ADDON_ZERO_TOLERANCE, 64eps(Float64) * scale)
end

function _addon_constant_limits(values, interval_values, logarithmic::Bool)
    if logarithmic
        all(>(0), interval_values) || throw(DomainError(
            interval_values,
            "logarithmic axes require strictly positive data"
        ))
        lower, upper = extrema(interval_values)
        first_exponent = floor(Int, log10(lower))
        last_exponent = ceil(Int, log10(upper))
        first_exponent == last_exponent && (first_exponent -= 1; last_exponent += 1)
        return 10.0^first_exponent, 10.0^last_exponent
    end
    all(value -> abs(value) <= _ADDON_ZERO_TOLERANCE, interval_values) &&
        return (-_ADDON_ZERO_TOLERANCE, _ADDON_ZERO_TOLERANCE)
    center = sum(extrema(values)) / 2
    base_halfspan = iszero(center) ? 1.0 : 0.05abs(center)
    interval_halfspan = maximum(abs(value - center) for value in interval_values)
    halfspan = max(base_halfspan, 2interval_halfspan)
    return center - halfspan, center + halfspan
end

function _addon_refresh_format!(axis, series, labels)
    for dim in (:x, :y)
        values = _addon_visible_values(series, dim)
        limits = axis.finallimits[]
        index = dim === :x ? 1 : 2
        lower = limits.origin[index]
        limit_values = (lower, lower + limits.widths[index])
        exponent = something(
            _addon_scientific_exponent(values),
            _addon_scientific_exponent(limit_values),
            0
        )
        scale = dim === :x ? axis.xscale[] : axis.yscale[]
        scale === Makie.log10 && continue
        label = dim === :x ? labels.x : labels.y
        if dim === :x
            axis.xtickformat[] = _addon_linear_tickformat(exponent)
            axis.xlabel[] = _addon_axis_label(label, exponent, :linear)
        else
            axis.ytickformat[] = _addon_linear_tickformat(exponent)
            axis.ylabel[] = _addon_axis_label(label, exponent, :linear)
        end
    end
    return axis
end

function _addon_reset!(axis, series, labels)
    autolimits!(axis)
    for dim in (:x, :y)
        values = _addon_visible_values(series, dim)
        isempty(values) && continue
        _addon_nearly_constant(values) || continue
        interval_values = _addon_visible_values(series, dim; include_uncertainty = true)
        logarithmic = (dim === :x ? axis.xscale[] : axis.yscale[]) === Makie.log10
        limits = _addon_constant_limits(values, interval_values, logarithmic)
        dim === :x ? xlims!(axis, limits...) : ylims!(axis, limits...)
    end
    _addon_refresh_format!(axis, series, labels)
    return axis
end

function _addon_axis!(
        position,
        xobservation,
        yobservation;
        title,
        xscale,
        yscale,
        xscales,
        yscales,
        xlabel = nothing,
        ylabel = nothing,
        attributes = (;)
)
    xaxis_label = xlabel === nothing ?
                  LineCableModels.Units.label(xobservation.quantity, xobservation.unit) :
                  String(xlabel)
    yaxis_label = ylabel === nothing ?
                  LineCableModels.Units.label(yobservation.quantity, yobservation.unit) :
                  String(ylabel)
    xexponent = something(_addon_scientific_exponent(xobservation.values), 0)
    yexponent = something(_addon_scientific_exponent(yobservation.values), 0)
    options = merge(
        (;
            title,
            tellwidth = false,
            tellheight = false,
            xscale = _addon_scale(xscale),
            yscale = _addon_scale(yscale),
            xticks = xscale === :log10 ? _addon_decade_ticks : Makie.automatic,
            yticks = yscale === :log10 ? _addon_decade_ticks : Makie.automatic,
            xtickformat = xscale === :log10 ? Makie.automatic :
                          _addon_linear_tickformat(xexponent),
            ytickformat = yscale === :log10 ? Makie.automatic :
                          _addon_linear_tickformat(yexponent),
            xlabel = _addon_axis_label(xaxis_label, xexponent, xscale),
            ylabel = _addon_axis_label(yaxis_label, yexponent, yscale)
        ),
        attributes)
    axis = Axis(position; options...)
    return axis, (; x = xaxis_label, y = yaxis_label)
end

function _addon_positions(count::Int, layout)
    if layout === nothing
        columns = max(1, ceil(Int, sqrt(count)))
        rows = cld(count, columns)
    else
        layout isa Tuple && length(layout) == 2 &&
        all(value -> value isa Integer && !(value isa Bool) && value > 0, layout) ||
            throw(ArgumentError("layout must be a tuple of two positive integers or nothing"))
        rows, columns = Int.(layout)
        rows * columns >= count || throw(DimensionMismatch(
            "layout provides $(rows * columns) slots for $count plots",
        ))
    end
    positions = Tuple((cld(index, columns), mod1(index, columns)) for index in 1:count)
    return positions, (rows, columns)
end

function _addon_panel_titles(panel_titles, expected::Int)
    panel_titles === nothing && return nothing
    panel_titles isa Tuple || panel_titles isa AbstractVector ||
        throw(ArgumentError(
            "panel_titles must be a tuple, vector, or nothing",
        ))
    length(panel_titles) == expected || throw(DimensionMismatch(
        "panel_titles must contain one entry per logical plot panel",
    ))
    return Tuple(String(value) for value in panel_titles)
end

function _addon_panel!(shell, position::Tuple{Int, Int})
    layout = GridLayout(
        3,
        3;
        width = Relative(1),
        height = Relative(1),
        tellwidth = false,
        tellheight = false
    )
    # Empty docks are layout-neutral; a legend's own padding supplies spacing
    # only after the caller activates one of these tracks.
    layout.default_rowgap = Fixed(0)
    layout.default_colgap = Fixed(0)
    rowgap!(layout, 0)
    colgap!(layout, 0)
    rowsize!(layout, 1, Fixed(0))
    rowsize!(layout, 2, Auto(false, 1))
    rowsize!(layout, 3, Fixed(0))
    colsize!(layout, 1, Fixed(0))
    colsize!(layout, 2, Auto(false, 1))
    colsize!(layout, 3, Fixed(0))
    shell.canvas[position...] = layout
    return (; logical_position = position, layout, content = layout[2, 2])
end

function _addon_responsive_axis_grid!(
        figure,
        grid,
        panels,
        axes,
        dimensions;
        hide_inner_y::Bool = false
)
    isempty(axes) && return axes
    length(panels) == length(axes) || throw(DimensionMismatch(
        "responsive plot panels must align with their axes",
    ))
    maximum_columns = dimensions[2]
    maximum_columns == 1 && return axes
    current_columns = Ref(0)
    updating = Ref(false)

    function reflow!(bounding_box)
        updating[] && return nothing
        width = Float64(bounding_box.widths[1])
        height = Float64(bounding_box.widths[2])
        width_columns = clamp(
            floor(Int, width / _ADDON_MIN_AXIS_CELL_WIDTH),
            1,
            maximum_columns
        )
        available_rows = max(1, floor(Int, height / _ADDON_MIN_AXIS_CELL_HEIGHT))
        height_columns = clamp(cld(length(axes), available_rows), 1, maximum_columns)
        # Pick the arrangement whose cell aspect best follows the available
        # canvas, while respecting the minimum usable cell dimensions.  This
        # keeps a short landscape window in columns and lets a tall window
        # reflow the same logical panels into rows.
        aspect_columns = round(Int, sqrt(length(axes) * width / max(height, 1.0)))
        columns = clamp(aspect_columns, height_columns, width_columns)
        columns == current_columns[] && return nothing
        updating[] = true
        try
            rows = cld(length(axes), columns)
            maximum_rows = cld(length(axes), 1)
            for column in 1:maximum_columns
                colsize!(
                    grid,
                    column,
                    column <= columns ? Auto(false, 1) : Fixed(0)
                )
            end
            for (index, (panel, axis)) in enumerate(zip(panels, axes))
                row = cld(index, columns)
                column = mod1(index, columns)
                grid[row, column] = panel.layout
                bottom_row = row == rows
                axis.xlabelvisible[] = bottom_row
                axis.xticklabelsvisible[] = bottom_row
                axis.xticksvisible[] = bottom_row
                if hide_inner_y
                    left_column = column == 1
                    axis.ylabelvisible[] = left_column
                    axis.yticklabelsvisible[] = left_column
                    axis.yticksvisible[] = left_column
                end
            end
            managed_rows = current_columns[] == 0 ? rows : maximum_rows
            for row in 1:managed_rows
                rowsize!(
                    grid,
                    row,
                    row <= rows ? Auto(false, 1) : Fixed(0)
                )
            end
            current_columns[] = columns
        finally
            updating[] = false
        end
        return nothing
    end

    on(figure.scene, grid.layoutobservables.computedbbox) do bounding_box
        reflow!(bounding_box)
    end
    reflow!(grid.layoutobservables.computedbbox[])
    return axes
end

function _addon_center_aspect_canvas!(shell)
    viewport = shell.figure.scene.viewport[]
    initial_width = max(1.0, Float64(minimum(viewport.widths)))
    shell.canvas.width = initial_width
    shell.canvas.tellwidth = true
    colsize!(shell.body, 2, Auto(true))
    shell.body.width = Auto()
    shell.body.tellwidth = true
    shell.body.halign = :center
    updating = Ref(false)

    function fit!(canvas_box)
        updating[] && return nothing
        suggested = shell.body.layoutobservables.suggestedbbox[]
        body_box = shell.body.layoutobservables.computedbbox[]
        all(isfinite, (
            canvas_box.widths...,
            suggested.widths...,
            body_box.widths...
        )) || return nothing
        noncanvas_width = max(
            0.0,
            Float64(body_box.widths[1]) - Float64(canvas_box.widths[1])
        )
        available_width = max(1.0, Float64(suggested.widths[1]) - noncanvas_width)
        target = min(Float64(canvas_box.widths[2]), available_width)
        current = shell.canvas.width
        current isa Real && isapprox(current, target; atol = 0.5) && return nothing
        updating[] = true
        try
            shell.canvas.width = target
        finally
            updating[] = false
        end
        return nothing
    end

    on(shell.figure.scene, shell.canvas.layoutobservables.computedbbox) do bounding_box
        fit!(bounding_box)
    end
    # A legend or colorbar can be added or moved after construction. Refit the
    # square canvas when the surrounding body changes, even if the canvas box
    # itself has not emitted yet.
    on(shell.figure.scene, shell.body.layoutobservables.computedbbox) do _
        fit!(shell.canvas.layoutobservables.computedbbox[])
    end
    fit!(shell.canvas.layoutobservables.computedbbox[])
    return shell
end

function _addon_legend_slot(body, position)
    position === :right && return body[2, 3], :vertical
    position === :left && return body[2, 1], :vertical
    position === :top && return body[1, 2], :horizontal
    position === :bottom && return body[3, 2], :horizontal
    if position isa Tuple && length(position) == 2 &&
       all(value -> value isa Integer && !(value isa Bool) && value > 0, position)
        position == (2, 2) && throw(ArgumentError(
            "grid position (2, 2) is reserved for the plot canvas",
        ))
        return body[Int(position[1]), Int(position[2])], :vertical
    end
    throw(ArgumentError(
        "legend_position must be :inside, :left, :right, :top, :bottom, a positive grid index tuple, or nothing",
    ))
end

function _addon_dock_indices(position)
    position === :right && return (2, 3)
    position === :left && return (2, 1)
    position === :top && return (1, 2)
    position === :bottom && return (3, 2)
    position isa Tuple && return Int.(position)
    return nothing
end

function _addon_activate_dock_tracks!(body, position)
    indices = _addon_dock_indices(position)
    indices === nothing && return body
    row, column = indices
    if row != 2
        rowsize!(body, row, Auto(true, 0))
        rowgap!(body, row < 2 ? row : row - 1, 8)
    end
    if column != 2
        colsize!(body, column, Auto(true, 0))
        colgap!(body, column < 2 ? column : column - 1, 8)
    end
    return body
end

function _addon_deactivate_dock_tracks!(body, position)
    indices = _addon_dock_indices(position)
    indices === nothing && return body
    row, column = indices
    if row != 2
        rowsize!(body, row, Fixed(0))
        rowgap!(body, row < 2 ? row : row - 1, 0)
    end
    if column != 2
        colsize!(body, column, Fixed(0))
        colgap!(body, column < 2 ? column : column - 1, 0)
    end
    return body
end

function _addon_remove_legend!(legend)
    legend === nothing && return nothing
    content = Makie.GridLayoutBase.gridcontent(legend)
    layout = content === nothing ? nothing : content.parent
    delete!(legend)
    if layout !== nothing && isempty(layout.content)
        layout_content = Makie.GridLayoutBase.gridcontent(layout)
        layout_content === nothing ||
            Makie.GridLayoutBase.remove_from_gridlayout!(layout_content)
    end
    return nothing
end

function _addon_set_legend_capacity!(legend, title, entries, ellipsis, capacity, state)
    total = length(entries)
    0 <= capacity <= total || throw(BoundsError(entries, capacity))
    capacity == state[] && return legend
    displayed = copy(entries[1:capacity])
    capacity < total && push!(displayed, ellipsis)
    legend.entrygroups[] = [(title, displayed)]
    state[] = capacity
    return legend
end

function _addon_responsive_legend!(figure, bounding_box, legend)
    title, built_entries = only(legend.entrygroups[])
    complete_entries = copy(built_entries[1:(end - 1)])
    ellipsis = last(built_entries)
    capacity = Ref(-1)
    fitting = Ref(false)
    extents = Dict{Tuple{Symbol, Int}, Float64}()

    function entry_extent(count::Int, orientation::Symbol)
        return get!(extents, (orientation, count)) do
            _addon_set_legend_capacity!(
                legend, title, complete_entries, ellipsis, count, capacity)
            dimension = orientation === :vertical ? 2 : 1
            value = legend.layoutobservables.autosize[][dimension]
            value === nothing ? 0.0 : Float64(value)
        end
    end
    function fit!(bounding_box)
        fitting[] && return nothing
        fitting[] = true
        try
            orientation = legend.orientation[]
            orientation in (:vertical, :horizontal) || return nothing
            dimension = orientation === :vertical ? 2 : 1
            available = max(0.0, Float64(bounding_box.widths[dimension]) - 2.0)
            total = length(complete_entries)
            if entry_extent(total, orientation) <= available
                _addon_set_legend_capacity!(
                    legend, title, complete_entries, ellipsis, total, capacity)
                return nothing
            end
            lower = 0
            upper = max(0, total - 1)
            best = 0
            while lower <= upper
                middle = (lower + upper) ÷ 2
                if entry_extent(middle, orientation) <= available
                    best = middle
                    lower = middle + 1
                else
                    upper = middle - 1
                end
            end
            _addon_set_legend_capacity!(
                legend, title, complete_entries, ellipsis, best, capacity)
        finally
            fitting[] = false
        end
        return nothing
    end

    _addon_set_legend_capacity!(
        legend, title, complete_entries, ellipsis, length(complete_entries), capacity)
    on(figure.scene, bounding_box) do bounds
        fit!(bounds)
    end
    on(figure.scene, legend.orientation) do _
        empty!(extents)
        fit!(bounding_box[])
    end
    return legend
end

function _addon_inside_aligns(anchor)
    if anchor isa Symbol
        value = String(anchor)
        length(value) == 2 || throw(ArgumentError(
            "inside legend anchors must contain a horizontal and vertical letter",
        ))
        horizontal = Dict('l' => :left, 'c' => :center, 'r' => :right)
        vertical = Dict('b' => :bottom, 'c' => :center, 't' => :top)
        haskey(horizontal, value[1]) || throw(ArgumentError(
            "inside legend anchors must begin with l, c, or r",
        ))
        haskey(vertical, value[2]) || throw(ArgumentError(
            "inside legend anchors must end with b, c, or t",
        ))
        return (; halign = horizontal[value[1]], valign = vertical[value[2]])
    end
    anchor isa Tuple && length(anchor) == 2 || throw(ArgumentError(
        "inside legend anchors must be symbols such as :rt or two-element tuples",
    ))
    return (; halign = anchor[1], valign = anchor[2])
end

function _addon_axes_viewport(axes, fallback)
    isempty(axes) && return fallback
    viewports = Tuple(axis.scene.viewport for axis in axes)
    return lift(viewports...) do bounds...
        left = minimum(bound.origin[1] for bound in bounds)
        bottom = minimum(bound.origin[2] for bound in bounds)
        right = maximum(bound.origin[1] + bound.widths[1] for bound in bounds)
        top = maximum(bound.origin[2] + bound.widths[2] for bound in bounds)
        Rect2f((left, bottom), (right - left, top - bottom))
    end
end

function _addon_relabel_legend!(labels, groups, order, requested)
    requested === nothing && return labels
    displayed = Any[group
                    for group in order if haskey(groups, group) && haskey(labels, group)]
    if requested isa AbstractDict
        for group in displayed
            current = labels[group]
            replacement = if haskey(requested, group)
                requested[group]
            elseif haskey(requested, current)
                requested[current]
            else
                continue
            end
            labels[group] = String(replacement)
        end
        return labels
    end
    requested isa Tuple || requested isa AbstractVector ||
        throw(ArgumentError(
            "legend_labels must be a tuple, vector, dictionary, or nothing",
        ))
    length(requested) == length(displayed) || throw(DimensionMismatch(
        "legend_labels must contain one entry for each displayed legend group",
    ))
    for (group, replacement) in zip(displayed, requested)
        labels[group] = String(replacement)
    end
    return labels
end

function _addon_legend!(
        figure,
        body,
        groups,
        order,
        labels;
        position,
        attributes,
        overflow::Symbol,
        title = nothing,
        anchor = :rt,
        inside_bbox = nothing,
        target = nothing,
        target_orientation = nothing
)
    position === nothing && return nothing
    attributes isa NamedTuple ||
        throw(ArgumentError("legend_attributes must be a NamedTuple"))
    overflow in (:ellipsis, :show_all) || throw(ArgumentError(
        "legend_overflow must be :ellipsis or :show_all",
    ))
    entries = Any[]
    displayed = String[]
    for group in order
        haskey(labels, group) || continue
        push!(entries, groups[group])
        push!(displayed, labels[group])
    end
    isempty(entries) && return nothing
    if position === :inside
        inside_bbox === nothing && throw(ArgumentError(
            "inside legends require a figure or panel plot-area bounding box",
        ))
        options = merge(
            (;
                bbox = inside_bbox,
                orientation = :vertical,
                _addon_inside_aligns(anchor)...,
                margin = (10, 10, 10, 10),
                tellwidth = false,
                tellheight = false
            ),
            attributes
        )
        if overflow === :show_all
            return Legend(
                figure,
                entries,
                displayed,
                title;
                options...
            )
        end
        ellipsis = LineElement(color = :transparent)
        legend = Legend(
            figure,
            Any[entries..., ellipsis],
            [displayed; "(...)"],
            title;
            options...
        )
        return _addon_responsive_legend!(
            figure,
            legend.layoutobservables.computedbbox,
            legend
        )
    end
    slot, resolved_orientation = _addon_legend_slot(body, position)
    target === nothing || (slot = target)
    default_orientation = target_orientation === nothing ?
                          resolved_orientation : target_orientation
    options = merge(
        (;
            orientation = default_orientation,
            halign = default_orientation === :vertical ? :left : :center,
            valign = default_orientation === :vertical ? :top : :center,
            tellwidth = position in (:left, :right) || position isa Tuple,
            tellheight = position in (:top, :bottom) || position isa Tuple
        ),
        attributes)
    legend_grid = GridLayout()
    slot[] = legend_grid
    target === nothing && _addon_activate_dock_tracks!(body, position)
    if overflow === :show_all
        return Legend(legend_grid[1, 1], entries, displayed, title; options...)
    end
    ellipsis = LineElement(color = :transparent)
    legend = Legend(
        legend_grid[1, 1],
        Any[entries..., ellipsis],
        [displayed; "(...)"],
        title;
        options...
    )
    return _addon_responsive_legend!(
        figure,
        legend_grid.layoutobservables.computedbbox,
        legend
    )
end

function _addon_plot_belongs_to_axis(plot, axis)
    return try
        getproperty(plot, :parent) === axis.scene
    catch
        false
    end
end

function _addon_panel_legend_data(
        panels,
        axes,
        groups,
        order,
        labels;
        panel_labels = nothing,
        panel_titles = nothing
)
    length(panels) == length(axes) || throw(DimensionMismatch(
        "plot panels must align with their axes",
    ))
    panel_labels === nothing || length(panel_labels) == length(axes) ||
        throw(
            DimensionMismatch("panel legend labels must align with their axes"),
        )
    panel_titles === nothing || length(panel_titles) == length(axes) ||
        throw(
            DimensionMismatch("panel legend titles must align with their axes"),
        )
    result = Dict{Tuple{Int, Int}, Any}()
    for (index, (panel, axis)) in enumerate(zip(panels, axes))
        scoped = Dict{Any, Vector{Any}}()
        scoped_order = Any[]
        scoped_labels = Dict(panel_labels === nothing ? labels : panel_labels[index])
        for key in order
            plots = Any[plot
                        for plot in groups[key] if _addon_plot_belongs_to_axis(plot, axis)]
            isempty(plots) && continue
            scoped[key] = plots
            push!(scoped_order, key)
        end
        title = Ref{Any}(panel_titles === nothing ? nothing : panel_titles[index])
        result[panel.logical_position] = (;
            panel,
            axis,
            groups = scoped,
            order = scoped_order,
            labels = scoped_labels,
            title
        )
    end
    return result
end

function _addon_panel_legend_pairs(value)
    value === nothing && return Pair[]
    value === () && return Pair[]
    value isa Pair && return Pair[value]
    value isa AbstractDict && return collect(pairs(value))
    value isa Tuple && all(item -> item isa Pair, value) && return collect(value)
    throw(ArgumentError(
        "panel_legends must be a pair, a tuple of pairs, a dictionary, or nothing",
    ))
end

function _addon_positive_panel_position(value)
    value isa Tuple && length(value) == 2 &&
    all(index -> index isa Integer && !(index isa Bool) && index > 0, value) ||
        throw(ArgumentError("panel legend keys must be positive `(row, column)` tuples"))
    return (Int(value[1]), Int(value[2]))
end

function _addon_without_legend_controls(options::NamedTuple)
    names = Tuple(filter(
        name -> name ∉ (:position, :overflow, :title, :anchor, :legend_labels),
        keys(options)
    ))
    return NamedTuple{names}(Tuple(getproperty(options, name) for name in names))
end

function _addon_legend_configuration(value; default_position, default_title = nothing)
    value isa Symbol && return (;
        position = value,
        overflow = :ellipsis,
        title = default_title,
        anchor = :rt,
        legend_labels = nothing,
        attributes = (;)
    )
    value isa NamedTuple || throw(ArgumentError(
        "a legend configuration must be a dock symbol or NamedTuple",
    ))
    position = get(value, :position, default_position)
    overflow = get(value, :overflow, :ellipsis)
    title = get(value, :title, default_title)
    anchor = get(value, :anchor, :rt)
    legend_labels = get(value, :legend_labels, nothing)
    attributes = _addon_without_legend_controls(value)
    return (; position, overflow, title, anchor, legend_labels, attributes)
end

function _addon_panel_legends!(figure, panel_data, requested)
    built = Dict{Tuple{Int, Int}, Any}()
    positions = Dict{Tuple{Int, Int}, Any}()
    for pair in _addon_panel_legend_pairs(requested)
        logical_position = _addon_positive_panel_position(first(pair))
        haskey(panel_data, logical_position) || throw(BoundsError(
            collect(keys(panel_data)), logical_position
        ))
        value = last(pair)
        (value === nothing || value === false) && continue
        data = panel_data[logical_position]
        configuration = _addon_legend_configuration(
            value;
            default_position = :right,
            default_title = data.title[]
        )
        value isa NamedTuple && haskey(value, :title) &&
            (data.title[] = configuration.title)
        _addon_relabel_legend!(
            data.labels,
            data.groups,
            data.order,
            configuration.legend_labels
        )
        legend = _addon_legend!(
            figure,
            data.panel.layout,
            data.groups,
            data.order,
            data.labels;
            position = configuration.position,
            attributes = configuration.attributes,
            overflow = configuration.overflow,
            title = configuration.title,
            anchor = configuration.anchor,
            inside_bbox = data.axis.scene.viewport
        )
        if legend !== nothing
            built[logical_position] = legend
            positions[logical_position] = configuration.position
        end
    end
    return (; legends = built, positions)
end

function _addon_colorbar!(position, scale; attributes)
    scale isa NamedTuple || throw(ArgumentError(
        "a material color scale must be a NamedTuple",
    ))
    all(name -> haskey(scale, name), (:colormap, :limits, :ticks, :label)) ||
        throw(ArgumentError(
            "a material color scale requires colormap, limits, ticks, and label",
        ))
    options = merge(
        (;
            colormap = scale.colormap,
            limits = scale.limits,
            ticks = scale.ticks,
            label = scale.label
        ),
        attributes)
    colorbar = Colorbar(position; options...)
    if !haskey(attributes, :alignmode)
        # Makie reserves space perpendicular to a colorbar, but not for labels
        # extending beyond its endpoints. Include the rendered text extents in
        # the native layout, independently of the bar's position or length.
        labels = colorbar.axis.elements[:ticklabels]
        bounds = Makie.fast_string_boundingboxes_obs(labels)
        onany(colorbar.blockscene, bounds, colorbar.vertical,
            colorbar.ticklabelsvisible; update = true) do boxes, vertical, visible
            dimension = vertical ? 2 : 1
            finite_boxes = filter(box -> isfinite(box.origin[dimension]) &&
                                  isfinite(box.widths[dimension]), boxes)
            before = visible ? ceil(maximum(
                box -> -box.origin[dimension], finite_boxes; init = 0.0
            )) : 0.0
            after = visible ? ceil(maximum(
                box -> box.origin[dimension] + box.widths[dimension],
                finite_boxes; init = 0.0
            )) : 0.0
            colorbar.alignmode[] = vertical ?
                                  Mixed(bottom = before, top = after) :
                                  Mixed(left = before, right = after)
        end
    end
    return colorbar
end

function LineCableModels.materialscale!(position, scheme; kwargs...)
    return _addon_colorbar!(position, scheme; attributes = (; kwargs...))
end

function _addon_colorbars!(
        body,
        scales;
        position,
        attributes,
        target = nothing,
        target_orientation = nothing
)
    isempty(scales) && return (; colorbars = (), layout = nothing)
    position === nothing && return (; colorbars = (), layout = nothing)
    attributes isa NamedTuple || throw(ArgumentError(
        "colorbar_attributes must be a NamedTuple",
    ))
    slot, resolved_orientation = _addon_legend_slot(body, position)
    target === nothing || (slot = target)
    default_orientation = target_orientation === nothing ?
                          resolved_orientation : target_orientation
    vertical = get(attributes, :vertical, default_orientation === :vertical)
    perpendicular_length = if vertical && default_orientation === :horizontal
        (; height = _ADDON_COLORBAR_DOCK_LENGTH)
    elseif !vertical && default_orientation === :vertical
        (; width = _ADDON_COLORBAR_DOCK_LENGTH)
    else
        (;)
    end
    options = merge((; vertical), perpendicular_length, attributes)
    grid = GridLayout()
    grid.default_rowgap = Fixed(10)
    grid.default_colgap = Fixed(8)
    slot[] = grid
    target === nothing && _addon_activate_dock_tracks!(body, position)
    compact_side_dock = !vertical && default_orientation === :vertical
    colorbars = map(enumerate(scales)) do (index, scale)
        if compact_side_dock
            Label(
                grid[index, 1],
                scale.label;
                halign = :right,
                valign = :center,
                fontsize = 14
            )
            compact_options = merge(
                options,
                (;
                    label = "",
                    labelvisible = false
                )
            )
            _addon_colorbar!(grid[index, 2], scale; attributes = compact_options)
        else
            colorbar_position = vertical ? grid[1, index] : grid[index, 1]
            _addon_colorbar!(colorbar_position, scale; attributes = options)
        end
    end
    if compact_side_dock
        colsize!(grid, 1, Auto(true))
        colsize!(grid, 2, Auto(true))
    end
    return (; colorbars = Tuple(colorbars), layout = grid)
end

function _addon_bind_visibility!(figure, axes, resets, groups, status)
    for plots in values(groups), plot in plots

        on(figure.scene, plot.visible) do _
            foreach(callback -> callback(), resets)
            status[] = "Axis limits fitted to visible series"
            return nothing
        end
    end
    return groups
end

function _addon_controls!(
        shell,
        axes,
        resets,
        xsetters,
        ysetters,
        legend,
        plot_reference;
        controls::Bool
)
    widgets = Dict{Symbol, Any}()
    controls || return widgets
    column = 1
    if !isempty(axes)
        reset = _addon_button!(shell.toolbar, column, _ADDON_REFRESH_ICON)
        column += 1
        widgets[:reset] = reset
        on(shell.figure.scene, reset.clicks) do _
            foreach(callback -> callback(), resets)
            shell.status[] = "Axis limits reset"
            return nothing
        end
    end
    save = _addon_button!(shell.toolbar, column, _ADDON_SAVE_ICON)
    column += 1
    widgets[:export_svg] = save
    on(shell.figure.scene, save.clicks) do _
        plot_reference[] === nothing || LineCableModels.export_svg(plot_reference[])
        return nothing
    end
    if !isempty(xsetters)
        active = all(axis -> axis.xscale[] === Makie.log10, axes)
        toggle = Toggle(shell.toolbar[1, column]; active)
        column += 1
        Label(shell.toolbar[1, column], "log x")
        column += 1
        widgets[:xlog] = toggle
        on(shell.figure.scene, toggle.active) do enabled
            foreach(setter -> setter(enabled ? :log10 : :linear), xsetters)
            foreach(callback -> callback(), resets)
            shell.status[] = enabled ? "x-axis scale set to log" :
                             "x-axis scale set to linear"
            return nothing
        end
    end
    if !isempty(ysetters)
        active = all(axis -> axis.yscale[] === Makie.log10, axes)
        toggle = Toggle(shell.toolbar[1, column]; active)
        column += 1
        Label(shell.toolbar[1, column], "log y")
        widgets[:ylog] = toggle
        on(shell.figure.scene, toggle.active) do enabled
            foreach(setter -> setter(enabled ? :log10 : :linear), ysetters)
            foreach(callback -> callback(), resets)
            shell.status[] = enabled ? "y-axis scale set to log" :
                             "y-axis scale set to linear"
            return nothing
        end
    end
    legend === nothing || (widgets[:legend] = legend)
    return widgets
end

function _addon_figure_title!(shell, title, attributes)
    title === nothing && return nothing
    attributes isa NamedTuple || throw(ArgumentError(
        "title_attributes must be a NamedTuple",
    ))
    options = merge(
        (; tellwidth = false, halign = :center, font = :bold, fontsize = 18),
        attributes
    )
    return Label(shell.root[0, 1], String(title); options...)
end

function _addon_finish!(
        shell,
        axes,
        resets,
        xsetters,
        ysetters,
        groups,
        order,
        group_labels;
        title,
        figure_title = nothing,
        title_attributes = (;),
        legend_position,
        legend_attributes,
        legend_overflow = :ellipsis,
        legend_title = nothing,
        legend_anchor = :rt,
        panels = (),
        panel_legends = (),
        panel_group_labels = nothing,
        panel_legend_titles = nothing,
        color_scales = (),
        colorbar_position = nothing,
        colorbar_attributes = (;),
        colorbar_target = nothing,
        colorbar_target_orientation = nothing,
        controls,
        display_plot,
        export_name,
        export_theme,
        open_export
)
    title_block = _addon_figure_title!(shell, figure_title, title_attributes)
    inside_bbox = _addon_axes_viewport(
        axes,
        shell.canvas.layoutobservables.computedbbox
    )
    legend_target = nothing
    dock_targets = Dict{Any, Any}()
    resolved_colorbar_target = colorbar_target
    resolved_colorbar_orientation = colorbar_target_orientation
    shared_orientation = nothing
    if legend_position !== nothing && legend_position !== :inside &&
       legend_position == colorbar_position &&
       !isempty(color_scales) && colorbar_target === nothing
        slot, shared_orientation = _addon_legend_slot(shell.body, legend_position)
        dock = shared_orientation === :vertical ?
               GridLayout(2, 1; height = Relative(1)) :
               GridLayout(1, 2; width = Relative(1))
        slot[] = dock
        _addon_activate_dock_tracks!(shell.body, legend_position)
        if shared_orientation === :vertical
            legend_target = dock[1, 1]
            resolved_colorbar_target = dock[2, 1]
            rowsize!(dock, 1, Auto(false, 1))
            # A length that is not intrinsic must share the available space,
            # not collapse to zero when the legend occupies the same dock.
            rowsize!(dock, 2, Auto(true, 1))
        else
            legend_target = dock[1, 1]
            resolved_colorbar_target = dock[1, 2]
            colsize!(dock, 1, Auto(false, 1))
            colsize!(dock, 2, Auto(true, 1))
        end
        resolved_colorbar_orientation = shared_orientation
        dock_targets[legend_position] = (;
            target = legend_target,
            orientation = shared_orientation
        )
    end
    legend = _addon_legend!(
        shell.figure,
        shell.body,
        groups,
        order,
        group_labels;
        position = legend_position,
        attributes = legend_attributes,
        overflow = legend_overflow,
        title = legend_title,
        anchor = legend_anchor,
        inside_bbox,
        target = legend_target,
        target_orientation = shared_orientation
    )
    panel_data = isempty(panels) ? Dict{Tuple{Int, Int}, Any}() :
                 _addon_panel_legend_data(
        panels,
        axes,
        groups,
        order,
        group_labels;
        panel_labels = panel_group_labels,
        panel_titles = panel_legend_titles
    )
    panel_legend_result = _addon_panel_legends!(
        shell.figure,
        panel_data,
        panel_legends
    )
    colorbar_result = _addon_colorbars!(
        shell.body,
        color_scales;
        position = colorbar_position,
        attributes = colorbar_attributes,
        target = resolved_colorbar_target,
        target_orientation = resolved_colorbar_orientation
    )
    _addon_bind_visibility!(shell.figure, axes, resets, groups, shell.status)
    reference = Ref{Any}(nothing)
    widgets = _addon_controls!(
        shell,
        axes,
        resets,
        xsetters,
        ysetters,
        legend,
        reference;
        controls
    )
    built = LineCableModels.UIPlot(
        shell.figure,
        Tuple(axes);
        title = title_block,
        controls = widgets,
        legend,
        panel_legends = panel_legend_result.legends,
        colorbars = colorbar_result.colorbars,
        addon_state = (;
            shell,
            groups,
            order,
            labels = group_labels,
            title = Ref{Any}(legend_title),
            panel_data,
            figure_legend_position = Ref{Any}(legend_position),
            panel_legend_positions = panel_legend_result.positions,
            legend_position,
            legend_attributes,
            legend_overflow,
            legend_anchor,
            inside_bbox,
            legend_target,
            shared_orientation,
            dock_targets,
            colorbar_position,
            colorbar_target,
            colorbar_layout = colorbar_result.layout
        ),
        export_name,
        export_theme,
        open_export
    )
    reference[] = built
    display_plot && _addon_display!(shell.figure, title)
    return built
end

function _addon_figure_legend_target!(data, position)
    position in (nothing, :inside) &&
        return (; target = nothing, orientation = nothing)
    if haskey(data.dock_targets, position)
        return data.dock_targets[position]
    end
    shares_colorbar = position == data.colorbar_position &&
                      data.colorbar_target === nothing &&
                      data.colorbar_layout !== nothing
    shares_colorbar || return (; target = nothing, orientation = nothing)

    slot, orientation = _addon_legend_slot(data.shell.body, position)
    dock = orientation === :vertical ?
           GridLayout(2, 1; height = Relative(1)) :
           GridLayout(1, 2; width = Relative(1))
    slot[] = dock
    _addon_activate_dock_tracks!(data.shell.body, position)
    if orientation === :vertical
        legend_target = dock[1, 1]
        colorbar_target = dock[2, 1]
        rowsize!(dock, 1, Auto(false, 1))
        rowsize!(dock, 2, Auto(true, 1))
    else
        legend_target = dock[1, 1]
        colorbar_target = dock[1, 2]
        colsize!(dock, 1, Auto(false, 1))
        colsize!(dock, 2, Auto(true, 1))
    end
    colorbar_target[] = data.colorbar_layout
    target = (; target = legend_target, orientation)
    data.dock_targets[position] = target
    return target
end

function LineCableModels.figurelegend!(
        plot::LineCableModels.UIPlot;
        position = :right,
        overflow::Symbol = :ellipsis,
        title = missing,
        anchor = :rt,
        legend_labels = nothing,
        kwargs...
)
    data = plot.addon_state
    data === nothing && throw(ArgumentError(
        "this plot does not retain controlled legend groups",
    ))
    previous_position = data.figure_legend_position[]
    _addon_remove_legend!(plot.legend)
    if previous_position !== nothing && previous_position !== :inside &&
       previous_position != position &&
       previous_position != data.colorbar_position
        _addon_deactivate_dock_tracks!(data.shell.body, previous_position)
    end
    _addon_relabel_legend!(data.labels, data.groups, data.order, legend_labels)
    ismissing(title) || (data.title[] = title)
    resolved_title = ismissing(title) ? data.title[] : title
    dock = _addon_figure_legend_target!(data, position)
    legend = _addon_legend!(
        data.shell.figure,
        data.shell.body,
        data.groups,
        data.order,
        data.labels;
        position,
        attributes = (; kwargs...),
        overflow,
        title = resolved_title,
        anchor,
        inside_bbox = data.inside_bbox,
        target = dock.target,
        target_orientation = dock.orientation
    )
    plot.legend = legend
    data.figure_legend_position[] = legend === nothing ? nothing : position
    if plot.controls isa AbstractDict
        legend === nothing ? delete!(plot.controls, :legend) :
        (plot.controls[:legend] = legend)
    end
    return legend
end

function LineCableModels.panellegend!(
        plot::LineCableModels.UIPlot,
        logical_position::Tuple{Int, Int};
        position = :right,
        overflow::Symbol = :ellipsis,
        title = missing,
        anchor = :rt,
        legend_labels = nothing,
        kwargs...
)
    data = plot.addon_state
    data === nothing && throw(ArgumentError(
        "this plot does not retain controlled legend groups",
    ))
    haskey(data.panel_data, logical_position) || throw(BoundsError(
        collect(keys(data.panel_data)), logical_position
    ))
    panel = data.panel_data[logical_position]
    previous_position = get(data.panel_legend_positions, logical_position, nothing)
    existing = get(plot.panel_legends, logical_position, nothing)
    _addon_remove_legend!(existing)
    if previous_position !== nothing && previous_position !== :inside &&
       previous_position != position
        _addon_deactivate_dock_tracks!(panel.panel.layout, previous_position)
    end
    _addon_relabel_legend!(panel.labels, panel.groups, panel.order, legend_labels)
    ismissing(title) || (panel.title[] = title)
    resolved_title = ismissing(title) ? panel.title[] : title
    legend = _addon_legend!(
        data.shell.figure,
        panel.panel.layout,
        panel.groups,
        panel.order,
        panel.labels;
        position,
        attributes = (; kwargs...),
        overflow,
        title = resolved_title,
        anchor,
        inside_bbox = panel.axis.scene.viewport
    )
    if legend === nothing
        delete!(plot.panel_legends, logical_position)
        delete!(data.panel_legend_positions, logical_position)
    else
        plot.panel_legends[logical_position] = legend
        data.panel_legend_positions[logical_position] = position
    end
    return legend
end

function LineCableModels.figuretitle!(
        plot::LineCableModels.UIPlot,
        title;
        kwargs...
)
    data = plot.addon_state
    data === nothing && throw(ArgumentError(
        "this plot does not retain an addon figure shell",
    ))
    plot.title === nothing || delete!(plot.title)
    plot.title = title === nothing ? nothing :
                 _addon_figure_title!(data.shell, title, (; kwargs...))
    return plot.title
end

function LineCableModels.paneltitle!(
        plot::LineCableModels.UIPlot,
        logical_position::Tuple{Int, Int},
        title
)
    data = plot.addon_state
    data === nothing && throw(ArgumentError(
        "this plot does not retain logical plot panels",
    ))
    haskey(data.panel_data, logical_position) || throw(BoundsError(
        collect(keys(data.panel_data)), logical_position
    ))
    axis = data.panel_data[logical_position].axis
    axis.title[] = title === nothing ? "" : String(title)
    return axis
end
