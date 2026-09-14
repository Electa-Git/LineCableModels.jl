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
            minimum_size,
            aspect_size=Tuple(Int.(viewport.widths))
        )
        Base.display(screen, figure)
    else
        display(figure)
    end
    return figure
end

function _addon_landscape_size(size)
    size isa Tuple{Int, Int} && all(>(0), size) || throw(ArgumentError(
        "fig_size must be a tuple of two positive integers",
    ))
    width, height = size
    return (max(width, cld(4height, 3)), height)
end

function _addon_shell(; size, controls::Bool, axis::NamedTuple=(;), figure::NamedTuple=(;), kwargs...)
    axis_keys = Makie.attribute_names(Axis)
    axis_attributes = merge((; (key=>value for (key,value) in kwargs if key in axis_keys)...), axis)
    series_attributes = (; (key=>value for (key,value) in kwargs if key ∉ axis_keys)...)
    size = _addon_landscape_size(size)
    figure = Figure(; merge((; size, figure_padding=(12,12,12,12)), figure)...)
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
    status = Observable("Ready")
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
    return (; figure, root, body, canvas, toolbar, status, axis_attributes, series_attributes)
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

# Full and residual matrix cells share a decoration budget. Their outer grids
# absorb protrusions, while blank positions remain ordinary layout tracks.
# Native layout updates are guarded; exports temporarily retain this budget.
function _addon_equal_matrix_cells!(pages)
    axes = [axis for page in pages for axis in page.axes]
    isempty(axes) && return nothing
    budget = zeros(Float32, 13) # axis, figure docks, panel docks, figure title
    updating = Ref(false)
    suspended = Ref(false)
    sides = (:left, :right, :bottom, :top)
    extent(legend, position, side) = legend === nothing || position != side ? 0f0 :
        something(legend.layoutobservables.autosize[][side in (:left,:right) ? 1 : 2], 0f0)
    function reserve!(grid, sizes)
        for (side, size) in zip(sides, sizes)
            size > 0 || continue
            _addon_activate_dock_tracks!(grid, side)
            if side in (:left,:right)
                colsize!(grid, side === :left ? 1 : 3, Fixed(size))
            else
                rowsize!(grid, side === :top ? 1 : 3, Fixed(size))
            end
        end
    end
    function fit!()
        (updating[] || suspended[]) && return nothing
        desired = [maximum(getfield(axis.layoutobservables.protrusions[], side)
            for axis in axes) for side in sides]
        append!(desired, [maximum(extent(page.legend,
            page.addon_state.figure_legend_position[], side) for page in pages) for side in sides])
        append!(desired, [maximum((extent(legend,
            page.addon_state.panel_legend_positions[key], side)
            for page in pages for (key,legend) in page.panel_legends); init=0f0) for side in sides])
        push!(desired, maximum(page.title === nothing ? 0f0 :
            something(page.title.layoutobservables.autosize[][2], 0f0) for page in pages))
        all(desired .<= budget) && return nothing
        updating[] = true
        try
            budget .= max.(budget, ceil.(desired))
            P = Makie.GridLayoutBase.Protrusion
            alignment = Mixed(left=P(budget[1]), right=P(budget[2]),
                bottom=P(budget[3]), top=P(budget[4]))
            for axis in axes
                axis.alignmode[] = alignment
            end
            for page in pages
                shell = page.addon_state.shell
                reserve!(shell.body, budget[5:8])
                for cell in page.addon_state.matrix_block.cells
                    reserve!(cell.layout, budget[9:12])
                end
                if budget[13] > 0
                    page.title === nothing && (shell.root[0,1] = GridLayout())
                    rowsize!(shell.root, 0, Fixed(budget[13]))
                end
            end
        finally
            updating[] = false
        end
        return nothing
    end
    for page in pages, axis in page.axes
        on(page.figure.scene, axis.layoutobservables.protrusions) do _
            fit!()
        end
    end
    for page in pages
        for block in (page.title, page.legend, values(page.panel_legends)...)
            block === nothing && continue
            on(page.figure.scene, block.layoutobservables.autosize) do _
                fit!()
            end
        end
        page.addon_state = merge(page.addon_state,
            (matrix_layout=(; suspended, refit=fit!),))
    end
    fit!()
    return nothing
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

function _addon_refit_matrix_block!(plot, block)
    data = plot.addon_state
    haskey(data, :matrix_layout) || return nothing
    fit! = data.matrix_layout.refit
    if block !== nothing
        on(block.blockscene, block.layoutobservables.autosize) do _
            fit!()
        end
    end
    fit!()
    return nothing
end

function _addon_scale(symbol::Symbol)
    symbol === :linear && return Makie.identity
    symbol === :log10 && return Makie.log10
    # Same signed-log scale, with no cancellation in its linear neighbourhood.
    symbol === :pseudolog10 && return Makie.ReversibleScale(
        x -> sign(x) * log1p(abs(x)) / log(10),
        x -> sign(x) * expm1(abs(x) * log(10));
        limits=(0.0f0, 3.0f0), name=:pseudolog10)
    throw(ArgumentError("unsupported axis scale :$symbol"))
end
_addon_scale(scale) = scale

function _addon_scientific_exponent(values)
    magnitudes = Float64[]
    for value in values
        numeric = LineCableModels.nominal(value)
        numeric isa Real || continue
        converted = Float64(numeric)
        isfinite(converted) && !iszero(converted) && push!(magnitudes, abs(converted))
    end
    isempty(magnitudes) && return nothing
    return 3fld(floor(Int, log10(maximum(magnitudes))), 3)
end

function _addon_linear_tickformat(exponent::Int)
    # Normalize subnormal magnitudes without underflowing the power of ten.
    shift = exponent < -307 ? 308 : 0
    scale = 10.0^(exponent + shift)
    return function (values)
        mantissas = [(Float64(value) * 10.0^shift) / scale for value in values]
        # Signed zeros are the same tick. A range transition can temporarily
        # underflow values through the preceding exponent before it is replaced.
        finite = sort!(unique(iszero(value) ? 0.0 : value
            for value in mantissas if isfinite(value)))
        isempty(finite) && return string.(mantissas)
        magnitude = maximum(abs, finite)
        digits = iszero(magnitude) ? 0 : max(0, 3 - floor(Int, log10(magnitude)))
        if length(finite) > 1
            # Retain distinct ticks when zooming into a narrow, offset interval.
            spacing = minimum(diff(finite))
            digits = max(digits, 1 - floor(Int, log10(spacing)))
        end
        digits = clamp(digits, 0, 17)
        return map(mantissas) do value
            label = @sprintf("%.*f", digits, value)
            digits > 0 && (label = rstrip(rstrip(label, '0'), '.'))
            label == "-0" ? "0" : label
        end
    end
end

function _addon_decade_ticks(vmin, vmax, count::Int)
    isfinite(vmin) && isfinite(vmax) && 0 < vmin <= vmax || return Float64[]
    span = log10(vmax)-log10(vmin)
    if span < 2
        isapprox(vmin,vmax;rtol=sqrt(eps(Float64)),atol=0) && return unique([vmin,vmax])
        ticks = Float64[]
        for exponent in floor(Int,log10(vmin)):floor(Int,log10(vmax))
            lower,upper = max(vmin,10.0^exponent),min(vmax,10.0^(exponent+1))
            lower < upper || continue
            budget = max(2,ceil(Int,count*(log10(upper)-log10(lower))/span)+1)
            # Locator arithmetic is local to one decade, including extreme SI
            # magnitudes. Only the returned positions use the published units.
            factor = 10.0^clamp(exponent,-307,307)
            values = Makie.get_tickvalues(Makie.LinearTicks(budget),lower/factor,upper/factor)
            append!(ticks,filter(x -> isfinite(x) && lower <= x <= upper,values.*factor))
            vmin <= 10.0^exponent <= vmax && push!(ticks,10.0^exponent)
        end
        isempty(ticks) && append!(ticks,(vmin,vmax))
        return sort!(unique!(ticks))
    end
    first_exponent = ceil(Int, log10(vmin))
    last_exponent = floor(Int, log10(vmax))
    step = max(1, cld(last_exponent - first_exponent, max(1, count - 1)))
    return 10.0 .^ (first_exponent:step:last_exponent)
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

function _addon_set_axis!(entries::AbstractVector, dim::Symbol, scale=nothing)
    dim in (:x, :y) || throw(ArgumentError("axis dimension must be :x or :y"))
    index = dim === :x ? 1 : 2
    # Resolve and validate the complete page before any native observable changes.
    # These are native axis bindings, not another interpretation of result data.
    targets = map(entries) do entry
        target = _addon_scale(scale === nothing ? entry.scale : scale)
        context = "axis :$dim ($(repr(entry.axis.title[])))"
        requested = entry.axis.limits[]
        requested = length(requested) == 4 ? (requested[1:2], requested[3:4]) : requested
        bounds = requested[index] === nothing ? () : requested[index]
        all(value -> value === nothing || isfinite(value), bounds) ||
            throw(DomainError(bounds, "$context requires finite explicit limits"))
        values = _addon_visible_values(entry.axis, dim)
        if target === Makie.log10 && entry.signed && !isempty(values) && !all(>(0),values)
            target = _addon_scale(:pseudolog10)
        end
        if target === Makie.log10
            all(>(0), values) && all(value -> value === nothing || value > 0, bounds) ||
                throw(DomainError(bounds, "logarithmic $context requires positive visible data, uncertainty bounds and explicit limits"))
        end
        if !isempty(values)
            lower, upper = extrema(values)
            if isapprox(lower, upper; rtol=sqrt(eps(Float64)), atol=0)
                lower, upper = _addon_constant_limits(values, values, target === Makie.log10)
            end
            explicit = isempty(bounds) ? (nothing, nothing) : bounds
            lower, upper = something(explicit[1], lower), something(explicit[2], upper)
            transformed = (target(lower), target(upper))
            all(isfinite, transformed) && transformed[1] < transformed[2] ||
                throw(DomainError((lower, upper), "$context requires distinct finite transformed limits"))
        end
        target
    end
    for (entry, target) in zip(entries, targets)
        axis = entry.axis
        previous = axis.targetlimits[]
        getproperty(axis, Symbol(dim, :scale))[] = target
        entry.reset(; xauto=dim === :x, yauto=dim === :y)
        # Makie refits both dimensions on transform changes. Restore the unrelated
        # view without turning its interactive limits into a caller request.
        other = 3 - index
        current = axis.targetlimits[]
        origin, widths = collect(current.origin), collect(current.widths)
        origin[other], widths[other] = previous.origin[other], previous.widths[other]
        axis.targetlimits[] = Makie.Rect2d(origin..., widths...)
    end
    return entries
end

function _addon_numeric_values(values)
    # Undefined observations remain missing in publications. Makie's numeric
    # line boundary uses NaN gaps, including an entirely undefined phase trace.
    nominal_values = map(value -> ismissing(value) ? NaN : LineCableModels.nominal(value), values)
    errors = LineCableModels.uncertainty.(values)
    return nominal_values, any(error -> !iszero(error), errors) ? errors : nothing
end

function _addon_line!(axis, xdata, ydata; dependent_plots, label, color = nothing, visible = true)
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
    append!(dependent_plots, (plot => first(plots) for plot in Iterators.drop(plots, 1)))
    return plots
end

function _addon_visible_values(series, dim::Symbol; include_uncertainty::Bool = false)
    values = Float64[]
    for item in series
        first(item.plots).visible[] || continue
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

function _addon_visible_values(axis::Axis, dim::Symbol)
    index = dim === :x ? 1 : 2
    bounds = Makie.data_limits(axis.scene, plot ->
        !to_value(get(plot, :visible, true)) ||
        !to_value(get(plot, Symbol(dim, :autolimits), true)) ||
        to_value(get(plot, :space, :data)) !== :data)
    lower, upper = bounds.origin[index], bounds.origin[index] + bounds.widths[index]
    return isfinite(lower) && isfinite(upper) ? [lower, upper] : Float64[]
end

function LineCableModels.plotwindow(
        callback::F;
        title::AbstractString,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        series_attributes = nothing,
        size::Tuple{Int, Int} = (800, 400),
        layout = nothing,
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        export_name::AbstractString = title,
        kwargs...
) where {F}
    _addon_activate_backend(backend)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls, kwargs...)
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
        # Only caller-supplied overrides apply to caller-constructed axes. Keep
        # native defaults/conversions intact; scale changes use common preflight.
        for axis in axes, (key,value) in pairs(shell.axis_attributes)
            key in (:xscale,:yscale) || setproperty!(axis,key,value)
        end
        requested_scales = !any(key -> key in (:xscale,:yscale), keys(shell.axis_attributes)) ? nothing :
            [(x=get(shell.axis_attributes,:xscale,axis.xscale[]),
                y=get(shell.axis_attributes,:yscale,axis.yscale[])) for axis in axes]
        resets = Function[_addon_reset!(axis) for axis in axes]
        native = series_attributes === nothing && isempty(shell.series_attributes) ? Any[] :
            Any[handle for axis in axes for handle in axis.scene.plots]
        order = [Symbol("series_$index") for index in eachindex(native)]
        groups = Dict(group => Any[handle] for (group, handle) in zip(order, native))
        _addon_finish!(
            shell,
            axes,
            resets,
            groups,
            order,
            Dict(group => string(group) for group in order);
            requested_scales,
            series_attributes,
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
        series_attributes = nothing,
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
        ylabel = nothing,
        xscale = :linear,
        yscale = :linear,
        kwargs...
) where {F}
    _addon_activate_backend(backend)
    resolved_panel_titles = _addon_panel_titles(panel_titles, 1)
    panel_title = resolved_panel_titles === nothing ? title :
                  only(resolved_panel_titles)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size = fig_size, controls, kwargs...)
        panel = _addon_panel!(shell, (1, 1))
        axis,scales = _addon_axis!(
            panel.content,
            xobservation,
            yobservation;
            title = panel_title,
            xscale,
            yscale,
            xlabel,
            ylabel,
            native_attributes=shell.axis_attributes
        )
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, String}()
        series = NamedTuple[]
        draw(axis, groups, order, labels, series)
        _addon_relabel_legend!(labels, groups, order, legend_labels)
        reset! = _addon_reset!(axis, series)
        _addon_finish!(
            shell,
            Any[axis],
            Function[reset!],
            groups,
            order,
            labels;
            requested_scales=(scales,),
            series_attributes,
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

function _addon_constant_limits(values, interval_values, logarithmic::Bool)
    if logarithmic
        all(>(0), interval_values) || throw(DomainError(
            interval_values,
            "logarithmic axes require strictly positive data"
        ))
        center = sum(log, extrema(values))/2
        halfspan = max(log(1.05),2maximum(value -> abs(log(value)-center),interval_values))
        lower,upper = exp(center-halfspan),exp(center+halfspan)
        isfinite(upper) && 0 < lower < upper && log10(lower) < log10(upper) ||
            throw(DomainError((lower,upper),"automatic logarithmic bounds are not representable"))
        return lower,upper
    end
    all(iszero, interval_values) && return (-1.0, 1.0)
    lower, upper = extrema(values)
    center = lower / 2 + upper / 2
    base_halfspan = max(0.05abs(center), eps(center))
    interval_halfspan = maximum(abs(value - center) for value in interval_values)
    halfspan = max(base_halfspan, 2interval_halfspan)
    return center - halfspan, center + halfspan
end

function _addon_axis_format!(axis)
    for (index, dim) in enumerate((:x, :y))
        scale = getproperty(axis, Symbol(dim, :scale))
        ticks = getproperty(axis, Symbol(dim, :ticks))
        tickformat = getproperty(axis, Symbol(dim, :tickformat))
        label = getproperty(axis, Symbol(dim, :label))
        labelsize = getproperty(axis, Symbol(dim, :ticklabelsize))
        labelfont = getproperty(axis, Symbol(dim, :ticklabelfont))
        rotation = getproperty(axis, Symbol(dim, :ticklabelrotation))
        conversion = getproperty(axis, Symbol(:dim, index, :_conversion))
        lineaxis = getproperty(axis, Symbol(dim, :axis))
        raw_label = Ref{Any}(label[])
        rendered_label = Ref{Any}(label[])
        installed_ticks = Ref{Any}(Makie.automatic)
        installed_format = Ref{Any}(Makie.automatic)
        installed_mode = Ref{Any}(nothing)
        updating = Ref(false)
        # Reuse one native text measurement per dimension, outside the data
        # scene. Font changes and long, precise labels affect available density.
        probe = text!(axis.blockscene, 0, 0; text="", markerspace=:data,
            visible=false, inspectable=false)
        function update!()
            updating[] && return nothing
            updating[] = true
            try
                limits = axis.finallimits[]
                current_scale, current_ticks = scale[], ticks[]
                current_format, current_label = tickformat[], label[]
                # A transform notification precedes Makie's own limit reset.
                # Do not send its old, possibly negative linear view to a log locator.
                current_scale === Makie.log10 && limits.origin[index] <= 0 && return nothing
                label_changed = current_label !== rendered_label[]
                label_changed && (raw_label[] = current_label)
                owned_ticks = current_ticks === installed_ticks[] || current_ticks === Makie.automatic
                owned_format = current_format === installed_format[] || current_format === Makie.automatic
                numeric = conversion[] === nothing
                lower, upper = limits.origin[index], limits.origin[index] + limits.widths[index]
                decades = current_scale === Makie.log10 && 0 < lower < upper &&
                    log10(upper) - log10(lower) >= 2
                exponent = something(_addon_scientific_exponent((lower, upper)), 0)
                signed_linear = current_scale === _addon_scale(:pseudolog10) &&
                    max(abs(lower), abs(upper)) < 1
                mode = if numeric && (current_scale === Makie.identity || signed_linear) &&
                        (owned_ticks || current_ticks isa AbstractVector{<:Real})
                    (:linear, exponent)
                elseif numeric && current_scale === Makie.log10 &&
                        (owned_ticks || current_ticks isa AbstractVector{<:Real})
                    owned_ticks && decades ? (:log10, 0) : (:linear, exponent)
                else
                    nothing
                end
                if owned_format
                    if mode != installed_mode[] || current_format === Makie.automatic
                        installed_format[] = if mode === nothing
                            Makie.automatic
                        elseif first(mode) === :linear
                            _addon_linear_tickformat(exponent)
                        else
                            # Native LineAxis may still hold its previous limits
                            # during a scale notification. Its temporary zero tick
                            # must not throw before the positive locator replaces it.
                            values -> [value <= 0 ? string(value) : Makie.rich("10", Makie.superscript(
                                replace(string(round(Int, log10(value))), "-" => "−");
                                offset=Makie.Vec2f(0.1, 0.0))) for value in values]
                        end
                        installed_mode[] = mode
                        tickformat[] === installed_format[] || (tickformat[] = installed_format[])
                    end
                end
                formatted = owned_format && mode !== nothing && first(mode) === :linear ?
                    _addon_axis_label(raw_label[], exponent, :linear) : raw_label[]
                if label_changed || !isequal(formatted, current_label)
                    rendered_label[] = formatted
                    label[] = formatted
                end
                if owned_ticks
                    pixels = axis.scene.viewport[].widths[index]
                    spacing = (index == 1 && current_scale === Makie.identity ? 5.5 : 3.0) * labelsize[]
                    count = clamp(floor(Int, pixels / max(1, spacing)), 3, 10)
                    probe.font[] = labelfont[]
                    probe.fontsize[] = labelsize[]
                    probe.rotation[] = rotation[]
                    # Native locator/formatter notification is synchronous. Fit
                    # its rendered strings, reducing only automatic tick density.
                    for _ in 1:9
                        selected = if !numeric
                            Makie.automatic
                        elseif current_scale === Makie.identity
                            # A native fit can briefly expose equal or adjacent
                            # Float64 endpoints before automatic padding runs.
                            # Preserve explicit narrow zooms too, without asking
                            # LinearTicks to subdivide an unrepresentable interval.
                            let count=count
                                (lo, hi) -> isapprox(lo, hi; rtol=sqrt(eps(Float64)), atol=0) ?
                                    unique([lo, hi]) : Makie.get_tickvalues(Makie.LinearTicks(count), lo, hi)
                            end
                        elseif current_scale === Makie.log10
                            _addon_decade_ticks(lower, upper, count)
                        elseif current_scale === _addon_scale(:pseudolog10)
                            # Native PseudologTicks dispatches on Makie's scale
                            # instance. Reuse its placement, not its cancelling
                            # transform, and pass numeric positions to this axis.
                            if isapprox(lower, upper; rtol=sqrt(eps(Float64)), atol=0)
                                unique([lower, upper])
                            else
                                locator = signed_linear ? Makie.LinearTicks(count) : Makie.PseudologTicks(count)
                                first(Makie.get_ticks(locator, Makie.pseudolog10,
                                    Makie.automatic, lower, upper))
                            end
                        else
                            Makie.automatic
                        end
                        if !isequal(selected, ticks[])
                            installed_ticks[] = selected
                            ticks[] = selected
                        end
                        selected === Makie.automatic && break
                        if current_scale === Makie.log10 && !decades
                            # Physical-value log ticks are NOT equally spaced on
                            # screen. Fit adjacent rendered labels, not count/width.
                            positions = Float64[]
                            extents = Float64[]
                            for (value,text) in zip(selected,lineaxis.ticklabels[])
                                value > 0 || continue
                                probe.text[] = text
                                push!(positions,pixels*(log10(value)-log10(lower))/(log10(upper)-log10(lower)))
                                push!(extents,Makie.boundingbox(probe,:data).widths[index])
                            end
                            if length(positions) == length(selected)
                                retained = Int[]
                                for i in eachindex(positions)
                                    if isempty(retained) || positions[i]-positions[last(retained)] >=
                                            (extents[i]+extents[last(retained)])/2+labelsize[]/2
                                        push!(retained,i)
                                    end
                                end
                                fitted_ticks = selected[retained]
                                if !isequal(ticks[],fitted_ticks)
                                    installed_ticks[] = fitted_ticks
                                    ticks[] = fitted_ticks
                                end
                            end
                            break
                        end
                        extent = 0.0
                        for text in lineaxis.ticklabels[]
                            probe.text[] = text
                            extent = max(extent, Makie.boundingbox(probe, :data).widths[index])
                        end
                        fitted = max(2, floor(Int, pixels / max(spacing, extent + labelsize[])))
                        fitted >= count && break
                        count = fitted
                    end
                end
            finally
                updating[] = false
            end
            return nothing
        end
        # Run before native tick conversion: a newly assigned labelled tuple
        # cannot be consumed with the previously installed numeric formatter.
        onany((_...) -> update!(), axis.scene, ticks, tickformat; priority=1)
        # Range/transform updates must instead follow native LineAxis propagation;
        # changing its formatter while it still has the old limits is unsafe.
        onany((_...) -> update!(), axis.scene, axis.finallimits, scale, label,
            axis.scene.viewport, labelsize, labelfont, rotation, conversion;
            priority=-3, update=true)
    end
    return axis
end

# Bind the native limit lifecycle once and return this axis's reset action.
function _addon_reset!(axis, series=())
    # Own numeric ticks before the first data fit, including its synchronous
    # native callbacks. Every recipe and caller-owned plotwindow uses this bind.
    _addon_axis_format!(axis)
    fitting = Ref(false)
    corrections = Any[nothing, nothing]
    function reset!(; xauto::Bool=true, yauto::Bool=true)
        fitting[] && return axis
        fitting[] = true
        try
            reset_limits!(axis; xauto, yauto)
            requested = axis.limits[]
            requested = length(requested) == 4 ? (requested[1:2], requested[3:4]) : requested
            for (index, dim) in enumerate((:x, :y))
                (index == 1 ? xauto : yauto) || continue
                corrections[index] = nothing
                explicit = requested[index] === nothing ? (nothing, nothing) : requested[index]
                all(value -> value !== nothing, explicit) && continue
                interval_values = _addon_visible_values(axis, dim)
                isempty(interval_values) && continue
                values = isempty(series) ? interval_values : _addon_visible_values(series, dim)
                isempty(values) && continue
                isapprox(extrema(values)...; rtol=sqrt(eps(Float64)), atol=0) || continue
                if !isempty(series)
                    expected = _addon_visible_values(series, dim; include_uncertainty=true)
                    # Native extra plots or independently hidden error bars own
                    # their actual extents, not the original observation array.
                    if !all(isapprox.(interval_values, collect(extrema(expected))))
                        isapprox(extrema(interval_values)...; rtol=sqrt(eps(Float64)), atol=0) || continue
                        values = interval_values
                    end
                end
                any(plot -> haskey(plot, :model) &&
                    any(j -> plot.model[][index, j] != (index == j), 1:4),
                    axis.scene.plots) && continue
                scale = getproperty(axis, Symbol(dim, :scale))[]
                limits = _addon_constant_limits(values, interval_values, scale === Makie.log10)
                lower = something(explicit[1], limits[1])
                upper = something(explicit[2], limits[2])
                origin, widths = collect(axis.targetlimits[].origin), collect(axis.targetlimits[].widths)
                original = (origin[index], widths[index])
                origin[index], widths[index] = lower, upper - lower
                corrections[index] = (; requested=requested[index], scale, original,
                    fitted=(origin[index], widths[index]))
                axis.targetlimits[] = Makie.Rect2d(origin..., widths...)
            end
        finally
            fitting[] = false
        end
        return axis
    end
    # Makie also refits before displaying a Figure. Preserve our degenerate-data
    # padding without storing automatic limits as user requests. Only the exact
    # native auto-fit for the same scale/request is corrected; zooms are untouched.
    # This callback uses two cached bounds, not a data scan on every view update.
    on(axis.scene, axis.targetlimits; priority=1) do view
        requested = axis.limits[]
        requested = length(requested) == 4 ? (requested[1:2], requested[3:4]) : requested
        origin, widths = collect(view.origin), collect(view.widths)
        for (index, dim) in enumerate((:x, :y))
            explicit = requested[index] === nothing ? (nothing, nothing) : requested[index]
            # Native nonlinear auto-fitting leaves an all-zero series at (0,0).
            # Repair that automatic singular view before the camera's reciprocal
            # scaling runs, including during the transform-triggered first fit.
            # This uses only the view bounds, never a new data scan on zoom.
            if iszero(widths[index]) && any(isnothing, explicit) &&
                    getproperty(axis, Symbol(:dim, index, :_conversion))[] === nothing
                bounds = _addon_constant_limits((origin[index],), (origin[index],),
                    getproperty(axis, Symbol(dim, :scale))[] === Makie.log10)
                lower, upper = something(explicit[1], bounds[1]), something(explicit[2], bounds[2])
                origin[index], widths[index] = lower, upper-lower
            end
            fitting[] && continue
            correction = corrections[index]
            correction === nothing && continue
            if isequal(requested[index], correction.requested) &&
                    getproperty(axis, Symbol(dim, :scale))[] === correction.scale &&
                    (origin[index], widths[index]) == correction.original
                origin[index], widths[index] = correction.fitted
            end
        end
        fitted = Makie.Rect2d(origin..., widths...)
        fitted == view || (axis.targetlimits[] = fitted)
        return nothing
    end
    on(_ -> reset!(), axis.scene, axis.limits; priority=-3)
    reset!()
    return reset!
end

function _addon_axis!(
        position,
        xobservation,
        yobservation;
        title,
        xscale,
        yscale,
        xlabel = nothing,
        ylabel = nothing,
        attributes = (;),
        native_attributes = (;)
)
    xaxis_label = xlabel === nothing ?
                  LineCableModels.Units.label(xobservation.quantity, xobservation.unit) :
                  String(xlabel)
    yaxis_label = ylabel === nothing ?
                  LineCableModels.Units.label(yobservation.quantity, yobservation.unit) :
                  String(ylabel)
    options = merge(
        (;
            title,
            tellwidth = false,
            tellheight = false,
            xscale = _addon_scale(xscale),
            yscale = _addon_scale(yscale),
            xlabel = xaxis_label,
            ylabel = yaxis_label
        ),
        attributes, native_attributes)
    scales = (x=_addon_scale(options.xscale),y=_addon_scale(options.yscale))
    # Draw on safe axes, then validate complete native extents in the common
    # finish before applying requested transforms to any axis on the page.
    axis = Axis(position; merge(options,(xscale=identity,yscale=identity))...)
    return axis,scales
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

function _addon_legend_sources!(legend, groups, dependents)
    legend === nothing && return nothing
    owners = IdDict{Any,Any}(dependents)
    for handles in values(groups), handle in handles
        get!(owners, handle, handle)
    end
    # Preserve Makie's glyphs, but target the registered owning plots. Composite
    # glyphs may refer to derived child attributes that are not writable inputs.
    # Deduplicate across the entire entry: several glyphs still mean one action.
    for (_, entries) in legend.entrygroups[], entry in entries
        seen = Base.IdSet{Any}()
        for element in entry.elements
            # Cairo's LineSegments renderer adds joinstyle=nothing to the
            # source graph. Native legend extraction then mistakes that cache
            # for a line style on recreation. Keep the fallback in the glyph,
            # without modifying the source graph or overriding a real style.
            if element isa LineElement && to_value(element.joinstyle) === nothing
                element.attributes[:joinstyle] = legend.joinstyle
            end
            targets = Makie.get_plots(element)
            resolved = Makie.Plot[]
            for plot in targets
                owner = plot
                while owner isa Makie.Plot && !haskey(owners, owner)
                    owner = owner.parent
                end
                owner = get(owners, owner, plot)
                while haskey(owners, owner) && owners[owner] !== owner
                    owner = owners[owner]
                end
                if owner ∉ seen
                    push!(seen, owner)
                    push!(resolved, owner)
                end
            end
            empty!(targets)
            append!(targets, resolved)
        end
    end
    # Rebuild native listeners after changing targets. Also initialise their
    # shades when a hidden entry is recreated or reappears after overflow.
    on(legend.blockscene, legend.entrygroups; priority=-1) do entrygroups
        for (_, entries) in entrygroups, entry in entries
            foreach(notify, Makie.get_plot_visibilities(entry))
        end
    end
    notify(legend.entrygroups)
    return legend
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
    fit!(bounding_box[])
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

function _addon_wrap_legend_label(label, width, measure)
    lines = String[]
    for paragraph in split(label, '\n'; keepempty=true)
        current = ""
        for word in split(paragraph)
            candidate = isempty(current) ? String(word) : "$current $word"
            if !isempty(current) && measure(candidate) > width
                push!(lines, current)
                current = ""
            end
            # Very long identifiers also need a lossless fallback. No formula
            # field is elided merely because its token has no spaces.
            for character in (isempty(current) ? String(word) : " $word")
                candidate = current * character
                if !isempty(current) && measure(candidate) > width
                    push!(lines, current)
                    current = ""
                end
                current *= character
            end
        end
        push!(lines, current)
    end
    return join(lines, '\n')
end

# Native Legend already owns the row-major grid, text rendering and click
# targets. Only its bank count and (when necessary) label line breaks change.
function _addon_grid_legend!(figure, bounds, legend)
    entries = last(only(legend.entrygroups[]))
    originals = [String(entry.label[]) for entry in entries]
    fitting = Ref(false)
    previous_metrics = Ref{Any}(nothing)
    # Match native Label: glyphs and positions must share data space for the
    # bounding box to include font extents rather than just the anchor point.
    probe = text!(legend.blockscene, 0, 0; text="", markerspace=:data,
        visible=false, inspectable=false)
    function fit!()
        fitting[] && return nothing
        available = Float64(bounds[].widths[1]) - sum(legend.margin[][1:2]) -
            sum(legend.padding[][1:2]) - 4
        available > 0 || return nothing
        fitting[] = true
        try
            widths = Float64[]
            for (entry, original) in zip(entries, originals)
                probe.font[] = entry.labelfont[]
                probe.fontsize[] = entry.labelsize[]
                measured = Dict{String,Float64}()
                measure = text -> get!(measured, text) do
                    probe.text[] = text
                    Float64(Makie.boundingbox(probe, :data).widths[1])
                end
                patch = Float64(entry.patchsize[][1]) + legend.patchlabelgap[]
                wrapped = measure(original) <= available-patch ? original :
                    _addon_wrap_legend_label(original, max(1.0, available-patch), measure)
                entry.label[] == wrapped || (entry.label[] = wrapped)
                push!(widths, patch + measure(wrapped))
            end
            columns = 1
            for count in length(entries):-1:1
                required = sum(maximum(widths[column:count:end]) for column in 1:count) +
                    (count-1)*legend.colgap[]
                if required <= available
                    columns = count
                    break
                end
            end
            metrics = (Tuple(widths), Tuple(entry.label[] for entry in entries),
                Tuple((entry.labelsize[],entry.labelfont[]) for entry in entries))
            if legend.nbanks[] != columns
                legend.nbanks[] = columns
            elseif previous_metrics[] != metrics
                # Native Legend relayout is triggered by bank/layout controls,
                # not by changed label extents alone.
                notify(legend.nbanks)
            end
            previous_metrics[] = metrics
        finally
            fitting[] = false
        end
        return nothing
    end
    on(legend.blockscene, bounds) do _
        fit!()
    end
    onany((_...) -> fit!(), legend.blockscene, legend.labelsize, legend.labelfont,
        legend.padding, legend.margin, legend.colgap, legend.patchlabelgap, legend.patchsize)
    for (index, entry) in enumerate(entries)
        on(legend.blockscene, entry.label) do label
            fitting[] && return
            originals[index] = String(label)
            fit!()
        end
    end
    fit!()
    return legend
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
        dependent_plots,
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
            legend = Legend(
                figure,
                entries,
                displayed,
                title;
                options...
            )
            return _addon_legend_sources!(legend, groups, dependent_plots)
        end
        ellipsis = LineElement(color = :transparent)
        legend = Legend(
            figure,
            Any[entries..., ellipsis],
            [displayed; "(...)"],
            title;
            options...
        )
        _addon_legend_sources!(legend, groups, dependent_plots)
        return _addon_responsive_legend!(
            figure,
            # Available space, not the legend's content-dependent size.
            legend.layoutobservables.suggestedbbox,
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
    if position in (:top, :bottom) &&
            !haskey(attributes, :orientation) && !haskey(attributes, :nbanks)
        # In Makie, vertical orientation with nbanks columns fills row first.
        # The legend's dock still determines tellheight/tellwidth, not this
        # native storage orientation.
        options = merge(options, (; orientation=:vertical, halign=:center))
        legend = Legend(legend_grid[1, 1], entries, displayed, title; options...)
        _addon_legend_sources!(legend, groups, dependent_plots)
        return _addon_grid_legend!(figure, legend_grid.layoutobservables.computedbbox, legend)
    end
    if overflow === :show_all
        legend = Legend(legend_grid[1, 1], entries, displayed, title; options...)
        return _addon_legend_sources!(legend, groups, dependent_plots)
    end
    ellipsis = LineElement(color = :transparent)
    legend = Legend(
        legend_grid[1, 1],
        Any[entries..., ellipsis],
        [displayed; "(...)"],
        title;
        options...
    )
    _addon_legend_sources!(legend, groups, dependent_plots)
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

function _addon_panel_legends!(figure, panel_data, requested, dependent_plots)
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
            dependent_plots,
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
            finite_boxes = filter(
                box -> isfinite(box.origin[dimension]) &&
                       isfinite(box.widths[dimension]),
                boxes)
            before = visible ?
                     ceil(maximum(
                box -> -box.origin[dimension], finite_boxes; init = 0.0
            )) : 0.0
            after = visible ?
                    ceil(maximum(
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
        previous = Ref(plot.visible[])
        on(figure.scene, plot.visible) do visible
            # A legend rebuild can notify without changing visibility. Keep
            # the current zoom and avoid fitting identical series repeatedly.
            visible == previous[] && return nothing
            previous[] = visible
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
        plot_reference[] === nothing && return nothing
        try
            output = LineCableModels.export_svg(plot_reference[])
            shell.status[] = "Saved SVG to $output"
        catch exception
            exception isa Union{ArgumentError, SystemError, Base.IOError} || rethrow()
            shell.status[] = sprint(showerror, exception)
        end
        return nothing
    end
    for (dim, setters) in ((:x, xsetters), (:y, ysetters))
        isempty(setters) && continue
        active = all(entry -> getproperty(entry.axis, Symbol(dim, :scale))[] in
            (Makie.log10, _addon_scale(:pseudolog10)), setters)
        toggle = Toggle(shell.toolbar[1, column]; active)
        column += 1
        caption = Label(shell.toolbar[1, column], "log $dim")
        column += 1
        widgets[Symbol(dim, :log)] = toggle
        changing = Ref(false)
        for entry in setters
            on(shell.figure.scene,getproperty(entry.axis,Symbol(dim,:scale));update=true) do _
                kinds = unique(getproperty(item.axis,Symbol(dim,:scale))[] === log10 ?
                    "log" : getproperty(item.axis,Symbol(dim,:scale))[] === identity ?
                    "linear" : getproperty(item.axis,Symbol(dim,:scale))[] === _addon_scale(:pseudolog10) ?
                    "signed log" : "custom" for item in setters)
                mode = length(kinds)==1 ? only(kinds) : "mixed log"
                caption.text[] = mode == "linear" ? "log $dim" : "$mode $dim"
                if !changing[]
                    changing[] = true
                    try
                        active = all(kind -> kind in ("log","signed log"),kinds)
                        toggle.active[] == active || (toggle.active[] = active)
                    finally
                        changing[] = false
                    end
                end
            end
        end
        on(shell.figure.scene, toggle.active) do enabled
            changing[] && return nothing
            changing[] = true
            try
                _addon_set_axis!(setters, dim, enabled ? :log10 : :linear)
                shell.status[] = enabled ? "Axis scale set to $(caption.text[])" : "$dim-axis scale set to linear"
            catch exception
                if exception isa Union{ArgumentError, DomainError}
                    toggle.active[] = !enabled
                    shell.status[] = sprint(showerror, exception)
                    return nothing
                end
                rethrow()
            finally
                changing[] = false
            end
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
        groups,
        order,
        group_labels;
        scale_controls::Bool=true,
        signed_ylog::Bool=false,
        requested_scales=nothing,
        dependent_plots = Pair{Makie.Plot,Makie.Plot}[],
        title,
        figure_title = nothing,
        title_attributes = (;),
        series_attributes = nothing,
        series_defaults = nothing,
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
    isempty(axes) && !isempty(shell.axis_attributes) && throw(ArgumentError(
        "native Axis attributes require a figure containing an Axis"))
    append!(dependent_plots, _addon_series_styles!(groups, order, series_attributes;
        defaults=series_defaults, shared=shell.series_attributes))
    setters = map(enumerate((:x,:y))) do (index,dim)
        entries = [(;axis,reset,signed=index==2 && signed_ylog) for (axis,reset) in zip(axes,resets)]
        if requested_scales !== nothing
            _addon_set_axis!([merge(entry,(scale=getproperty(scales,dim),))
                for (entry,scales) in zip(entries,requested_scales)],dim)
        end
        scale_controls || return NamedTuple[]
        filter(entries) do entry
            getproperty(entry.axis,Symbol(:dim,index,:_conversion))[] === nothing &&
                getproperty(entry.axis,Symbol(dim,:scale))[] in
                    (identity,log10,_addon_scale(:pseudolog10)) &&
                !isempty(_addon_visible_values(entry.axis,dim))
        end
    end
    xsetters,ysetters = setters
    for (dependent, owner) in dependent_plots
        dependent.visible[] = owner.visible[]
        previous = Ref(owner.visible[])
        on(shell.figure.scene, owner.visible) do visible
            # Legend relayout re-emits the current state to initialise shading;
            # it is not a series action and must preserve native component edits.
            visible == previous[] && return nothing
            previous[] = visible
            dependent.visible[] = visible
        end
    end
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
        dependent_plots,
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
        panel_legends,
        dependent_plots
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
            dependent_plots,
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
        dependent_plots = data.dependent_plots,
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
    _addon_refit_matrix_block!(plot, legend)
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
        dependent_plots = data.dependent_plots,
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
    _addon_refit_matrix_block!(plot, legend)
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
    _addon_refit_matrix_block!(plot, plot.title)
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
