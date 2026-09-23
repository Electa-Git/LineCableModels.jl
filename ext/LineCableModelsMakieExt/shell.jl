const _ADDON_BUTTON_SIZE = 32
const _ADDON_BUTTON_BACKGROUND = Makie.RGBf(0.94, 0.94, 0.94)
const _ADDON_ICON_COLOR = Makie.RGBAf(0.15, 0.15, 0.15, 1.0)
const _ADDON_COLORBAR_DOCK_LENGTH = 140
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

function _addon_theme(; export_theme::Symbol = :default)
    export_theme in (:default, :publication) || throw(ArgumentError(
        "export_theme must be :default or :publication",
    ))
    return Theme(
        backgroundcolor = :grey90,
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

function _addon_shell(;
        size, controls::Bool, axis::NamedTuple = (;), figure::NamedTuple = (;),
        widgets = (), guide_gap = 8, guide_spacing = (;),
        colorbar_position = _omitted, colorbar_attributes = (;), colorbar_group_attributes = (;), kwargs...)
    guide_gap=_addon_guide_gap(guide_gap)
    guide_spacing=_addon_guide_spacing(guide_spacing)
    _addon_colorbar_group_attributes(colorbar_group_attributes)
    axis_keys = (propertynames(Axis)..., :palette)
    axis_attributes = merge(
        (;
            (key=>value for (key, value) in kwargs if key in axis_keys)...), axis)
    series_attributes = (; (key=>value for (key, value) in kwargs if key ∉ axis_keys)...)
    size isa Tuple{Int, Int} && all(>(0), size) ||
        throw(ArgumentError("figure size must be a tuple of two positive integers"))
    figure = Figure(; merge((; size, figure_padding = (12, 12, 12, 12)), figure)...)
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
    chrome=Any[]
    if controls
        root[2, 1] = body
        root[1, 1] = toolbar
        status_label=Label(root[3, 1], status; halign = :left, fontsize = 11)
        append!(chrome, (toolbar, status_label))
        rowsize!(root, 1, Auto(true))
        rowsize!(root, 2, Auto(false, 1))
        rowsize!(root, 3, Auto(true))
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
    return (; figure, reference_size = Tuple(figure.scene.viewport[].widths),
        root, body, canvas, toolbar, status, chrome, axis_attributes,
        series_attributes, widgets, guide_gap, guide_spacing,
        colorbar_position, colorbar_attributes, colorbar_group_attributes)
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

function _addon_scale(symbol::Symbol)
    symbol === :linear && return Base.identity
    symbol === :log10 && return Base.log10
    # Same signed-log scale, with no cancellation in its linear neighbourhood.
    symbol === :pseudolog10 && return Makie.ReversibleScale(
        x -> sign(x) * log1p(abs(x)) / log(10),
        x -> sign(x) * expm1(abs(x) * log(10));
        limits = (0.0f0, 3.0f0), name = :pseudolog10)
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
        isapprox(vmin, vmax; rtol = sqrt(eps(Float64)), atol = 0) &&
            return unique([vmin, vmax])
        ticks = Float64[]
        for exponent in floor(Int, log10(vmin)):floor(Int, log10(vmax))
            lower, upper = max(vmin, 10.0^exponent), min(vmax, 10.0^(exponent+1))
            lower < upper || continue
            budget = max(2, ceil(Int, count*(log10(upper)-log10(lower))/span)+1)
            # Locator arithmetic is local to one decade, including extreme SI
            # magnitudes. Only the returned positions use the published units.
            factor = 10.0^clamp(exponent, -307, 307)
            values = Makie.get_tickvalues(Makie.LinearTicks(budget), lower/factor, upper/factor)
            append!(ticks, filter(x -> isfinite(x) && lower <= x <= upper, values .*
                                                                           factor))
            vmin <= 10.0^exponent <= vmax && push!(ticks, 10.0^exponent)
        end
        isempty(ticks) && append!(ticks, (vmin, vmax))
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

function _addon_set_axis!(entries::AbstractVector, dim::Symbol, scale = nothing)
    dim in (:x, :y) || throw(ArgumentError("axis dimension must be :x or :y"))
    index = dim === :x ? 1 : 2
    # Resolve and validate the complete page before any native observable changes.
    # These are native axis bindings, not another interpretation of result data.
    targets = map(entries) do entry
        requested_scale = scale === nothing ? entry.scale : scale
        target = _addon_scale(requested_scale)
        context = "axis :$dim ($(repr(entry.axis.title[])))"
        requested = entry.axis.limits[]
        requested = length(requested) == 4 ? (requested[1:2], requested[3:4]) : requested
        bounds = requested[index] === nothing ? () : requested[index]
        all(value -> value === nothing || isfinite(value), bounds) ||
            throw(DomainError(bounds, "$context requires finite explicit limits"))
        values = _addon_visible_values(entry.axis, dim, entry.series)
        if requested_scale === :log10 &&
           (any(<=(0), values) || any(value -> value!==nothing && value<=0, bounds))
            target = _addon_scale(:pseudolog10)
        end
        if target === Base.log10
            all(>(0), values) && all(value -> value === nothing || value > 0, bounds) ||
                throw(DomainError(bounds,
                    "logarithmic $context requires positive visible data, uncertainty bounds and explicit limits"))
        end
        inverse=Makie.inverse_transform(target)
        inverse===nothing &&
            throw(ArgumentError("$context requires a native scale with an inverse"))
        for value in bounds
            value===nothing && continue
            transformed=target(value)
            isfinite(transformed) && isfinite(inverse(transformed)) || throw(DomainError(
                value, "$context requires finite transformed explicit bounds and inverse values"))
        end
        empty_bounds=isempty(values) ? defaultlimits(requested[index], target) : nothing
        if !isempty(values)
            lower, upper = extrema(values)
            if isapprox(lower, upper; rtol = sqrt(eps(Float64)), atol = 0)
                lower, upper = _addon_constant_limits(values, values, target === Base.log10)
            end
            explicit = isempty(bounds) ? (nothing, nothing) : bounds
            lower, upper = something(explicit[1], lower), something(explicit[2], upper)
            transformed = (target(lower), target(upper))
            all(isfinite, transformed) && transformed[1] < transformed[2] ||
                throw(DomainError((lower, upper), "$context requires distinct finite transformed limits"))
        end
        if empty_bounds!==nothing
            transformed=target.(empty_bounds)
            all(isfinite, transformed) && transformed[1]<transformed[2] &&
            all(isfinite, inverse.(transformed)) || throw(DomainError(empty_bounds,
                "$context requires finite distinct native default limits"))
        end
        (; scale = target, empty_bounds)
    end
    previous=[(scale = getproperty(entry.axis, Symbol(dim, :scale))[],
                  view = entry.axis.targetlimits[])
              for entry in entries]
    try
        for (entry, target, saved) in zip(entries, targets, previous)
            axis=entry.axis
            if target.empty_bounds!==nothing
                # Native empty-data fitting keeps the previous target limits.
                # Seed its own scale defaults through a linear transition so an
                # empty log axis never inherits the linear interval [0,10].
                getproperty(axis, Symbol(dim, :scale))[]=identity
                current=axis.targetlimits[]
                origin, widths=collect(current.origin), collect(current.widths)
                origin[index]=target.empty_bounds[1]
                widths[index]=target.empty_bounds[2]-target.empty_bounds[1]
                axis.targetlimits[]=Makie.Rect2d(origin..., widths...)
            end
            getproperty(axis, Symbol(dim, :scale))[]=target.scale
            entry.reset(; xauto = dim===:x, yauto = dim===:y)
            # Scale application preserves the orthogonal interactive view.
            other=3-index
            current=axis.targetlimits[]
            origin, widths=collect(current.origin), collect(current.widths)
            origin[other], widths[other]=saved.view.origin[other], saved.view.widths[other]
            axis.targetlimits[]=Makie.Rect2d(origin..., widths...)
        end
    catch
        for (entry, saved) in zip(entries, previous)
            getproperty(entry.axis, Symbol(dim, :scale))[]=identity
            entry.axis.targetlimits[]=saved.view
            getproperty(entry.axis, Symbol(dim, :scale))[]=saved.scale
            entry.axis.targetlimits[]=saved.view
        end
        rethrow()
    end
    return entries
end

function _addon_numeric_values(values)
    # Undefined observations remain missing in retained products. Makie's numeric
    # line boundary uses NaN gaps, including an entirely undefined phase trace.
    nominal_values = map(value -> ismissing(value) ? NaN : LineCableModels.nominal(value), values)
    errors = LineCableModels.uncertainty.(values)
    return nominal_values, any(error -> !iszero(error), errors) ? errors : nothing
end

function _addon_line!(axis, xdata, ydata; dependent_plots, label, color = nothing,
        visible = true, phase = (1, 1), endpoints = false, marker_coordinates = nothing,
        errorbar_sampling = :all, yerror = nothing, interval_support = nothing)
    x, xerror = _addon_numeric_values(xdata)
    y, inferred_error = _addon_numeric_values(ydata)
    yerror = yerror === nothing ? inferred_error :
             (all(iszero, yerror) ? nothing : Float64.(yerror))
    attributes = color === nothing ? (; linewidth = 2) : (; linewidth = 2, color)
    plots = Any[lines!(axis, x, y; label, visible, attributes...)]
    line = first(plots)
    uncertain_indices = findall(eachindex(x)) do index
        isfinite(x[index]) && isfinite(y[index]) &&
            any(
                error -> error !== nothing && isfinite(error[index]) &&
                         !iszero(error[index]), (xerror, yerror))
    end
    glyphs = if marker_coordinates !== nothing || errorbar_sampling === :staggered
        lift(line[1], axis.scene.viewport) do points, viewport
            _addon_glyph_indices(length(points), viewport.widths[1], phase;
                endpoints, uncertain_indices = filter(<=(length(points)), uncertain_indices),
                errorbar_sampling)
        end
    else
        nothing
    end
    if marker_coordinates !== nothing
        marker_coordinates[line] = lift(line[1], glyphs) do points, indices
            points[indices.markers]
        end
    end
    error_color = color === nothing ? :black : color
    for (direction, error) in ((:y, yerror), (:x, xerror))
        error === nothing && continue
        values = only(Makie.convert_arguments(Makie.Errorbars, x, y, error))
        coordinates = if errorbar_sampling === :staggered
            # A native line edit does not rewrite independent uncertainty bars.
            lift(axis.scene.viewport) do viewport
                indices=_addon_glyph_indices(length(values), viewport.widths[1], phase;
                    endpoints, uncertain_indices, errorbar_sampling)
                values[indices.intervals]
            end
        else
            values
        end
        bars = errorbars!(axis, coordinates; color = error_color, direction,
            whiskerwidth = 3, linewidth = 1, visible)
        interval_support===nothing || (interval_support[bars]=coordinates)
        on(axis.scene, line.color) do value
            bars.color[] = value
        end
        push!(plots, bars)
    end
    append!(dependent_plots, (plot => first(plots) for plot in Iterators.drop(plots, 1)))
    return plots
end

function _addon_visible_values(series, dim::Symbol; include_uncertainty::Bool = false)
    values = Float64[]
    for item in series
        first(item.plots).visible[] || continue
        data = dim === :x ? item.xdata : item.ydata
        data === nothing && continue
        if !include_uncertainty
            bounds=Makie.data_limits(first(item.plots))
            index=dim===:x ? 1 : 2
            lower=bounds.origin[index]
            upper=lower+bounds.widths[index]
            isfinite(lower) && isfinite(upper) && append!(values, (lower, upper))
            continue
        elseif haskey(item, :full_support)
            append!(values, getproperty(item.full_support, dim))
            continue
        end
        errors = get(item, dim === :x ? :xerror : :yerror, nothing)
        for (index, sample) in enumerate(data)
            nominal_value = LineCableModels.nominal(sample)
            nominal_value isa Real || continue
            numeric = Float64(nominal_value)
            isfinite(numeric) || continue
            interval = abs(Float64(errors === nothing ?
                                   LineCableModels.uncertainty(sample) : errors[index]))
            if isfinite(interval) && !iszero(interval)
                push!(values, numeric - interval, numeric + interval)
            else
                push!(values, numeric)
            end
        end
    end
    return values
end

# Primary categorical and element-index points are all data samples, never
# decorative curve glyphs. Native conversion owns categorical coordinates.
function _addon_points!(axis, xdata, ydata; dependent_plots, label, color = nothing,
        visible = true, phase = (1, 1), endpoints = false, marker_coordinates = nothing,
        errorbar_sampling = :all, yerror = nothing, interval_support = nothing)
    y, inferred=_addon_numeric_values(ydata)
    errors=yerror===nothing ? inferred : yerror
    attributes=color===nothing ? (;) : (; color)
    points=scatter!(axis, xdata, y; label, visible, attributes...)
    plots=Any[points]
    if errors!==nothing && !all(iszero, errors)
        bars=errorbars!(axis, xdata, y, errors; visible, color = something(color, :black))
        on(axis.scene, points.color) do value
            bars.color[]=value
        end
        push!(plots, bars)
        push!(dependent_plots, bars=>points)
    end
    return plots
end

function _addon_visible_values(axis::Axis, dim::Symbol, series = ())
    index = dim === :x ? 1 : 2
    bounds = Makie.data_limits(axis.scene,
        plot -> !to_value(get(plot, :visible, true)) ||
                !to_value(get(plot, Symbol(dim, :autolimits), true)) ||
                to_value(get(plot, :space, :data)) !== :data)
    lower, upper = bounds.origin[index], bounds.origin[index] + bounds.widths[index]
    values = isfinite(lower) && isfinite(upper) ? [lower, upper] : Float64[]
    for item in series
        get(item, :sampled_intervals, false) || continue
        # Undrawn intervals still constrain scientific axes. Independently hidden
        # bars, disabled autolimits and caller-added native plots remain respected.
        support=get(item, :interval_support, nothing)
        any(
            plot -> plot isa Makie.Errorbars && plot.direction[] === dim &&
                    plot.visible[] && to_value(get(plot, Symbol(dim, :autolimits), true)) &&
                    (support===nothing || !haskey(support, plot) ||
                     isequal(plot[1][], to_value(support[plot]))),
            item.plots) || continue
        append!(values, _addon_visible_values((item,), dim; include_uncertainty = true))
    end
    return isempty(values) ? values : collect(extrema(values))
end

function LineCableModels.plotwindow(
        callback::F;
        title::AbstractString,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        series_attributes = nothing,
        size::Tuple{Int, Int} = (800, 400),
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
        callback(shell.canvas)
        axes = Any[content for content in shell.figure.content if content isa Axis]
        # Only caller-supplied overrides apply to caller-constructed axes. Keep
        # native defaults/conversions intact; scale changes use common preflight.
        for axis in axes, (key, value) in pairs(shell.axis_attributes)

            key in (:xscale, :yscale) || setproperty!(axis, key, value)
        end
        requested_scales = !any(key -> key in (:xscale, :yscale), keys(shell.axis_attributes)) ?
                           nothing :
                           [(x = get(shell.axis_attributes, :xscale, axis.xscale[]),
                                y = get(shell.axis_attributes, :yscale, axis.yscale[]))
                            for axis in axes]
        resets = Function[_addon_reset!(axis) for axis in axes]
        native = Any[handle for axis in axes for handle in axis.scene.plots]
        order = [Symbol("series_$index") for index in eachindex(native)]
        groups = Dict(group => Any[handle] for (group, handle) in zip(order, native))
        _addon_finish!(
            shell,
            axes,
            resets,
            groups,
            order,
            Dict(group => something(to_value(get(handle, :label, nothing)), string(group))
            for (group, handle) in zip(order, native));
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
        axis, scales = _addon_axis!(
            panel.content,
            xobservation,
            yobservation;
            title = panel_title,
            xscale,
            yscale,
            xlabel,
            ylabel,
            native_attributes = shell.axis_attributes
        )
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, Any}()
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
            requested_scales = (scales,),
            axis_series = (series,),
            series_attributes,
            title,
            figure_title,
            title_attributes,
            legend_position,
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
        halfspan = max(log(1.05), 2maximum(value -> abs(log(value)-center), interval_values))
        lower, upper = exp(center-halfspan), exp(center+halfspan)
        isfinite(upper) && 0 < lower < upper && log10(lower) < log10(upper) ||
            throw(DomainError((lower, upper), "automatic logarithmic bounds are not representable"))
        return lower, upper
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
        previous_inputs = Ref{Any}(nothing)
        updating = Ref(false)
        # Reuse one native text measurement per dimension, outside the data
        # scene. Font changes and long, precise labels affect available density.
        probe = text!(axis.blockscene, 0, 0; text = "", markerspace = :data,
            visible = false, inspectable = false)
        function update!()
            updating[] && return nothing
            updating[] = true
            try
                limits = axis.finallimits[]
                current_scale, current_ticks = scale[], ticks[]
                current_format, current_label = tickformat[], label[]
                # A transform notification precedes Makie's own limit reset.
                # Do not send its old, possibly negative linear view to a log locator.
                current_scale === Base.log10 && limits.origin[index] <= 0 && return nothing
                label_changed = current_label !== rendered_label[]
                label_changed && (raw_label[] = current_label)
                owned_ticks = current_ticks === installed_ticks[] ||
                              current_ticks === Makie.automatic
                owned_format = current_format === installed_format[] ||
                               current_format === Makie.automatic
                numeric = conversion[] === nothing
                lower, upper = limits.origin[index],
                limits.origin[index] + limits.widths[index]
                inputs=(lower, upper, current_scale, axis.scene.viewport[].widths[index],
                    labelsize[], labelfont[], rotation[], conversion[], current_label)
                # Native layout can notify every peer axis without changing
                # these inputs. Restarting density fitting then briefly installs
                # denser ticks and recursively changes all peer protrusions.
                # Explicit native overrides/automatic resets still run below.
                owned_ticks && owned_format && !label_changed &&
                    current_ticks===installed_ticks[] &&
                    current_format===installed_format[] &&
                    isequal(inputs, previous_inputs[]) && return nothing
                decades = current_scale === Base.log10 && 0 < lower < upper &&
                          log10(upper) - log10(lower) >= 2
                exponent = something(_addon_scientific_exponent((lower, upper)), 0)
                signed_linear = current_scale === _addon_scale(:pseudolog10) &&
                                max(abs(lower), abs(upper)) < 1
                mode = if numeric && (current_scale === Base.identity || signed_linear) &&
                          (owned_ticks || current_ticks isa AbstractVector{<:Real})
                    (:linear, exponent)
                elseif numeric && current_scale === Base.log10 &&
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
                            values -> [value <= 0 ? string(value) :
                                       Makie.rich("10",
                                           Makie.superscript(
                                               replace(string(round(Int, log10(value))), "-" => "−");
                                               offset = Makie.Vec2f(0.1, 0.0)))
                                       for value in values]
                        end
                        installed_mode[] = mode
                        tickformat[] === installed_format[] ||
                            (tickformat[] = installed_format[])
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
                    spacing = (index == 1 && current_scale === Base.identity ? 5.5 : 3.0) *
                              labelsize[]
                    count = clamp(floor(Int, pixels / max(1, spacing)), 3, 10)
                    probe.font[] = labelfont[]
                    probe.fontsize[] = labelsize[]
                    probe.rotation[] = rotation[]
                    # Native locator/formatter notification is synchronous. Fit
                    # its rendered strings, reducing only automatic tick density.
                    for _ in 1:9
                        selected = if !numeric
                            Makie.automatic
                        elseif current_scale === Base.identity
                            # A native fit can briefly expose equal or adjacent
                            # Float64 endpoints before automatic padding runs.
                            # Preserve explicit narrow zooms too, without asking
                            # LinearTicks to subdivide an unrepresentable interval.
                            let count=count
                                (lo, hi) -> isapprox(lo, hi; rtol = sqrt(eps(Float64)), atol = 0) ?
                                            unique([lo, hi]) :
                                            Makie.get_tickvalues(Makie.LinearTicks(count), lo, hi)
                            end
                        elseif current_scale === Base.log10
                            _addon_decade_ticks(lower, upper, count)
                        elseif current_scale === _addon_scale(:pseudolog10)
                            # Native PseudologTicks dispatches on Makie's scale
                            # instance. Reuse its placement, not its cancelling
                            # transform, and pass numeric positions to this axis.
                            if isapprox(lower, upper; rtol = sqrt(eps(Float64)), atol = 0)
                                unique([lower, upper])
                            else
                                locator = signed_linear ? Makie.LinearTicks(count) :
                                          Makie.PseudologTicks(count)
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
                        if current_scale === Base.log10 && !decades
                            # Physical-value log ticks are NOT equally spaced on
                            # screen. Fit adjacent rendered labels, not count/width.
                            positions = Float64[]
                            extents = Float64[]
                            for (value, text) in zip(selected, lineaxis.ticklabels[])
                                value > 0 || continue
                                probe.text[] = text
                                push!(positions, pixels*(log10(value)-log10(lower))/(log10(upper)-log10(lower)))
                                push!(extents, Makie.boundingbox(probe, :data).widths[index])
                            end
                            if length(positions) == length(selected)
                                retained = Int[]
                                for i in eachindex(positions)
                                    if isempty(retained) ||
                                       positions[i]-positions[last(retained)] >=
                                       (extents[i]+extents[last(retained)])/2+labelsize[]/2
                                        push!(retained, i)
                                    end
                                end
                                fitted_ticks = selected[retained]
                                if !isequal(ticks[], fitted_ticks)
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
                        fitted = max(2, floor(Int, pixels /
                                                   max(spacing, extent + labelsize[])))
                        fitted >= count && break
                        count = fitted
                    end
                end
                previous_inputs[]=inputs
            finally
                updating[] = false
            end
            return nothing
        end
        # Run before native tick conversion: a newly assigned labelled tuple
        # cannot be consumed with the previously installed numeric formatter.
        onany((_...) -> update!(), axis.scene, ticks, tickformat; priority = 1)
        # Range/transform updates must instead follow native LineAxis propagation;
        # changing its formatter while it still has the old limits is unsafe.
        onany((_...) -> update!(), axis.scene, axis.finallimits, scale, label,
            axis.scene.viewport, labelsize, labelfont, rotation, conversion;
            priority = -3, update = true)
    end
    return axis
end

# Bind the native limit lifecycle once and return this axis's reset action.
function _addon_reset!(axis, series = ())
    # Own numeric ticks before the first data fit, including its synchronous
    # native callbacks. Every recipe and caller-owned plotwindow uses this bind.
    _addon_axis_format!(axis)
    fitting = Ref(false)
    corrections = Any[nothing, nothing]
    function reset!(; xauto::Bool = true, yauto::Bool = true, preserve_view::Bool = false)
        fitting[] && return axis
        view = preserve_view ? axis.targetlimits[] : nothing
        fitting[] = true
        try
            reset_limits!(axis; xauto, yauto)
            requested = axis.limits[]
            requested = length(requested) == 4 ? (requested[1:2], requested[3:4]) :
                        requested
            for (index, dim) in enumerate((:x, :y))
                (index == 1 ? xauto : yauto) || continue
                corrections[index] = nothing
                explicit = requested[index] === nothing ? (nothing, nothing) :
                           requested[index]
                all(value -> value !== nothing, explicit) && continue
                rendered = _addon_visible_values(axis, dim)
                interval_values = _addon_visible_values(axis, dim, series)
                isempty(interval_values) && continue
                values = isempty(series) ? interval_values :
                         _addon_visible_values(series, dim)
                isempty(values) && continue
                constant = isapprox(extrema(values)...; rtol = sqrt(eps(Float64)), atol = 0)
                constant || interval_values != rendered || continue
                if constant && !isempty(series)
                    expected = _addon_visible_values(series, dim; include_uncertainty = true)
                    # Native extra plots or independently hidden error bars own
                    # their actual extents, not the original observation array.
                    if !all(isapprox.(interval_values, collect(extrema(expected))))
                        isapprox(extrema(interval_values)...; rtol = sqrt(eps(Float64)), atol = 0) ||
                            continue
                        values = interval_values
                    end
                end
                any(
                    plot -> haskey(plot, :model) &&
                            any(j -> plot.model[][index, j] != (index == j), 1:4),
                    axis.scene.plots) && continue
                scale = getproperty(axis, Symbol(dim, :scale))[]
                limits = if constant
                    _addon_constant_limits(values, interval_values, scale === Base.log10)
                else
                    lower, upper = scale.(extrema(interval_values))
                    low_margin, high_margin = getproperty(axis, Symbol(dim, :autolimitmargin))[]
                    span = upper - lower
                    Makie.inverse_transform(scale).((lower - low_margin * span,
                        upper + high_margin * span))
                end
                lower = something(explicit[1], limits[1])
                upper = something(explicit[2], limits[2])
                origin, widths = collect(axis.targetlimits[].origin),
                collect(axis.targetlimits[].widths)
                original = (origin[index], widths[index])
                origin[index], widths[index] = lower, upper - lower
                corrections[index] = (; requested = requested[index], scale, original,
                    fitted = (origin[index], widths[index]))
                axis.targetlimits[] = Makie.Rect2d(origin..., widths...)
            end
        finally
            try
                view === nothing || (axis.targetlimits[] = view)
            finally
                fitting[] = false
            end
        end
        return axis
    end
    # Makie also refits before displaying a Figure. Preserve our degenerate-data
    # padding without storing automatic limits as user requests. Only the exact
    # native auto-fit for the same scale/request is corrected; zooms are untouched.
    # This callback uses two cached bounds, not a data scan on every view update.
    on(axis.scene, axis.targetlimits; priority = 1) do view
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
                    getproperty(axis, Symbol(dim, :scale))[] === Base.log10)
                lower, upper = something(explicit[1], bounds[1]),
                something(explicit[2], bounds[2])
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
    on(_ -> reset!(), axis.scene, axis.limits; priority = -3)
    if any(
        item -> get(item, :sampled_intervals, false) &&
                any(plot -> plot isa Makie.Errorbars, item.plots),
        series)
        # A resize changes the drawn interval subset. Refresh the cached native
        # auto-fit after glyph/layout updates, without replacing the user's view.
        # A subsequent native display fit then still restores complete bounds.
        on(_ -> reset!(; preserve_view = true), axis.scene, axis.scene.viewport; priority = -4)
    end
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
                  Units.label(xobservation.quantity, xobservation.unit) :
                  xlabel
    yaxis_label = ylabel === nothing ?
                  Units.label(yobservation.quantity, yobservation.unit) :
                  ylabel
    options = merge(
        (;
            title,
            tellwidth = false,
            tellheight = false,
            xscale,
            yscale,
            xlabel = xaxis_label,
            ylabel = yaxis_label
        ),
        attributes, native_attributes)
    scales = (x = options.xscale, y = options.yscale)
    # Draw on safe axes, then validate complete native extents in the common
    # finish before applying requested transforms to any axis on the page.
    axis = Axis(position; merge(options, (xscale = identity, yscale = identity))...)
    return axis, scales
end

function _addon_panel_titles(panel_titles, expected::Int; defaults = nothing)
    panel_titles === nothing && return nothing
    panel_titles isa Function && return Tuple(panel_titles(i) for i in 1:expected)
    panel_titles isa AbstractDict &&
        return Tuple(get(panel_titles, i, defaults===nothing ? nothing : defaults[i])
        for i in 1:expected)
    panel_titles isa Tuple || panel_titles isa AbstractVector ||
        throw(ArgumentError(
            "panel_titles must be a tuple, vector, dictionary, function, or nothing",
        ))
    length(panel_titles) == expected || throw(DimensionMismatch(
        "panel_titles must contain one entry per logical plot panel",
    ))
    return Tuple(panel_titles)
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

function _addon_remove_legend!(legend)
    legend === nothing && return nothing
    content = GridLayoutBase.gridcontent(legend)
    layout = content === nothing ? nothing : content.parent
    delete!(legend)
    if layout !== nothing && isempty(layout.content)
        layout_content = GridLayoutBase.gridcontent(layout)
        layout_content === nothing ||
            GridLayoutBase.remove_from_gridlayout!(layout_content)
    end
    return nothing
end

function _addon_legend_sources!(legend, groups, dependents)
    legend === nothing && return nothing
    owners = IdDict{Any, Any}(dependents)
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
            # Makie's LegendElement extension contract requires this mutable
            # vector to identify the plots represented by a glyph.
            targets = element.plots
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
    visibilities = map(collect(keys(owners))) do plot
        # on returns an ObserverFunction with a documented `observable` field.
        # Release the temporary subscription; native legend listeners own the
        # updates. notify(plot.visible) is a no-op for Makie Computed inputs.
        subscription = on(identity, plot.visible)
        observable = subscription.observable
        off(subscription)
        observable
    end
    on(legend.blockscene, legend.entrygroups; priority = -1) do _
        foreach(notify, visibilities)
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
    legend.entrygroups[] = [(first(only(legend.entrygroups[])), displayed)]
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
    on(legend.blockscene, bounding_box) do bounds
        fit!(bounds)
    end
    on(legend.blockscene, legend.orientation) do _
        empty!(extents)
        fit!(bounding_box[])
    end
    fit!(bounding_box[])
    return legend
end

function _addon_wrap_legend_label(label, width, measure)
    lines = String[]
    for paragraph in split(label, '\n'; keepempty = true)
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
function _addon_grid_legend!(bounds, legend, position; automatic = true, fitting_geometry)
    entries = last(only(legend.entrygroups[]))
    originals = Any[entry.label[] for entry in entries]
    fitting = Ref(false)
    managed = Ref(automatic)
    previous_metrics = Ref{Any}(nothing)
    # Match native Label: glyphs and positions must share data space for the
    # bounding box to include font extents rather than just the anchor point.
    probe = text!(legend.blockscene, 0, 0; text = "", markerspace = :data,
        visible = false, inspectable = false)
    function fit!()
        (fitting[] || !managed[]) && return nothing
        horizontal = position[] in (:top, :bottom)
        available = Float64(bounds[].widths[1]) - sum(legend.margin[][1:2]) -
                    sum(legend.padding[][1:2]) - 4
        available > 0 || return nothing
        fitting[] = true
        previous_geometry=fitting_geometry[]
        fitting_geometry[]=true
        try
            widths = Float64[]
            for (entry, original) in zip(entries, originals)
                probe.font[] = entry.labelfont[]
                probe.fontsize[] = entry.labelsize[]
                measured = Dict{Any, Float64}()
                measure = text -> get!(measured, text) do
                    probe.text[] = text
                    Float64(Makie.boundingbox(probe, :data).widths[1])
                end
                patch = Float64(entry.patchsize[][1]) + legend.patchlabelgap[]
                wrapped = !horizontal || !(original isa String) ||
                          measure(original) <= available-patch ? original :
                          _addon_wrap_legend_label(original, max(1.0, available-patch), measure)
                entry.label[] == wrapped || (entry.label[] = wrapped)
                push!(widths, patch + measure(wrapped))
            end
            columns = 1
            for count in (horizontal ? (length(entries):-1:1) : (1:-1:1))
                required = sum(maximum(widths[column:count:end]) for column in 1:count) +
                           (count-1)*legend.colgap[]
                if required <= available
                    columns = count
                    break
                end
            end
            metrics = (Tuple(widths), Tuple(entry.label[] for entry in entries),
                Tuple((entry.labelsize[], entry.labelfont[]) for entry in entries))
            legend.orientation[]===:vertical || (legend.orientation[]=:vertical)
            if legend.nbanks[] != columns
                legend.nbanks[] = columns
            elseif previous_metrics[] != metrics
                # Native Legend relayout is triggered by bank/layout controls,
                # not by changed label extents alone.
                notify(legend.nbanks)
            end
            previous_metrics[] = metrics
        finally
            fitting_geometry[]=previous_geometry
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
            originals[index] = label
            fit!()
        end
    end
    for setting in (legend.nbanks, legend.orientation)
        on(legend.blockscene, setting; priority = 1) do _
            fitting[] || (managed[]=false)
        end
    end
    fit!()
    return fit!
end

function _addon_axes_viewport(axes, fallback; panels = (), cells = ())
    isempty(axes) && return fallback
    viewports=Tuple(axis.scene.viewport for axis in axes)
    if isempty(cells)
        return lift(viewports...) do bounds...
            left=minimum(bound.origin[1] for bound in bounds)
            bottom=minimum(bound.origin[2] for bound in bounds)
            right=maximum(bound.origin[1]+bound.widths[1] for bound in bounds)
            top=maximum(bound.origin[2]+bound.widths[2] for bound in bounds)
            Rect2f((left, bottom), (right-left, top-bottom))
        end
    end
    panel_boxes=Tuple(panel.layout.layoutobservables.computedbbox for panel in panels)
    cell_boxes=Tuple(cell.layout.layoutobservables.computedbbox for cell in cells)
    return lift(viewports..., panel_boxes..., cell_boxes...) do values...
        n=length(axes)
        frames=values[1:n]
        selected=values[(n + 1):2n]
        domain=values[(2n + 1):end]
        lower=ntuple(
            d -> minimum(frame.origin[d]-panel.origin[d]
            for (frame, panel) in zip(frames, selected)),
            2)
        upper=ntuple(
            d -> minimum(panel.origin[d]+panel.widths[d]-frame.origin[d]-frame.widths[d]
            for (frame, panel) in zip(frames, selected)),
            2)
        origin=ntuple(d -> minimum(box.origin[d] for box in domain)+lower[d], 2)
        stop=ntuple(d -> maximum(box.origin[d]+box.widths[d] for box in domain)-upper[d], 2)
        Rect2f(origin, stop .- origin)
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
            labels[group] = replacement
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
        labels[group] = replacement
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
    displayed = Any[]
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
                halign = :right, valign = :top,
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
    target===nothing &&
        throw(ArgumentError("native legends require an assigned guide slot"))
    default_orientation=target_orientation
    options = merge(
        (;
            orientation = default_orientation,
            halign = default_orientation === :vertical ? :left : :center,
            valign = default_orientation === :vertical ? :top : :center,
            tellwidth = position in (:left, :right) || position isa Tuple,
            tellheight = position in (:top, :bottom) || position isa Tuple
        ),
        attributes)
    if overflow === :show_all
        legend = Legend(target, entries, displayed, title; options...)
        return _addon_legend_sources!(legend, groups, dependent_plots)
    end
    ellipsis = LineElement(color = :transparent)
    legend = Legend(
        target,
        Any[entries..., ellipsis],
        [displayed; "(...)"],
        title;
        options...
    )
    _addon_legend_sources!(legend, groups, dependent_plots)
    return _addon_responsive_legend!(
        figure,
        inside_bbox,
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
    result = Dict{Any, Any}()
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

function _addon_panel_identity(value)
    value isa Integer && !(value isa Bool) && value>0 && return Int(value)
    value isa Tuple && length(value) == 2 &&
    all(index -> index isa Integer && !(index isa Bool) && index > 0, value) ||
        throw(ArgumentError("panel identities must be positive indices or original `(row, column)` tuples"))
    return (Int(value[1]), Int(value[2]))
end

function _addon_without_legend_controls(options::NamedTuple)
    names = Tuple(filter(
        name -> name ∉ (:position, :overflow, :title, :legend_labels),
        keys(options)
    ))
    return NamedTuple{names}(Tuple(getproperty(options, name) for name in names))
end

function _addon_legend_configuration(value; default_position, default_title = nothing)
    value isa Symbol && return (;
        position = value,
        overflow = :show_all,
        title = default_title,
        legend_labels = nothing,
        attributes = (;)
    )
    value isa NamedTuple || throw(ArgumentError(
        "a legend configuration must be a dock symbol or NamedTuple",
    ))
    position = get(value, :position, default_position)
    overflow = get(value, :overflow, :show_all)
    title = get(value, :title, default_title)
    haskey(value, :anchor) &&
        throw(ArgumentError("anchor was removed; use native halign and valign"))
    legend_labels = get(value, :legend_labels, nothing)
    attributes = _addon_without_legend_controls(value)
    return (; position, overflow, title, legend_labels, attributes)
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
        caption = colorbar.axis.elements[:labeltext]
        managed=Ref{Any}(colorbar.alignmode[])
        onany(colorbar.blockscene, fast_string_boundingboxes_obs(labels),
            fast_string_boundingboxes_obs(caption), colorbar.vertical, colorbar.ticklabelsvisible,
            colorbar.labelvisible, colorbar.layoutobservables.computedbbox, colorbar.spinewidth,
            colorbar.layoutobservables.protrusions; update = true) do boxes,
        captions, vertical, visible, labelvisible, bounds, stroke, protrusions
            colorbar.alignmode[]==managed[] || return nothing
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
            if labelvisible
                # The property caption is centered along the bar. Its overhang
                # belongs to this item's footprint, just like endpoint ticks.
                half_length=bounds.widths[dimension]/2
                for box in captions
                    isfinite(box.origin[dimension]) && isfinite(box.widths[dimension]) ||
                        continue
                    before=max(before, ceil(-box.origin[dimension]-half_length))
                    after=max(after, ceil(box.origin[dimension]+box.widths[dimension]-half_length))
                end
            end
            border=max(0.0, Float64(stroke)/2)
            before=max(before, border)
            after=max(after, border)
            # Enclose the visible border as well as native perpendicular text.
            # These are measured protrusions, not sibling-spacing margins.
            # Native LineAxis places text beyond a full spine width, while
            # its reported text protrusion omits that offset.
            side(name) = GridLayoutBase.Protrusion(max(
                getproperty(protrusions, name) +
                (getproperty(protrusions, name)>0 ? max(0.0, stroke) : 0.0),
                border))
            managed[] = vertical ?
                        Mixed(bottom = before, top = after, left = side(:left), right = side(:right)) :
                        Mixed(left = before, right = after, bottom = side(:bottom), top = side(:top))
            colorbar.alignmode[]==managed[] || (colorbar.alignmode[]=managed[])
        end
    end
    return colorbar
end

function LineCableModels.materialscale!(position, scheme; kwargs...)
    return _addon_colorbar!(position, scheme; attributes = (; kwargs...))
end

function _addon_colorbars!(slot, scales; attributes, orientation, main = false)
    attributes isa NamedTuple ||
        throw(ArgumentError("colorbar_attributes must be a NamedTuple"))
    vertical=get(attributes, :vertical, orientation===:vertical)
    length=main ? Auto() : _ADDON_COLORBAR_DOCK_LENGTH
    defaults=vertical ? (; vertical, height = length) : (; vertical, width = length)
    grid=GridLayout(; alignmode = Outside())
    slot[]=grid
    scene=GridLayoutBase.top_parent(grid).scene
    previous_scenes=copy(scene.children)
    items=try
        map(enumerate(scales)) do (index, scale)
            item=GridLayout(; alignmode = Outside())
            grid[index, 1]=item
            bar=_addon_colorbar!(item[1, 1], scale; attributes = merge(defaults, attributes))
            companion=Label(item[1, 2], bar.label; halign = :right, valign = :center,
                fontsize = bar.labelsize, font = bar.labelfont, color = bar.labelcolor)
            (; layout = item, bar, companion, visible = Ref(bar.blockscene.visible[]),
                labelvisible = Ref(bar.labelvisible[]),
                managed = (vertical = Ref(!haskey(attributes, :vertical)),
                    width = Ref(!haskey(attributes, :width)), height = Ref(!haskey(attributes, :height))))
        end
    catch
        _addon_delete_subtree!(grid)
        # A native constructor can fail before registering its block in the
        # grid. Release that incomplete scene and its owned subscriptions too.
        for child in copy(scene.children)
            any(previous -> previous===child, previous_scenes) || empty!(child)
        end
        rethrow()
    end
    return (; colorbars = Tuple(item.bar for item in items), layout = grid, items)
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

function _addon_figure_title!(shell, title, attributes)
    title === nothing && return nothing
    attributes isa NamedTuple || throw(ArgumentError(
        "title_attributes must be a NamedTuple",
    ))
    options = merge(
        (; tellwidth = true, halign = :center, font = :bold, fontsize = 18),
        attributes
    )
    return Label(shell.root[0, 1], title; options...)
end

function _addon_finish!(
        shell,
        axes,
        resets,
        groups,
        order,
        group_labels;
        requested_scales = nothing,
        axis_series = nothing,
        dependent_plots = Pair{Makie.Plot, Makie.Plot}[],
        title,
        figure_title = nothing,
        title_attributes = (;),
        series_attributes = nothing,
        series_defaults = nothing,
        marker_coordinates = nothing,
        legend_position,
        legend_attributes,
        legend_overflow = :ellipsis,
        legend_title = nothing,
        panels = (),
        frame_cells = (),
        panel_legends = (),
        panel_group_labels = nothing,
        panel_legend_titles = nothing,
        color_scales = (),
        colorbar_position = nothing,
        colorbar_attributes = (;),
        controls,
        display_plot,
        export_name,
        export_theme,
        open_export
)
    _addon_colorbar_group_attributes(shell.colorbar_group_attributes, length(color_scales))
    _addon_validate_native_guide(Colorbar, merge(colorbar_attributes, shell.colorbar_attributes))
    isempty(axes) && !isempty(shell.axis_attributes) &&
        throw(ArgumentError(
            "native Axis attributes require a figure containing an Axis"))
    append!(dependent_plots,
        _addon_series_styles!(groups, order, series_attributes;
            defaults = series_defaults, shared = shell.series_attributes, marker_coordinates))
    entries = [(; axis, reset, series = axis_series === nothing ? () : axis_series[i])
               for (i, (axis, reset)) in enumerate(zip(axes, resets))]
    setters = map(enumerate((:x, :y))) do (index, dim)
        if requested_scales !== nothing
            _addon_set_axis!(
                [merge(entry, (scale = getproperty(scales, dim),))
                 for (entry, scales) in zip(entries, requested_scales)],
                dim)
        end
        filter(entries) do entry
            _addon_numeric_axis(entry.axis, dim) &&
                getproperty(entry.axis, Symbol(dim, :scale))[] in
                (identity, log10, _addon_scale(:pseudolog10)) &&
                !isempty(_addon_visible_values(entry.axis, dim, entry.series))
        end
    end
    xsetters, ysetters = setters
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
        shell.canvas.layoutobservables.computedbbox;
        panels, cells = frame_cells
    )
    panel_data = if isempty(panels)
        Dict{Any, Any}(i => begin
                           selected=Dict{Any, Vector{Any}}(key=>filter(
                                                               plot -> plot in
                                                                       axis.scene.plots,
                                                               plots)
                           for (key, plots) in groups
                           if any(plot -> plot in axis.scene.plots, plots))
                           (; axis, panel = (logical_position = i, layout = nothing),
                               groups = selected,
                               order = filter(key -> haskey(selected, key), order), labels = copy(group_labels), title = Ref{Any}(nothing))
                       end for (i, axis) in enumerate(axes))
    else
        _addon_panel_legend_data(panels, axes, groups, order, group_labels;
            panel_labels = panel_group_labels, panel_titles = panel_legend_titles)
    end
    _addon_bind_visibility!(shell.figure, axes, resets, groups, shell.status)
    built = LineCableModels.UIPlot(
        shell.figure,
        Tuple(axes);
        title = title_block,
        status = shell.status,
        addon_state = (;
            shell,
            axis_bindings = entries, controls_enabled = controls,
            widget_bindings = Dict{Symbol, Any}(), widget_order = Symbol[],
            groups,
            dependent_plots,
            order,
            labels = group_labels,
            panel_data,
            inside_bbox,
            guides = Dict{Any, Any}(), guide_order = Any[], guide_docks = Any[],
            composing_guides = Ref(false), guide_gap = shell.guide_gap, guide_spacing = Ref(shell.guide_spacing),
            fitting_geometry = Ref(false), frame_padding = IdDict{Any, Any}(), presentation_ready = Ref(false),
            panel_page = isempty(panels) ? nothing :
                         (index = (1, 1), dimensions = size(shell.canvas),
                coordinates = Tuple(panel.logical_position for panel in panels)),
            color_scales,
            colorbar_group_attributes = Ref{Any}(shell.colorbar_group_attributes)
        ),
        export_name,
        export_theme,
        open_export
    )
    figure_guide=_addon_guide_state(:legend, nothing, legend_position, legend_attributes;
        overflow = legend_overflow, title = legend_title)
    built.addon_state.guides[(:legend, nothing)]=figure_guide
    push!(built.addon_state.guide_order, (:legend, nothing))
    for (identity, value) in _addon_panel_legend_pairs(panel_legends)
        identity=_addon_panel_identity(identity)
        haskey(panel_data, identity) ||
            throw(ArgumentError("panel $(repr(identity)) is absent from this figure"))
        value===nothing || value===false ||
            begin
                config=_addon_legend_configuration(value; default_position = :right, default_title = nothing)
                _addon_relabel_legend!(
                    panel_data[identity].labels, panel_data[identity].groups,
                    panel_data[identity].order, config.legend_labels)
                key=(:legend, identity)
                built.addon_state.guides[key]=_addon_guide_state(
                    :legend, identity, config.position, config.attributes;
                    overflow = config.overflow, title = config.title)
                push!(built.addon_state.guide_order, key)
            end
    end
    resolved_position=shell.colorbar_position===_omitted ? colorbar_position :
                      shell.colorbar_position
    bars=_addon_guide_state(:colorbars, nothing, resolved_position,
        merge(colorbar_attributes, shell.colorbar_attributes))
    built.addon_state.guides[(:colorbars, nothing)]=bars
    push!(built.addon_state.guide_order, (:colorbars, nothing))
    _addon_compose_guides!(built)
    on(shell.figure.scene, shell.figure.scene.viewport) do _
        if !built.addon_state.fitting_geometry[]
            _addon_release_frames!(built)
            _addon_fit_panel_aspects!(built)
        end
        _addon_compose_guides!(built)
        nothing
    end
    _addon_controls!(built, xsetters, ysetters)
    if controls
        for widget in shell.widgets
            widget(built)
        end
    end
    _addon_fit_panel_aspects!(built)
    _addon_compose_guides!(built)
    for axis in axes
        _addon_watch_presentation!(built,
            axis,
            (:title, :titlesize, :titlefont, :subtitle, :subtitlesize,
                :xlabelsize, :ylabelsize, :xticklabelsize, :yticklabelsize))
    end
    title_block===nothing ||
        _addon_watch_presentation!(built, title_block, (:text, :fontsize, :font))
    built.addon_state.presentation_ready[]=true
    _addon_edit_presentation!(() -> nothing, built)
    display_plot && _addon_display!(shell.figure, title)
    return built
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
    return _addon_edit_presentation!(plot) do
        plot.title === nothing || delete!(plot.title)
        plot.title = title === nothing ? nothing :
                     _addon_figure_title!(data.shell, title, (; kwargs...))
        plot.title===nothing ||
            _addon_watch_presentation!(plot, plot.title, (:text, :fontsize, :font))
        _addon_compose_guides!(plot)
        return plot.title
    end
end

function LineCableModels.paneltitle!(
        plot::LineCableModels.UIPlot,
        logical_position,
        title
)
    data = plot.addon_state
    data === nothing && throw(ArgumentError(
        "this plot does not retain logical plot panels",
    ))
    logical_position=_addon_panel_identity(logical_position)
    haskey(data.panel_data, logical_position) ||
        throw(ArgumentError("panel $(repr(logical_position)) is absent from this figure"))
    return _addon_edit_presentation!(plot) do
        axis = data.panel_data[logical_position].axis
        axis.title[] = title === nothing ? "" : title
        return axis
    end
end
