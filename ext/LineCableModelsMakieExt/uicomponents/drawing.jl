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
    errors = uncertainty.(values)
    return nominal_values, any(error -> !iszero(error), errors) ? errors : nothing
end

function _line_errors!(plots, axis, x, y, xerror, yerror, visible, color)
    if yerror !== nothing
        push!(
            plots,
            errorbars!(
                axis,
                x,
                y,
                yerror;
                color,
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
                color,
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
    error_color = get(attributes, :color, :black)
    return _line_errors!(plots, axis, x, y, xerror, yerror, visible, error_color)
end
