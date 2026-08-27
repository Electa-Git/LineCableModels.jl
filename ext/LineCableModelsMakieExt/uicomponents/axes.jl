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
