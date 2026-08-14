struct LineParameterPlotSpec <: PlotBuilder.AbstractPlotSpec end

const _PLOT_QUANTITY = Dict(
    :R => (:resistance, :ohm, :base),
    :X => (:reactance, :ohm, :base),
    :L => (:inductance, :henry, :milli),
    :G => (:conductance, :siemens, :base),
    :B => (:susceptance, :siemens, :base),
    :C => (:capacitance, :farad, :micro),
    :Z_re => (:resistance, :ohm, :base),
    :Z_im => (:reactance, :ohm, :base),
    :Z_abs => ((:impedance, :abs), :ohm, :base),
    :Z_angle => ((:impedance, :angle), :degree, :base),
    :Y_re => (:conductance, :siemens, :base),
    :Y_im => (:susceptance, :siemens, :base),
    :Y_abs => ((:admittance, :abs), :siemens, :base),
    :Y_angle => ((:admittance, :angle), :degree, :base)
)

function _quantity_prefix(quantity_units, component::Symbol, semantic, fallback::Symbol)
    quantity_units === nothing && return fallback
    quantity_units isa Symbol && return quantity_units
    if quantity_units isa NamedTuple || quantity_units isa AbstractDict
        haskey(quantity_units, component) && return quantity_units[component]
        semantic isa Symbol && haskey(quantity_units, semantic) &&
            return quantity_units[semantic]
        return fallback
    end
    throw(ArgumentError("quantity_units must be a Symbol, NamedTuple, dictionary, or nothing"))
end

function _component_unit(
        component::Symbol,
        parameter_basis::Symbol,
        length_unit::Symbol,
        quantity_units
)
    semantic, unit_name, fallback_prefix = _PLOT_QUANTITY[component]
    tag = UnitHandler.QuantityTag{semantic}()
    prefix = _quantity_prefix(quantity_units, component, semantic, fallback_prefix)
    native = UnitHandler.default_unit(tag, parameter_basis)
    target = if unit_name === :degree
        UnitHandler.units(prefix, unit_name)
    elseif parameter_basis === :per_length
        UnitHandler.units(prefix, unit_name; per = (length_unit, :meter))
    else
        UnitHandler.units(prefix, unit_name)
    end
    return tag, target, UnitHandler.scale_factor(native, target)
end

function _indices(selector, count::Int)
    selector === nothing && return collect(1:count)
    selector isa Colon && return collect(1:count)
    selector isa Integer && return [Int(selector)]
    selector isa AbstractRange && return collect(Int, selector)
    selector isa AbstractVector && return collect(Int, selector)
    throw(ArgumentError("conductor selections must be integers, ranges, vectors, `:`, or nothing"))
end

function _conductor_pairs(object, selector)
    row_count, column_count, _ = size(object)
    row_selector, column_selector = selector === nothing ? (nothing, nothing) : selector
    rows = _indices(row_selector, row_count)
    columns = _indices(column_selector, column_count)
    all(index -> index in 1:row_count, rows) || throw(BoundsError(object, rows))
    all(index -> index in 1:column_count, columns) || throw(BoundsError(object, columns))
    return [(i, j) for i in rows for j in columns]
end

function _components(object, mode::Symbol, coord::Symbol)
    mode in (:ZY, :RLCG) || throw(ArgumentError("mode must be :ZY or :RLCG"))
    coord in (:cart, :polar) || throw(ArgumentError("coord must be :cart or :polar"))
    if mode === :RLCG
        object isa SeriesImpedance && return (:R, :L)
        object isa ShuntAdmittance && return (:G, :C)
    elseif object isa SeriesImpedance
        return coord === :cart ? (:Z_re, :Z_im) : (:Z_abs, :Z_angle)
    elseif object isa ShuntAdmittance
        return coord === :cart ? (:Y_re, :Y_im) : (:Y_abs, :Y_angle)
    end
    throw(ArgumentError("unsupported line-parameter plot object $(typeof(object))"))
end

function _component_values(component::Symbol, object, frequency_values)
    values = object.values
    component in (:R, :Z_re) && return real.(values)
    component in (:X, :Z_im) && return imag.(values)
    component === :L && return imag.(values) ./ reshape(2π .* frequency_values, 1, 1, :)
    component in (:G, :Y_re) && return real.(values)
    component in (:B, :Y_im) && return imag.(values)
    component === :C && return imag.(values) ./ reshape(2π .* frequency_values, 1, 1, :)
    component in (:Z_abs, :Y_abs) && return abs.(values)
    component in (:Z_angle, :Y_angle) && return angle.(values) .* (180 / π)
    throw(ArgumentError("unsupported line-parameter component :$component"))
end

function _finite_exponent(curves)
    maximum_value = 0.0
    for curve in curves, sample in curve

        nominal = abs(Measurements.value(sample))
        isfinite(nominal) && (maximum_value = max(maximum_value, nominal))
    end
    iszero(maximum_value) && return 0
    exponent = floor(Int, log10(maximum_value))
    return abs(exponent) < 3 ? 0 : exponent
end

function _axis_label(quantity, unit, exponent::Int)
    label = UnitHandler.get_label(quantity)
    unit_label = UnitHandler.get_label(unit)
    base = isempty(unit_label) ? label : "$label [$unit_label]"
    return exponent == 0 ? base : "$base  × 10^$exponent"
end

function _line_pages(
        object,
        frequency_values;
        mode::Symbol,
        coord::Symbol,
        freq_unit::Symbol,
        length_unit::Symbol,
        quantity_units,
        con,
        fig_size::Tuple{Int, Int},
        xscale::Symbol,
        yscale::Symbol
)
    length(frequency_values) > 1 || return PlotBuilder.PageSpec[]
    size(object, 3) == length(frequency_values) || throw(
        DimensionMismatch("frequency count does not match line-parameter samples"),
    )
    any(iszero, frequency_values) && mode === :RLCG &&
        throw(
            DomainError(frequency_values, "RLCG plotting is undefined at zero frequency"),
        )
    xscale in (:linear, :log10) || throw(ArgumentError("xscale must be :linear or :log10"))
    yscale in (:linear, :log10) || throw(ArgumentError("yscale must be :linear or :log10"))

    frequency_quantity = UnitHandler.QuantityTag{:freq}()
    frequency_target = UnitHandler.units(freq_unit, :hertz)
    scaled_frequency = frequency_values .* UnitHandler.scale_factor(
        UnitHandler.default_unit(frequency_quantity),
        frequency_target
    )
    x_exponent = _finite_exponent((scaled_frequency,))
    displayed_frequency = scaled_frequency ./ 10.0^x_exponent
    pairs = _conductor_pairs(object, con)
    pages = PlotBuilder.PageSpec[]

    for component in _components(object, mode, coord)
        quantity, target_unit, conversion = _component_unit(
            component,
            basis(object),
            length_unit,
            quantity_units
        )
        values = _component_values(component, object, frequency_values)
        curves = [collect(view(values, i, j, :)) .* conversion for (i, j) in pairs]
        active = [index
                  for index in eachindex(curves)
                  if
                  any(value -> abs(Measurements.value(value)) > eps(Float64), curves[index])]
        active_curves = curves[active]
        active_pairs = pairs[active]
        y_exponent = _finite_exponent(active_curves)
        series = PlotBuilder.SeriesSpec[]
        symbol = UnitHandler.get_symbol(quantity)
        for (curve, (i, j)) in zip(active_curves, active_pairs)
            push!(
                series,
                PlotBuilder.SeriesSpec(
                    :line,
                    displayed_frequency,
                    curve ./ 10.0^y_exponent,
                    nothing,
                    "$symbol[$i,$j]";
                    attributes = (; linewidth = 2)
                )
            )
        end
        title = UnitHandler.get_label(quantity)
        xaxis = PlotBuilder.AxisSpec(
            :x,
            frequency_quantity,
            frequency_target,
            _axis_label(frequency_quantity, frequency_target, x_exponent),
            xscale
        )
        yaxis = PlotBuilder.AxisSpec(
            :y,
            quantity,
            target_unit,
            _axis_label(quantity, target_unit, y_exponent),
            yscale
        )
        view_spec = PlotBuilder.ViewSpec(
            xaxis,
            yaxis,
            nothing,
            title,
            series,
            (; component)
        )
        push!(
            pages,
            PlotBuilder.PageSpec(
                title,
                fig_size,
                :single,
                PlotBuilder.ViewSpec[view_spec],
                (;
                    component,
                    x_exponent,
                    y_exponent,
                    controls = PlotBuilder.control_definitions(),
                    configuration = (;
                        mode,
                        coord,
                        freq_unit,
                        length_unit,
                        quantity_units,
                        conductors = con
                    )
                )
            )
        )
    end
    return pages
end

function PlotBuilder.make_render(
        ::Type{LineParameterPlotSpec},
        object::Union{SeriesImpedance, ShuntAdmittance};
        frequencies,
        mode::Symbol = :ZY,
        coord::Symbol = :cart,
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        con = nothing,
        fig_size::Tuple{Int, Int} = (800, 400),
        xscale::Symbol = :linear,
        yscale::Symbol = :linear
)
    pages = _line_pages(
        object,
        collect(frequencies);
        mode,
        coord,
        freq_unit,
        length_unit,
        quantity_units,
        con,
        fig_size,
        xscale,
        yscale
    )
    return PlotBuilder.RenderSpec(LineParameterPlotSpec, pages)
end

function PlotBuilder.make_render(
        ::Type{LineParameterPlotSpec},
        parameters::LineParameters;
        kwargs...
)
    impedance = PlotBuilder.make_render(
        LineParameterPlotSpec,
        parameters.Z;
        frequencies = parameters.f,
        kwargs...
    )
    admittance = PlotBuilder.make_render(
        LineParameterPlotSpec,
        parameters.Y;
        frequencies = parameters.f,
        kwargs...
    )
    return PlotBuilder.RenderSpec(
        LineParameterPlotSpec,
        vcat(impedance.figures, admittance.figures)
    )
end
