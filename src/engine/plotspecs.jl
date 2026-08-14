struct LineParameterPlotSpec <: PlotBuilder.AbstractPlotSpec end

function _component_unit(
        component::Symbol,
        parameter_basis::Symbol,
        length_unit::Symbol,
        quantity_units
)
    resolved = UnitHandler.line_component_unit(
        component,
        parameter_basis;
        length_unit,
        quantity_units
    )
    return resolved.quantity, resolved.units, resolved.scale
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
    selector === nothing ||
        (selector isa Tuple && length(selector) == 2) ||
        throw(
            ArgumentError("conductor selection must be a tuple (rows, columns) or nothing"),
        )
    row_selector, column_selector = selector === nothing ? (nothing, nothing) : selector
    rows = _indices(row_selector, row_count)
    columns = _indices(column_selector, column_count)
    all(index -> index in 1:row_count, rows) || throw(BoundsError(object, rows))
    all(index -> index in 1:column_count, columns) || throw(BoundsError(object, columns))
    return [(i, j) for i in rows for j in columns]
end

_line_parameter_kind(::SeriesImpedance) = :series
_line_parameter_kind(::ShuntAdmittance) = :shunt

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

function _axis_label(quantity, unit)
    label = UnitHandler.get_label(quantity)
    unit_label = UnitHandler.get_label(unit)
    return isempty(unit_label) ? label : "$label [$unit_label]"
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
        yscale::Symbol,
        export_theme::Symbol,
        open_export::Bool
)
    PlotBuilder._validate_export_theme(export_theme)
    all(isfinite, frequency_values) || throw(ArgumentError("frequencies must be finite"))
    if length(frequency_values) <= 1
        @warn "Frequency vector has $(length(frequency_values)) sample(s); nothing to plot."
        return PlotBuilder.PageSpec[]
    end
    size(object, 3) == length(frequency_values) || throw(
        DimensionMismatch("frequency count does not match line-parameter samples"),
    )
    any(iszero, frequency_values) && mode === :RLCG &&
        throw(
            DomainError(frequency_values, "RLCG plotting is undefined at zero frequency"),
        )
    xscale in (:linear, :log10) || throw(ArgumentError("xscale must be :linear or :log10"))
    yscale in (:linear, :log10) || throw(ArgumentError("yscale must be :linear or :log10"))
    xscale === :log10 && any(<=(0), frequency_values) &&
        throw(
            DomainError(frequency_values, "logarithmic frequency axes require positive frequencies"),
        )

    frequency_quantity = UnitHandler.QuantityTag{:freq}()
    frequency_target = UnitHandler.units(freq_unit, :hertz)
    scaled_frequency = frequency_values .* UnitHandler.scale_factor(
        UnitHandler.default_unit(frequency_quantity),
        frequency_target
    )
    x_exponent = _finite_exponent((scaled_frequency,))
    pairs = _conductor_pairs(object, con)
    pages = PlotBuilder.PageSpec[]

    components = UnitHandler.line_components(_line_parameter_kind(object), mode, coord)
    for component in components
        quantity, target_unit, conversion = _component_unit(
            component,
            basis(object),
            length_unit,
            quantity_units
        )
        values = UnitHandler.line_component_values(
            component,
            object.values,
            frequency_values
        )
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
                    scaled_frequency,
                    curve,
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
            _axis_label(frequency_quantity, frequency_target),
            xscale
        )
        yaxis = PlotBuilder.AxisSpec(
            :y,
            quantity,
            target_unit,
            _axis_label(quantity, target_unit),
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
                    export_theme,
                    open_export,
                    controls = PlotBuilder.control_definitions(xlog = true, ylog = true),
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
        yscale::Symbol = :linear,
        export_theme::Symbol = :default,
        open_export::Bool = true
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
        yscale,
        export_theme,
        open_export
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
