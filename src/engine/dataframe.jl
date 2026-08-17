import DataFrames: DataFrame, metadata!

const _LP_FREQ_COL = :frequency

function _frequency_vector(object, provided)
    if provided === nothing
        return Float64.(collect(axes(object, 3)))
    end
    values = Float64.(collect(provided))
    length(values) == size(object, 3) || throw(
        DimensionMismatch("frequency vector length does not match the parameter depth"),
    )
    all(isfinite, values) || throw(ArgumentError("frequencies must be finite"))
    return values
end

const _DATAFRAME_COLUMN = Dict(
    :R => :R,
    :L => :L,
    :G => :G,
    :C => :C,
    :Z_re => :real,
    :Z_im => :imag,
    :Z_abs => :magnitude,
    :Z_angle => :angle,
    :Y_re => :real,
    :Y_im => :imag,
    :Y_abs => :magnitude,
    :Y_angle => :angle
)

function _clip_field(value::Real, tolerance)
    isfinite(value) || return value
    return abs(value) <= tolerance ? zero(value) : value
end

_clip_field(value, _) = value

function _dataframe_unit_label(component, object_basis, length_unit, quantity_units)
    quantity, target,
    factor = _component_unit(
        component,
        object_basis,
        length_unit,
        quantity_units
    )
    return quantity, target, factor, UnitHandler.get_label(target)
end

function _matrix_dataframes(
        object,
        frequency_values,
        component_names;
        frequency_unit,
        length_unit,
        quantity_units,
        tolerance
)
    isempty(component_names) && throw(
        ArgumentError("at least one quantity is required for each matrix family"),
    )
    frequency_quantity = UnitHandler.QuantityTag{:freq}()
    frequency_target = UnitHandler.units(frequency_unit, :hertz)
    frequency_factor = UnitHandler.scale_factor(
        UnitHandler.default_unit(frequency_quantity),
        frequency_target
    )
    displayed_frequency = frequency_values .* frequency_factor
    component_arrays = Dict(
        component => _line_component_values(Val(component), object, frequency_values)
    for
    component in component_names
    )
    row_count, column_count, _ = size(object)
    frames = Matrix{DataFrame}(undef, row_count, column_count)

    for row in 1:row_count, column in 1:column_count

        frame = DataFrame(_LP_FREQ_COL => displayed_frequency)
        unit_map = Dict{Symbol, String}(
            _LP_FREQ_COL => UnitHandler.get_label(frequency_target),
        )
        for component in component_names
            _, target,
            factor,
            unit_label = _dataframe_unit_label(
                component,
                basis(object),
                length_unit,
                quantity_units
            )
            values = collect(view(component_arrays[component], row, column, :)) .* factor
            column_name = _DATAFRAME_COLUMN[component]
            frame[!, column_name] = _clip_field.(values, tolerance)
            unit_map[column_name] = unit_label
        end
        metadata!(frame, "units", unit_map, style = :note)
        frames[row, column] = frame
    end
    return frames
end

"""
$(TYPEDSIGNATURES)

Convert each matrix entry to a frequency-indexed `DataFrame`. The result is an
`n × n` matrix of frames. The optional accessor tuple selects the displayed
quantities. By default, series impedance is represented by its real and
imaginary parts, and shunt admittance is represented likewise.

Container [`basis`](@ref) determines whether units are per length or total.
"""
function DataFrame(
        parameters::Union{SeriesImpedance, ShuntAdmittance},
        quantities::Tuple = ();
        freqs = nothing,
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        tol::Real = sqrt(eps(Float64))
)
    isfinite(tol) && tol >= 0 || throw(
        ArgumentError("tol must be finite and nonnegative"),
    )
    frequency_values = _frequency_vector(parameters, freqs)
    components = _resolve_line_components(parameters, quantities)
    return _matrix_dataframes(
        parameters,
        frequency_values,
        components;
        frequency_unit = freq_unit,
        length_unit,
        quantity_units,
        tolerance = float(tol)
    )
end

"""
$(TYPEDSIGNATURES)

Return `(series, shunt)`, two matrices of frequency-indexed `DataFrame`s,
using the frequencies, basis, and domain stored by `parameters`. The optional
accessor tuple must select at least one series and one shunt quantity.
"""
function DataFrame(
        parameters::LineParameters,
        quantities::Tuple = ();
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        tol::Real = sqrt(eps(Float64))
)
    isfinite(tol) && tol >= 0 || throw(
        ArgumentError("tol must be finite and nonnegative"),
    )
    components = _resolve_line_components(parameters, quantities)
    series_components = Tuple(
        component for component in components if component in _SERIES_COMPONENTS
    )
    shunt_components = Tuple(
        component for component in components if component in _SHUNT_COMPONENTS
    )
    isempty(series_components) && throw(ArgumentError(
        "LineParameters presentation requires a series accessor; use DataFrame(Y(parameters), ...) for a shunt-only table",
    ))
    isempty(shunt_components) && throw(ArgumentError(
        "LineParameters presentation requires a shunt accessor; use DataFrame(Z(parameters), ...) for a series-only table",
    ))
    frequency_values = frequencies(parameters)
    common = (;
        frequency_unit = freq_unit,
        length_unit,
        quantity_units,
        tolerance = float(tol)
    )
    series = _matrix_dataframes(
        Z(parameters),
        frequency_values,
        series_components;
        common...
    )
    shunt = _matrix_dataframes(
        Y(parameters),
        frequency_values,
        shunt_components;
        common...
    )
    return series, shunt
end
