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

function _clip_field(value::Measurements.Measurement, tolerance)
    nominal = abs(Measurements.value(value)) <= tolerance ? 0.0 : Measurements.value(value)
    uncertainty_value = abs(Measurements.uncertainty(value)) <= tolerance ? 0.0 :
                        Measurements.uncertainty(value)
    return Measurements.measurement(nominal, uncertainty_value)
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
        frequency_values;
        mode,
        coord,
        frequency_unit,
        length_unit,
        quantity_units,
        tolerance
)
    frequency_quantity = UnitHandler.QuantityTag{:freq}()
    frequency_target = UnitHandler.units(frequency_unit, :hertz)
    frequency_factor = UnitHandler.scale_factor(
        UnitHandler.default_unit(frequency_quantity),
        frequency_target
    )
    displayed_frequency = frequency_values .* frequency_factor
    component_names = UnitHandler.line_components(
        _line_parameter_kind(object),
        mode,
        coord
    )
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
    DataFrame(parameters::Union{SeriesImpedance,ShuntAdmittance}; kwargs...)

Convert each matrix entry to a frequency-indexed `DataFrame`. The result is an
`n × n` matrix of frames. `mode=:RLCG` returns physical components and
`mode=:ZY` returns Cartesian or polar components selected by `coord`.

Container [`basis`](@ref) determines whether units are per length or total.
"""
function DataFrame(
        parameters::Union{SeriesImpedance, ShuntAdmittance};
        freqs = nothing,
        mode::Symbol = :RLCG,
        coord::Symbol = :cart,
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        tol::Real = sqrt(eps(Float64))
)
    isfinite(tol) && tol >= 0 || throw(
        ArgumentError("tol must be finite and nonnegative"),
    )
    frequency_values = _frequency_vector(parameters, freqs)
    return _matrix_dataframes(
        parameters,
        frequency_values;
        mode,
        coord,
        frequency_unit = freq_unit,
        length_unit,
        quantity_units,
        tolerance = float(tol)
    )
end

"""
    DataFrame(parameters::LineParameters; kwargs...)

Return `(series, shunt)`, two matrices of frequency-indexed `DataFrame`s,
using the frequencies, basis, and domain stored by `parameters`.
"""
function DataFrame(
        parameters::LineParameters;
        mode::Symbol = :RLCG,
        coord::Symbol = :cart,
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        tol::Real = sqrt(eps(Float64))
)
    common = (
        freqs = frequencies(parameters),
        mode = mode,
        coord = coord,
        freq_unit = freq_unit,
        length_unit = length_unit,
        quantity_units = quantity_units,
        tol = tol
    )
    return DataFrame(Z(parameters); common...), DataFrame(Y(parameters); common...)
end
