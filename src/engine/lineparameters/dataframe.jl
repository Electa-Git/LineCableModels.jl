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
    return quantity, target, factor, Units.label(target)
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
    frequency_quantity = Units.QuantityTag{:frequency}()
    frequency_target = Units.units(frequency_unit, :hertz)
    frequency_factor = Units.scale_factor(
        Units.native_unit(frequency_quantity),
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
            _LP_FREQ_COL => Units.label(frequency_target),
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
        parameters.Z,
        frequency_values,
        series_components;
        common...
    )
    shunt = _matrix_dataframes(
        parameters.Y,
        frequency_values,
        shunt_components;
        common...
    )
    return series, shunt
end

function _comparison_floor(zero_atol::NamedTuple, quantity::Symbol)
    names = propertynames(zero_atol)
    length(names) == 2 && all(name -> name in names, (:Z, :Y)) || throw(ArgumentError(
        "zero_atol must contain exactly the Z and Y fields",
    ))
    value = getproperty(zero_atol, quantity)
    value isa Real || throw(ArgumentError(
        "zero_atol.$quantity must be a real number",
    ))
    isfinite(value) && value >= 0 || throw(ArgumentError(
        "zero_atol.$quantity must be finite and nonnegative",
    ))
    return value
end

function _comparison_rows(quantity::Symbol, error::RMSError, floor::Real)
    T = eltype(error.absolute)
    rows = Int[]
    columns = Int[]
    absolute = T[]
    relative = Union{Missing, T}[]
    relative_status = Symbol[]
    for index in CartesianIndices(error.absolute)
        below_floor = floor > 0 && error.absolute[index] <= floor
        push!(rows, index[1])
        push!(columns, index[2])
        push!(absolute, error.absolute[index])
        push!(relative, below_floor ? missing : error.relative[index])
        push!(relative_status, below_floor ? :below_absolute_floor : :reported)
    end
    return DataFrame(
        quantity = fill(quantity, length(rows)),
        row = rows,
        column = columns,
        rms_absolute = absolute,
        rms_relative = relative,
        relative_status = relative_status
    )
end

"""
$(TYPEDSIGNATURES)

Convert the element-wise RMS errors in a [`LineParametersBenchmark`](@ref) to a long-form `DataFrame`.

The `zero_atol` keyword declares separate absolute display floors for Z and Y.
When an absolute RMS error is at or below its floor, `rms_relative` is
`missing` and `relative_status` is `:below_absolute_floor`. The display rule
prevents a large relative value on a physically negligible term from
dominating the table. The stored [`RMSError`](@ref) matrices remain unchanged.

# Arguments

- `comparison`: Element-wise Z/Y comparison.

# Keywords

- `zero_atol`: Named tuple containing nonnegative `Z` and `Y` display floors. `Z` uses the impedance units of the compared result, while `Y` uses its admittance units. A zero floor disables masking for that quantity.

# Returns

- A `DataFrame` with quantity, matrix indices, absolute RMS, displayed relative RMS, and relative-reporting state.

# Errors

- `ArgumentError` when `zero_atol` does not contain exactly `Z` and `Y`, or when either floor is invalid.
"""
function DataFrame(
        comparison::LineParametersBenchmark;
        zero_atol::NamedTuple = (Z = 0.0, Y = 0.0)
)
    impedance_floor = _comparison_floor(zero_atol, :Z)
    admittance_floor = _comparison_floor(zero_atol, :Y)
    frame = vcat(
        _comparison_rows(
            :Z,
            RMSError(
                observe(comparison, Z, absolute_error),
                observe(comparison, Z, relative_error)
            ),
            impedance_floor
        ),
        _comparison_rows(
            :Y,
            RMSError(
                observe(comparison, Y, absolute_error),
                observe(comparison, Y, relative_error)
            ),
            admittance_floor
        )
    )
    metadata!(
        frame,
        "relative RMS",
        "rms_relative is missing when rms_absolute is at or below the corresponding zero_atol display floor; the comparison object retains the raw value",
        style = :note
    )
    metadata!(frame, "zero_atol", zero_atol, style = :note)
    return frame
end
