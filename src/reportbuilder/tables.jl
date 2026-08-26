"""
$(TYPEDEF)

Publish a [`DataModel.CableConstants`](@ref) result as one R/L/C table.
"""
struct CableConstantsTableDefinition <: AbstractReportDefinition end

"""
$(TYPEDEF)

Define the quantities, frequency display, unit overrides, and clipping
tolerance for one long-form line-parameter table.

$(TYPEDFIELDS)
"""
struct LineParametersTableDefinition{Q <: Tuple, F, U} <: AbstractReportDefinition
    "Scientific accessors to publish. An empty tuple selects the source defaults."
    quantities::Q
    "Frequency samples for a standalone matrix, or `nothing` for `LineParameters`."
    freqs::F
    "SI prefix used to display frequency."
    frequency_unit::Symbol
    "Length prefix used for per-length quantities."
    length_unit::Symbol
    "Optional display-unit overrides aligned with the published requests."
    quantity_units::U
    "Absolute threshold used by [`clip`](@ref)."
    tolerance::Float64
end

"""
$(TYPEDEF)

Publish a [`Engine.LineParametersBenchmark`](@ref) as one comparison table.

$(TYPEDFIELDS)
"""
struct BenchmarkTableDefinition{T <: NamedTuple} <: AbstractReportDefinition
    "Absolute Z and Y floors below which relative RMS errors are omitted."
    zero_atol::T
end

const _TableOnlyReport = Union{
    CableConstantsTableDefinition,
    LineParametersTableDefinition,
    BenchmarkTableDefinition
}
illustrate(::_TableOnlyReport, source, published, table) = nothing
encode(::_TableOnlyReport, source, published, table, ::Nothing) = nothing
write(::_TableOnlyReport, source, published, table, ::Nothing, ::Nothing) = nothing

entitle(::CableConstantsTableDefinition, source::DataModel.CableConstants) = source
entitle(::LineParametersTableDefinition, source::Engine.LineParameters) = source
entitle(::LineParametersTableDefinition, source::Engine.SeriesImpedance) = source
entitle(::LineParametersTableDefinition, source::Engine.ShuntAdmittance) = source
entitle(::BenchmarkTableDefinition, source::Engine.LineParametersBenchmark) = source

function select(::CableConstantsTableDefinition, source::DataModel.CableConstants)
    requests = (R = R, L = L, C = C)
    targets = map(requests) do selector
        Units.native_unit(Units.quantity(selector), basis(source))
    end
    return observables(source, requests; units = targets)
end

function tabulate(::CableConstantsTableDefinition, source, published::NamedTuple)
    payloads = values(published)
    return DataFrame(
        parameter = [Units.symbol(payload.quantity) for payload in payloads],
        value = [payload.values for payload in payloads],
        unit = [Units.label(payload.unit) for payload in payloads]
    )
end

function _line_definition(
        quantities::Tuple,
        freqs,
        frequency_unit::Symbol,
        length_unit::Symbol,
        quantity_units,
        tolerance::Real
)
    isfinite(tolerance) && tolerance >= 0 || throw(
        ArgumentError("tol must be finite and nonnegative"),
    )
    return LineParametersTableDefinition(
        quantities,
        freqs,
        frequency_unit,
        length_unit,
        quantity_units,
        float(tolerance)
    )
end

function _frequency_values(::Engine.LineParameters, provided)
    provided === nothing || throw(ArgumentError(
        "LineParameters already owns its frequency samples",
    ))
    return nothing
end

function _standalone_frequency_values(source, provided)
    provided === nothing &&
        throw(ArgumentError(
            "standalone impedance and admittance reports require explicit frequencies",
        ))
    values = Float64.(collect(provided))
    length(values) == size(source, 3) || throw(
        DimensionMismatch("frequency vector length does not match the parameter depth"),
    )
    all(isfinite, values) || throw(ArgumentError("frequencies must be finite"))
    return values
end

function _frequency_values(source::Engine.SeriesImpedance, provided)
    return _standalone_frequency_values(source, provided)
end

function _frequency_values(source::Engine.ShuntAdmittance, provided)
    return _standalone_frequency_values(source, provided)
end

function _frequency_payload(
        source::Engine.LineParameters,
        ::Nothing,
        prefix::Symbol
)
    target = Units.units(prefix, :hertz)
    return observables(
        source,
        (frequency = (frequencies, Colon()),);
        units = (frequency = target,)
    ).frequency
end

function _standalone_frequency_payload(values, prefix::Symbol)
    target = Units.units(prefix, :hertz)
    quantity = Units.quantity(frequencies)
    factor = Units.scale_factor(Units.native_unit(quantity), target)
    return (; values = values .* factor, quantity, unit = target)
end

function _frequency_payload(
        ::Engine.SeriesImpedance,
        values,
        prefix::Symbol
)
    return _standalone_frequency_payload(values, prefix)
end

function _frequency_payload(
        ::Engine.ShuntAdmittance,
        values,
        prefix::Symbol
)
    return _standalone_frequency_payload(values, prefix)
end

_default_selectors(::Engine.LineParameters) = (Z, Y)
_default_selectors(::Engine.SeriesImpedance) = (Z,)
_default_selectors(::Engine.ShuntAdmittance) = (Y,)

_line_column_keys(::typeof(Z), count) = (:real, :imag)
_line_column_keys(::typeof(Y), count) = (:real, :imag)
_line_column_keys(::typeof(real), count) = ntuple(_ -> :real, count)
_line_column_keys(::typeof(imag), count) = ntuple(_ -> :imag, count)
_line_column_keys(::typeof(abs), count) = ntuple(_ -> :magnitude, count)
_line_column_keys(::typeof(angle), count) = ntuple(_ -> :angle, count)
function _line_column_keys(selector::Function, count)
    return ntuple(_ -> nameof(selector), count)
end

function _line_family_requests(entries, parent)
    selected = Tuple(entry for entry in entries if line_parent(last(entry)) === parent)
    names = first.(selected)
    length(unique(names)) == length(names) || throw(
        ArgumentError("line-parameter accessors select duplicate table columns"),
    )
    return (; (key => request for (key, request) in selected)...)
end

function _line_table_requests(source, selectors::Tuple, frequencies)
    resolved = line_requests(source, selectors; frequencies)
    selected = isempty(selectors) ? _default_selectors(source) : selectors
    entries = Tuple(
        key => request
    for selector in selected
    for (key, request) in zip(
        _line_column_keys(
            selector,
            length(line_requests(source, (selector,); frequencies))
        ),
        line_requests(source, (selector,); frequencies)
    )
    )
    Tuple(last.(entries)) == resolved || error(
        "line request expansion changed while assigning table columns",
    )
    return (
        series = _line_family_requests(entries, Z),
        shunt = _line_family_requests(entries, Y)
    )
end

function _publish_family(source, requests, definition::LineParametersTableDefinition)
    isempty(requests) && return nothing
    targets = unit_targets(
        requests,
        basis(source);
        length_prefix = definition.length_unit,
        overrides = definition.quantity_units
    )
    return observables(source, requests; units = targets)
end

function select(definition::LineParametersTableDefinition, source)
    frequencies = _frequency_values(source, definition.freqs)
    selectors = isempty(definition.quantities) ? _default_selectors(source) :
                definition.quantities
    requests = _line_table_requests(source, selectors, frequencies)
    series_requests = requests.series
    shunt_requests = requests.shunt
    isempty(series_requests) && isempty(shunt_requests) &&
        throw(ArgumentError(
            "the selected accessors are not reportable for $(typeof(source))",
        ))
    _validate_line_families(source, series_requests, shunt_requests)
    frequency = _frequency_payload(source, frequencies, definition.frequency_unit)
    series = _publish_family(source, series_requests, definition)
    shunt = _publish_family(source, shunt_requests, definition)
    return (; frequency, series, shunt, requests)
end

function _validate_line_families(
        ::Engine.LineParameters,
        series_requests,
        shunt_requests
)
    isempty(series_requests) && throw(ArgumentError(
        "LineParameters presentation requires a series accessor",
    ))
    isempty(shunt_requests) && throw(ArgumentError(
        "LineParameters presentation requires a shunt accessor",
    ))
    return nothing
end

_validate_line_families(::Engine.SeriesImpedance, series_requests, shunt_requests) = nothing
_validate_line_families(::Engine.ShuntAdmittance, series_requests, shunt_requests) = nothing

"""
$(TYPEDSIGNATURES)

Replace finite real report values within `tolerance` of zero with exact zero.
"""
function clip(value::Real, tolerance)
    isfinite(value) || return value
    return abs(value) <= tolerance ? zero(value) : value
end

clip(value::Complex, _) = value
clip(value::Missing, _) = value

function _line_families(selected)
    return ((:series, selected.series), (:shunt, selected.shunt))
end

function _first_line_payload(selected)
    for (_, published) in _line_families(selected)
        published === nothing || return first(values(published))
    end
    error("line publication produced no scientific values")
end

function tabulate(definition::LineParametersTableDefinition, source, selected)
    frequency_values = selected.frequency.values
    Tfrequency = eltype(frequency_values)
    Tvalue = eltype(_first_line_payload(selected).values)
    family_column = Symbol[]
    row_column = Int[]
    column_column = Int[]
    frequency_column = Tfrequency[]
    quantity_column = Symbol[]
    value_column = Tvalue[]
    unit_column = String[]
    units_metadata = Dict{Tuple{Symbol, Symbol}, String}()
    headings_metadata = Dict{Tuple{Symbol, Symbol}, String}()
    requests_metadata = Dict{Tuple{Symbol, Symbol}, Any}()

    for (family, published) in _line_families(selected)
        published === nothing && continue
        first_payload = first(values(published))
        row_count, column_count, frequency_count = size(first_payload.values)
        frequency_count == length(frequency_values) || throw(DimensionMismatch(
            "frequency count does not match line-parameter samples",
        ))
        for (key, payload) in pairs(published)
            units_metadata[(family, key)] = Units.label(payload.unit)
            headings_metadata[(family, key)] = Units.label(payload.quantity, payload.unit)
            requests_metadata[(family, key)] = getproperty(
                getproperty(selected.requests, family),
                key
            )
        end
        for row in 1:row_count, column in 1:column_count

            for frequency_index in eachindex(frequency_values)
                for (key, payload) in pairs(published)
                    push!(family_column, family)
                    push!(row_column, row)
                    push!(column_column, column)
                    push!(frequency_column, frequency_values[frequency_index])
                    push!(quantity_column, key)
                    push!(
                        value_column,
                        clip(
                            payload.values[row, column, frequency_index],
                            definition.tolerance
                        )
                    )
                    push!(unit_column, Units.label(payload.unit))
                end
            end
        end
    end

    table = DataFrame(
        family = family_column,
        row = row_column,
        column = column_column,
        frequency = frequency_column,
        quantity = quantity_column,
        value = value_column,
        unit = unit_column
    )
    metadata!(table, "basis", basis(source), style = :note)
    metadata!(
        table,
        "frequency_unit",
        Units.label(selected.frequency.unit),
        style = :note
    )
    metadata!(table, "units", units_metadata, style = :note)
    metadata!(table, "headings", headings_metadata, style = :note)
    metadata!(table, "requests", requests_metadata, style = :note)
    metadata!(
        table,
        "row_order",
        (:family, :row, :column, :frequency, :quantity),
        style = :note
    )
    return table
end

function _comparison_floor(zero_atol::NamedTuple, quantity::Symbol)
    names = propertynames(zero_atol)
    length(names) == 2 && all(name -> name in names, (:Z, :Y)) || throw(
        ArgumentError("zero_atol must contain exactly the Z and Y fields"),
    )
    value = getproperty(zero_atol, quantity)
    value isa Real || throw(ArgumentError("zero_atol.$quantity must be a real number"))
    isfinite(value) && value >= 0 || throw(ArgumentError(
        "zero_atol.$quantity must be finite and nonnegative",
    ))
    return value
end

function select(::BenchmarkTableDefinition, comparison::Engine.LineParametersBenchmark)
    requests = (
        Z_absolute = (Z, Engine.absolute_error),
        Z_relative = (Z, Engine.relative_error),
        Y_absolute = (Y, Engine.absolute_error),
        Y_relative = (Y, Engine.relative_error)
    )
    targets = map(requests) do request
        Units.native_unit(Units.quantity(request...), basis(comparison))
    end
    return observables(comparison, requests; units = targets)
end

function _comparison_rows(quantity, absolute, relative, floor)
    T = eltype(absolute.values)
    rows = Int[]
    columns = Int[]
    absolute_values = T[]
    relative_values = Union{Missing, T}[]
    relative_status = Symbol[]
    for index in CartesianIndices(absolute.values)
        below_floor = floor > 0 && absolute.values[index] <= floor
        push!(rows, index[1])
        push!(columns, index[2])
        push!(absolute_values, absolute.values[index])
        push!(relative_values, below_floor ? missing : relative.values[index])
        push!(relative_status, below_floor ? :below_absolute_floor : :reported)
    end
    return DataFrame(
        quantity = fill(quantity, length(rows)),
        row = rows,
        column = columns,
        rms_absolute = absolute_values,
        rms_relative = relative_values,
        relative_status = relative_status
    )
end

function tabulate(definition::BenchmarkTableDefinition, source, published)
    impedance_floor = _comparison_floor(definition.zero_atol, :Z)
    admittance_floor = _comparison_floor(definition.zero_atol, :Y)
    table = vcat(
        _comparison_rows(
            :Z,
            published.Z_absolute,
            published.Z_relative,
            impedance_floor
        ),
        _comparison_rows(
            :Y,
            published.Y_absolute,
            published.Y_relative,
            admittance_floor
        )
    )
    metadata!(
        table,
        "relative RMS",
        "rms_relative is missing when rms_absolute is at or below the corresponding zero_atol display floor; the comparison retains the raw value",
        style = :note
    )
    metadata!(table, "zero_atol", definition.zero_atol, style = :note)
    metadata!(table,
        "units",
        (
            Z = Units.label(published.Z_absolute.unit),
            Y = Units.label(published.Y_absolute.unit),
            relative = Units.label(published.Z_relative.unit)
        ),
        style = :note)
    return table
end

"""
$(TYPEDSIGNATURES)

Return the published resistance, inductance, and capacitance of `constants` as
a table in native per-length units.
"""
function DataFrame(constants::DataModel.CableConstants)::DataFrame
    return report(CableConstantsTableDefinition(), constants).table::DataFrame
end

"""
$(TYPEDSIGNATURES)

Return one long-form table for a series-impedance or shunt-admittance result.
Rows are ordered by family, matrix row, matrix column, frequency, and resolved
quantity.
"""
function DataFrame(
        parameters::Union{Engine.SeriesImpedance, Engine.ShuntAdmittance},
        quantities::Tuple = ();
        freqs = nothing,
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        tol::Real = sqrt(eps(Float64))
)
    definition = _line_definition(
        quantities,
        freqs,
        freq_unit,
        length_unit,
        quantity_units,
        tol
    )
    return report(definition, parameters).table::DataFrame
end

"""
$(TYPEDSIGNATURES)

Return one long-form table containing the series and shunt values in
`parameters`. Rows are ordered by family, matrix row, matrix column, frequency,
and resolved quantity.
"""
function DataFrame(
        parameters::Engine.LineParameters,
        quantities::Tuple = ();
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        tol::Real = sqrt(eps(Float64))
)
    definition = _line_definition(
        quantities,
        nothing,
        freq_unit,
        length_unit,
        quantity_units,
        tol
    )
    return report(definition, parameters).table::DataFrame
end

"""
$(TYPEDSIGNATURES)

Return the element-wise absolute and relative RMS errors in `comparison` as a
long-form table.
"""
function DataFrame(
        comparison::Engine.LineParametersBenchmark;
        zero_atol::NamedTuple = (Z = 0.0, Y = 0.0)
)
    return report(BenchmarkTableDefinition(zero_atol), comparison).table::DataFrame
end
