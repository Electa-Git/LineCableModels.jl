struct CableConstantsTable <: AbstractReportDefinition end

struct LineParametersTable{Q <: Tuple, F, U} <: AbstractReportDefinition
    quantities::Q
    freqs::F
    frequency_unit::Symbol
    length_unit::Symbol
    quantity_units::U
    tolerance::Float64
end

struct BenchmarkTable{T <: NamedTuple} <: AbstractReportDefinition
    zero_atol::T
end

const _TableOnlyReport = Union{CableConstantsTable, LineParametersTable, BenchmarkTable}
illustrate(::_TableOnlyReport, source, published, table) = nothing
encode(::_TableOnlyReport, source, published, table, ::Nothing) = nothing
write(::_TableOnlyReport, source, published, table, ::Nothing, ::Nothing) = nothing

entitle(::CableConstantsTable, source::DataModel.CableConstants) = source
entitle(::LineParametersTable, source::Union{
    Engine.LineParameters,
    Engine.SeriesImpedance,
    Engine.ShuntAdmittance
}) = source
entitle(::BenchmarkTable, source::Engine.LineParametersBenchmark) = source

function select(::CableConstantsTable, source::DataModel.CableConstants)
    requests = (R = R, L = L, C = C)
    targets = map(requests) do selector
        Units.native_unit(Units.quantity(selector), basis(source))
    end
    return observables(source, requests; units = targets)
end

function tabulate(::CableConstantsTable, source, published::NamedTuple)
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
    return LineParametersTable(
        quantities,
        freqs,
        frequency_unit,
        length_unit,
        quantity_units,
        float(tolerance)
    )
end

function _frequency_values(source, provided)
    if source isa Engine.LineParameters
        provided === nothing || throw(ArgumentError(
            "LineParameters already owns its frequency samples",
        ))
        return nothing
    end
    if provided === nothing
        throw(ArgumentError(
            "standalone impedance and admittance reports require explicit frequencies",
        ))
    end
    values = Float64.(collect(provided))
    length(values) == size(source, 3) || throw(
        DimensionMismatch("frequency vector length does not match the parameter depth"),
    )
    all(isfinite, values) || throw(ArgumentError("frequencies must be finite"))
    return values
end

function _frequency_payload(source, values, prefix::Symbol)
    target = Units.units(prefix, :hertz)
    if source isa Engine.LineParameters
        return observables(
            source,
            (frequency = (frequencies, Colon()),);
            units = (frequency = target,)
        ).frequency
    end
    quantity = Units.quantity(frequencies)
    factor = Units.scale_factor(Units.native_unit(quantity), target)
    return (; values = values .* factor, quantity, unit = target)
end

_line_request(::Val{:series}, ::typeof(R), frequencies) = ((:R, R),)
_line_request(::Val{:series}, ::typeof(X), frequencies) = ((:X, X),)
_line_request(::Val{:series}, ::typeof(L), frequencies) =
    ((:L, frequencies === nothing ? L : (L, frequencies)),)
_line_request(::Val{:series}, ::typeof(real), frequencies) = ((:real, R),)
_line_request(::Val{:series}, ::typeof(imag), frequencies) = ((:imag, X),)
_line_request(::Val{:series}, ::typeof(abs), frequencies) = ((:magnitude, (Z, abs)),)
_line_request(::Val{:series}, ::typeof(angle), frequencies) = ((:angle, (Z, angle)),)
_line_request(::Val{:series}, ::typeof(Z), frequencies) =
    ((:real, R), (:imag, X))

_line_request(::Val{:shunt}, ::typeof(G), frequencies) = ((:G, G),)
_line_request(::Val{:shunt}, ::typeof(B), frequencies) = ((:B, B),)
_line_request(::Val{:shunt}, ::typeof(C), frequencies) =
    ((:C, frequencies === nothing ? C : (C, frequencies)),)
_line_request(::Val{:shunt}, ::typeof(real), frequencies) = ((:real, G),)
_line_request(::Val{:shunt}, ::typeof(imag), frequencies) = ((:imag, B),)
_line_request(::Val{:shunt}, ::typeof(abs), frequencies) = ((:magnitude, (Y, abs)),)
_line_request(::Val{:shunt}, ::typeof(angle), frequencies) = ((:angle, (Y, angle)),)
_line_request(::Val{:shunt}, ::typeof(Y), frequencies) =
    ((:real, G), (:imag, B))

_line_request(::Val, selector, frequencies) = ()

function _line_requests(family::Val, selectors::Tuple, frequencies)
    entries = Tuple(entry for selector in selectors
    for entry in _line_request(family, selector, frequencies))
    keys = first.(entries)
    length(unique(keys)) == length(keys) || throw(
        ArgumentError("line-parameter accessors select duplicate table columns"),
    )
    return (; (key => request for (key, request) in entries)...)
end

_default_selectors(::Engine.LineParameters) = (Z, Y)
_default_selectors(::Engine.SeriesImpedance) = (Z,)
_default_selectors(::Engine.ShuntAdmittance) = (Y,)

function _request_quantity(request)
    request isa Function && return Units.quantity(request)
    supported_pair = length(request) >= 2 && request[2] isa Function
    return supported_pair ? Units.quantity(request[1], request[2]) :
           Units.quantity(first(request))
end

function _quantity_prefix(quantity_units, key::Symbol, fallback::Symbol)
    quantity_units === nothing && return fallback
    quantity_units isa Symbol && return quantity_units
    quantity_units isa NamedTuple || quantity_units isa AbstractDict || throw(
        ArgumentError("quantity_units must be a prefix, keyed collection, or nothing"),
    )
    return haskey(quantity_units, key) ? quantity_units[key] : fallback
end

function _request_target(request, key, result_basis, length_unit, quantity_units)
    quantity = _request_quantity(request)
    default = Units.display_unit(quantity, result_basis; length_prefix = length_unit)
    isempty(default.numerator) && return default
    selected = _quantity_prefix(
        quantity_units,
        key,
        first(default.numerator).prefix
    )
    selected isa Units.UnitExpr && return selected
    selected isa Symbol || throw(
        ArgumentError("quantity-unit overrides must be prefixes or UnitExpr values"),
    )
    return Units.display_unit(
        quantity,
        result_basis;
        length_prefix = length_unit,
        prefix = selected
    )
end

function _publish_family(source, requests, definition::LineParametersTable)
    isempty(requests) && return nothing
    target_values = map(keys(requests), values(requests)) do key, request
        _request_target(
            request,
            key,
            basis(source),
            definition.length_unit,
            definition.quantity_units
        )
    end
    targets = NamedTuple{keys(requests)}(target_values)
    return observables(source, requests; units = targets)
end

function select(definition::LineParametersTable, source)
    frequencies = _frequency_values(source, definition.freqs)
    selectors = isempty(definition.quantities) ? _default_selectors(source) :
                definition.quantities
    series_requests = source isa Engine.ShuntAdmittance ? (;) :
                      _line_requests(Val(:series), selectors, frequencies)
    shunt_requests = source isa Engine.SeriesImpedance ? (;) :
                     _line_requests(Val(:shunt), selectors, frequencies)
    isempty(series_requests) && isempty(shunt_requests) && throw(ArgumentError(
        "the selected accessors are not reportable for $(typeof(source))",
    ))
    source isa Engine.LineParameters && isempty(series_requests) && throw(ArgumentError(
        "LineParameters presentation requires a series accessor",
    ))
    source isa Engine.LineParameters && isempty(shunt_requests) && throw(ArgumentError(
        "LineParameters presentation requires a shunt accessor",
    ))
    frequency = _frequency_payload(source, frequencies, definition.frequency_unit)
    series = _publish_family(source, series_requests, definition)
    shunt = _publish_family(source, shunt_requests, definition)
    return (; frequency, series, shunt)
end

function _clip_field(value::Real, tolerance)
    isfinite(value) || return value
    return abs(value) <= tolerance ? zero(value) : value
end

_clip_field(value, _) = value

function _matrix_tables(frequency, published::NamedTuple, tolerance)
    first_payload = first(values(published))
    row_count, column_count, _ = size(first_payload.values)
    frames = Matrix{DataFrame}(undef, row_count, column_count)
    for row in 1:row_count, column in 1:column_count
        table = DataFrame(frequency = frequency.values)
        unit_labels = Dict{Symbol, String}(:frequency => Units.label(frequency.unit))
        headings = Dict{Symbol, String}(
            :frequency => Units.label(frequency.quantity, frequency.unit),
        )
        for (key, payload) in pairs(published)
            selected = collect(view(payload.values, row, column, :))
            table[!, key] = _clip_field.(selected, tolerance)
            unit_labels[key] = Units.label(payload.unit)
            headings[key] = Units.label(payload.quantity, payload.unit)
        end
        metadata!(table, "units", unit_labels, style = :note)
        metadata!(table, "headings", headings, style = :note)
        frames[row, column] = table
    end
    return frames
end

function tabulate(definition::LineParametersTable, source, selected)
    series = selected.series === nothing ? nothing :
             _matrix_tables(selected.frequency, selected.series, definition.tolerance)
    shunt = selected.shunt === nothing ? nothing :
            _matrix_tables(selected.frequency, selected.shunt, definition.tolerance)
    source isa Engine.LineParameters && return series, shunt
    return series === nothing ? shunt : series
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

function select(::BenchmarkTable, comparison::Engine.LineParametersBenchmark)
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

function tabulate(definition::BenchmarkTable, source, published)
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
    metadata!(table, "units", (
        Z = Units.label(published.Z_absolute.unit),
        Y = Units.label(published.Y_absolute.unit),
        relative = Units.label(published.Z_relative.unit)
    ), style = :note)
    return table
end

"""
$(TYPEDSIGNATURES)

Return the published resistance, inductance, and capacitance of `constants` as
a table in native per-length units.
"""
function DataFrame(constants::DataModel.CableConstants)::DataFrame
    return report(CableConstantsTable(), constants).table
end

"""
$(TYPEDSIGNATURES)

Return one frequency-indexed table per matrix term of a series-impedance or
shunt-admittance result. Selected quantities are published before tabulation.
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
    return report(definition, parameters).table
end

"""
$(TYPEDSIGNATURES)

Return separate matrices of series and shunt tables for `parameters`.
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
    return report(definition, parameters).table
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
    return report(BenchmarkTable(zero_atol), comparison).table
end
