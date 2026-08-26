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
    "Scientific selectors or explicit observable requests to publish."
    requests::Q
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
        requests::Tuple,
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
        requests,
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

_default_table_requests(::Engine.LineParameters) = (R, X, G, B)
_default_table_requests(::Engine.SeriesImpedance) = (R, X)
_default_table_requests(::Engine.ShuntAdmittance) = (G, B)

function _table_request_identity(request::Tuple)
    length(request) == 5 && return (request[1], request[2])
    return first(request)
end

function _table_request_quantity(request)
    identity = _table_request_identity(request)
    return identity isa Tuple ? Units.quantity(identity...) : Units.quantity(identity)
end

function _table_request_indices(request::Tuple)
    offset = _table_request_identity(request) isa Tuple ? 2 : 1
    length(request) == offset + 3 || throw(ArgumentError(
        "line tables require row, column, and frequency indices",
    ))
    return request[(offset + 1):end]
end

function _table_requests(source, requested::Tuple)
    selected = isempty(requested) ? _default_table_requests(source) : requested
    requests = map(selected) do request
        request isa Function && return (request, Colon(), Colon(), Colon())
        request isa Tuple || throw(ArgumentError(
            "line tables accept selectors or explicit observable request tuples",
        ))
        _table_request_indices(request)
        return request
    end
    identities = map(_table_request_identity, requests)
    length(unique(identities)) == length(identities) || throw(ArgumentError(
        "line tables do not accept duplicate scientific quantities",
    ))
    return requests
end

function _table_indices(selector)
    selector isa Colon && return selector
    selector isa Integer && return [Int(selector)]
    selector isa AbstractRange && return collect(Int, selector)
    selector isa AbstractVector && return collect(Int, selector)
    throw(ArgumentError("observable indices must be integers, ranges, vectors, or `:`"))
end

function _table_coordinates(request)
    row, column, sample = _table_request_indices(request)
    rows = _table_indices(row)
    columns = _table_indices(column)
    samples = _table_indices(sample)
    return rows, columns, samples
end

function _materialized_table_request(source, frequencies, request, coordinates)
    rows, columns, samples = coordinates
    identity = _table_request_identity(request)
    prefix = identity isa Tuple ? identity : (identity,)
    source isa Engine.SeriesImpedance && identity === L &&
        (prefix = (L, frequencies))
    source isa Engine.ShuntAdmittance && identity === C &&
        (prefix = (C, frequencies))
    return (prefix..., rows, columns, samples)
end

function _resolved_coordinates(coordinates, payload)
    dimensions = size(payload.values)
    return map(coordinates, dimensions) do coordinate, count
        coordinate isa Colon ? collect(1:count) : coordinate
    end
end

function select(definition::LineParametersTableDefinition, source)
    frequency_values = _frequency_values(source, definition.freqs)
    requests = _table_requests(source, definition.requests)
    coordinates = map(_table_coordinates, requests)
    sample_indices = last.(coordinates)
    all(==(first(sample_indices)), sample_indices) || throw(DimensionMismatch(
        "all line-table requests must select the same frequency indices",
    ))
    names = Tuple(Symbol("request_$index") for index in eachindex(requests))
    materialized = map(requests, coordinates) do request, request_coordinates
        _materialized_table_request(
            source,
            frequency_values,
            request,
            request_coordinates
        )
    end
    named_requests = NamedTuple{names}(materialized)
    targets = unit_targets(
        named_requests,
        basis(source);
        length_prefix = definition.length_unit,
        overrides = definition.quantity_units
    )
    if source isa Engine.LineParameters
        frequency_target = Units.units(definition.frequency_unit, :hertz)
        combined_requests = merge(
            (frequency = (frequencies, first(sample_indices)),),
            named_requests
        )
        combined_units = merge((frequency = frequency_target,), targets)
        combined = observables(source, combined_requests; units = combined_units)
        frequency = combined.frequency
        published = NamedTuple{names}(Tuple(getproperty(combined, name) for name in names))
    else
        frequency = _standalone_frequency_payload(
            frequency_values[first(sample_indices)],
            definition.frequency_unit
        )
        published = observables(source, named_requests; units = targets)
    end
    resolved_coordinates = map(_resolved_coordinates, coordinates, values(published))
    families = map(payload -> Units.family(payload.quantity), values(published))
    return (;
        frequency,
        published,
        requests,
        coordinates = resolved_coordinates,
        families
    )
end

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

_family_name(::Val{:series}) = :series
_family_name(::Val{:shunt}) = :shunt

function tabulate(definition::LineParametersTableDefinition, source, selected)
    frequency_values = selected.frequency.values
    Tfrequency = eltype(frequency_values)
    Tvalue = promote_type((eltype(payload.values) for payload in values(selected.published))...)
    family_column = Symbol[]
    row_column = Int[]
    column_column = Int[]
    frequency_column = Tfrequency[]
    quantity_column = Symbol[]
    value_column = Tvalue[]
    unit_column = String[]
    family_rank = Int[]
    frequency_rank = Int[]
    request_rank = Int[]
    units_metadata = Dict{Tuple{Symbol, Symbol}, String}()
    headings_metadata = Dict{Tuple{Symbol, Symbol}, String}()
    requests_metadata = Dict{Tuple{Symbol, Symbol}, Any}()

    for (request_index, payload) in enumerate(values(selected.published))
        family = _family_name(selected.families[request_index])
        rows, columns, _ = selected.coordinates[request_index]
        row_count, column_count, frequency_count = size(payload.values)
        frequency_count == length(frequency_values) || throw(DimensionMismatch(
            "frequency count does not match line-parameter samples",
        ))
        quantity = Symbol(Units.symbol(payload.quantity))
        metadata_key = (family, quantity)
        units_metadata[metadata_key] = Units.label(payload.unit)
        headings_metadata[metadata_key] = Units.label(payload.quantity, payload.unit)
        requests_metadata[metadata_key] = selected.requests[request_index]
        for local_row in 1:row_count, local_column in 1:column_count

            for frequency_index in eachindex(frequency_values)
                push!(family_column, family)
                push!(row_column, rows[local_row])
                push!(column_column, columns[local_column])
                push!(frequency_column, frequency_values[frequency_index])
                push!(quantity_column, quantity)
                push!(family_rank, family === :series ? 1 : 2)
                push!(frequency_rank, frequency_index)
                push!(request_rank, request_index)
                push!(
                    value_column,
                    clip(
                        payload.values[local_row, local_column, frequency_index],
                        definition.tolerance
                    )
                )
                push!(unit_column, Units.label(payload.unit))
            end
        end
    end

    order = sortperm(eachindex(family_column); by = index -> (
        family_rank[index],
        row_column[index],
        column_column[index],
        frequency_rank[index],
        request_rank[index]
    ))
    table = DataFrame(
        family = family_column[order],
        row = row_column[order],
        column = column_column[order],
        frequency = frequency_column[order],
        quantity = quantity_column[order],
        value = value_column[order],
        unit = unit_column[order]
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
