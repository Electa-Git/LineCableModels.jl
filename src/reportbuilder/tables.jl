"""
$(TYPEDEF)

Publish a [`DataModel.CableConstants`](@ref) result as one R/L/C table.

$(TYPEDFIELDS)
"""
struct CableConstantsTableDefinition <: AbstractReportDefinition
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end
CableConstantsTableDefinition() = CableConstantsTableDefinition(true)

"""
$(TYPEDEF)

Define the observable requests and display units for one wide line-parameter
table.

$(TYPEDFIELDS)
"""
struct LineParametersTableDefinition{Q <: Tuple, U} <: AbstractReportDefinition
    "Explicit observable requests in output-column order."
    requests::Q
    "SI prefix used to display frequency."
    frequency_unit::Symbol
    "Length prefix used for per-length quantities."
    length_unit::Symbol
    "Optional display-unit overrides resolved from the requests."
    quantity_units::U
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end

"""
$(TYPEDEF)

Publish a [`Engine.LineParametersBenchmark`](@ref) as one wide comparison
table.

$(TYPEDFIELDS)
"""
struct BenchmarkTableDefinition <: AbstractReportDefinition
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end
BenchmarkTableDefinition() = BenchmarkTableDefinition(true)

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
entitle(::BenchmarkTableDefinition, source::Engine.LineParametersBenchmark) = source

function select(definition::CableConstantsTableDefinition, source::DataModel.CableConstants)
    requests = (R, L, C)
    targets = map(requests) do selector
        Units.native_unit(Units.quantity(selector), basis(source))
    end
    return observables(source, requests; units = targets, clip = definition.clip)
end

function tabulate(::CableConstantsTableDefinition, source, published::Tuple)
    names = _observation_names(published)
    columns = map(payload -> [payload.values], published)
    table = DataFrame(NamedTuple{names}(columns))
    metadata!(table, "basis", basis(source), style = :note)
    metadata!(table, "row_order", (), style = :note)
    return _observation_columns!(table, names, published)
end

function _line_definition(
        requests::Tuple,
        frequency_unit::Symbol,
        length_unit::Symbol,
        quantity_units,
        clip::Bool
)
    isempty(requests) && throw(ArgumentError(
        "line tables require at least one explicit observable request",
    ))
    all(request -> request isa Tuple, requests) || throw(ArgumentError(
        "line tables require requests constructed with @observe",
    ))
    return LineParametersTableDefinition(
        requests,
        frequency_unit,
        length_unit,
        quantity_units,
        clip
    )
end

function _resolve_table_observations(source::Engine.LineParameters, requests::Tuple)
    resolved = map(request -> observation_request(source, request), requests)
    all(request -> length(request.indices) == 3, resolved) || throw(ArgumentError(
        "line tables require row, column, and frequency indices",
    ))
    published_quantities = map(request -> request.quantity, resolved)
    length(unique(published_quantities)) == length(published_quantities) ||
        throw(ArgumentError(
            "line tables do not accept duplicate scientific quantities",
        ))

    dimensions = size(Z(source))
    coordinates = map(resolved) do request
        map(observation_indices, request.indices, dimensions)
    end
    all(==(first(coordinates)), coordinates) || throw(DimensionMismatch(
        "all line-table requests must select the same row, column, and frequency indices",
    ))
    return resolved, first(coordinates)
end

function select(definition::LineParametersTableDefinition, source::Engine.LineParameters)
    resolved, coordinates = _resolve_table_observations(source, definition.requests)
    materialized = map(resolved) do request
        materialize_observation(request, coordinates)
    end
    quantity_targets = unit_targets(
        materialized,
        basis(source);
        length_prefix = definition.length_unit,
        overrides = definition.quantity_units
    )
    frequency_target = Units.units(definition.frequency_unit, :hertz)
    requests = ((frequencies, last(coordinates)), materialized...)
    targets = (frequency_target, quantity_targets...)
    published = observables(source, requests; units = targets, clip = definition.clip)
    frequency = first(published)
    observations = Base.tail(published)
    expected = Tuple(length(indices) for indices in coordinates)
    all(payload -> size(payload.values) == expected, observations) || throw(
        DimensionMismatch("published line quantities do not share one coordinate grid"),
    )
    length(frequency.values) == last(expected) || throw(DimensionMismatch(
        "frequency count does not match published line quantities",
    ))
    return (; frequency, observations, coordinates)
end

function _line_columns(selected)
    rows, columns, samples = selected.coordinates
    frequency_values = selected.frequency.values
    frequency_column = [frequency_values[k]
                        for k in eachindex(samples) for _ in rows for _ in columns]
    row_column = [row for _ in samples for row in rows for _ in columns]
    column_column = [column for _ in samples for _ in rows for column in columns]
    values = map(selected.observations) do payload
        [payload.values[local_row, local_column, local_frequency]
         for local_frequency in eachindex(samples)
         for local_row in eachindex(rows)
         for local_column in eachindex(columns)]
    end
    return frequency_column, row_column, column_column, values
end

function _tabulate_line_parameters(source, selected)
    frequency, rows, columns, quantity_values = _line_columns(selected)
    quantity_names = _observation_names(selected.observations)
    table = DataFrame(merge(
        (; frequency, row = rows, column = columns),
        NamedTuple{quantity_names}(quantity_values)
    ))
    metadata!(table, "basis", basis(source), style = :note)
    metadata!(table, "row_order", (:frequency, :row, :column), style = :note)
    observation_names = (:frequency, quantity_names...)
    published = (selected.frequency, selected.observations...)
    return _observation_columns!(table, observation_names, published)
end

function tabulate(::LineParametersTableDefinition, source, selected)
    return _tabulate_line_parameters(source, selected)
end

function select(
        definition::BenchmarkTableDefinition,
        comparison::Engine.LineParametersBenchmark
)
    requests = (
        (Z, Engine.absolute_error),
        (Z, Engine.relative_error),
        (Y, Engine.absolute_error),
        (Y, Engine.relative_error)
    )
    return observables(comparison, requests; clip = definition.clip)
end

function tabulate(::BenchmarkTableDefinition, source, published::Tuple)
    dimensions = size(first(published).values)
    all(payload -> size(payload.values) == dimensions, published) || throw(
        DimensionMismatch("benchmark quantities do not share one matrix shape"),
    )
    rows = [row for row in 1:dimensions[1] for _ in 1:dimensions[2]]
    columns = [column for _ in 1:dimensions[1] for column in 1:dimensions[2]]
    names = _observation_names(published)
    values = map(payload -> collect(vec(transpose(payload.values))), published)
    table = DataFrame(merge(
        (; row = rows, column = columns),
        NamedTuple{names}(values)
    ))
    metadata!(table, "basis", basis(source), style = :note)
    metadata!(table, "row_order", (:row, :column), style = :note)
    return _observation_columns!(table, names, published)
end

"""
$(TYPEDSIGNATURES)

Return the resistance, inductance, and capacitance of `constants` as one wide
table in native per-length units.
"""
function DataFrame(
        constants::DataModel.CableConstants;
        clip::Bool = true
)::DataFrame
    return report(CableConstantsTableDefinition(clip), constants).table::DataFrame
end

"""
$(TYPEDSIGNATURES)

Return one wide line-parameter table. The first three columns are `frequency`,
`row`, and `column`; each explicit observable request adds one physical-quantity
column.
"""
function DataFrame(
        parameters::Engine.LineParameters,
        requests::Tuple;
        freq_unit::Symbol = :base,
        length_unit::Symbol = :kilo,
        quantity_units = nothing,
        clip::Bool = true
)
    definition = _line_definition(
        requests,
        freq_unit,
        length_unit,
        quantity_units,
        clip
    )
    return report(definition, parameters).table::DataFrame
end

"""
$(TYPEDSIGNATURES)

Return one wide table containing the absolute and relative RMS errors owned by
`comparison`.
"""
function DataFrame(
        comparison::Engine.LineParametersBenchmark;
        clip::Bool = true
)
    return report(BenchmarkTableDefinition(clip), comparison).table::DataFrame
end
