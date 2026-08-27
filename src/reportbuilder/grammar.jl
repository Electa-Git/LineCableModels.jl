"""
$(TYPEDEF)

Supertype for policies consumed by [`report`](@ref).
"""
abstract type AbstractReportDefinition end

"""
$(TYPEDEF)

Hold the table, optional PlotBuilder artifact, and written output produced by
[`report`](@ref).

$(TYPEDFIELDS)
"""
struct ReportArtifact{T, I, O}
    "Human-facing table or structured collection of tables."
    table::T
    "Optional backend-neutral plot artifact."
    illustration::I
    "Written destination or handle, or `nothing` for an in-memory report."
    output::O
end

"""
$(TYPEDEF)

Request one generic wide table from a completed scientific source.

The request and unit tuples follow [`observables`](@ref). `illustration` may be
a PlotBuilder definition type; `plot_options` are passed to `make_render`.

$(TYPEDFIELDS)
"""
struct TableReportDefinition{R <: Tuple, U <: Tuple, P, O <: NamedTuple} <:
       AbstractReportDefinition
    "Positional scientific observation requests."
    requests::R
    "Optional display-unit overrides aligned with `requests`."
    units::U
    "PlotBuilder definition type, or `nothing` for a table-only report."
    illustration::P
    "Keyword options passed to PlotBuilder."
    plot_options::O
    "Whether detached display residue is replaced with exact zero."
    clip::Bool
end

function TableReportDefinition(
        requests::Tuple;
        units::Tuple = (),
        illustration = nothing,
        plot_options::NamedTuple = (;),
        clip::Bool = true
)
    return TableReportDefinition(requests, units, illustration, plot_options, clip)
end

"Reject unsupported definition/source pairs before report construction."
function entitle end

"Publish the observations required by a report definition."
function select end

"Construct the report-owned table representation."
function tabulate end

"Construct an optional PlotBuilder artifact."
function illustrate end

"Encode a supported report representation."
function encode end

"Write an encoded report when the definition requests an external artifact."
function write end

@required AbstractReportDefinition begin
    entitle(::AbstractReportDefinition, source)
    select(::AbstractReportDefinition, source)
    tabulate(::AbstractReportDefinition, source, published)
    illustrate(::AbstractReportDefinition, source, published, table)
    encode(::AbstractReportDefinition, source, published, table, illustration)
    write(::AbstractReportDefinition, source, published, table, illustration, encoded)
end

"Return the completed report artifact."
function finish end

function entitle(definition::TableReportDefinition, source)
    validate_observables(source, definition.requests, definition.units)
    return source
end

function select(definition::TableReportDefinition, source)
    return observables(
        source,
        definition.requests;
        units = definition.units,
        clip = definition.clip
    )
end

function _table_column(values)
    values isa Number && return [values]
    values isa AbstractVector && return collect(values)
    values isa AbstractArray && return collect(vec(values))
    return [values]
end

function _observation_names(published::Tuple)
    names = Tuple(Symbol(Units.symbol(payload.quantity)) for payload in published)
    all(name -> !isempty(string(name)), names) || throw(ArgumentError(
        "every published quantity must define a nonempty table symbol",
    ))
    length(unique(names)) == length(names) || throw(ArgumentError(
        "published quantities must have distinct table symbols",
    ))
    return names
end

function _observation_contract(names::Tuple, published::Tuple)
    records = map(published) do payload
        (; quantity = payload.quantity, unit = payload.unit)
    end
    return NamedTuple{names}(records)
end

function _observation_columns!(table::DataFrame, names::Tuple, published::Tuple)
    metadata!(
        table,
        "observation_columns",
        _observation_contract(names, published),
        style = :note
    )
    return table
end

"""
$(TYPEDSIGNATURES)

Return the scientific quantity and display unit owned by every observed column
of `table`.

# Returns

- A named tuple keyed by DataFrame column name. Each value contains `quantity`
  and `unit`.

# Errors

- Throws when `table` was not constructed from an owned observation
  publication.
"""
function observation_columns(table::DataFrame)
    return metadata(table, "observation_columns")
end

function tabulate(::TableReportDefinition, source, published::Tuple)
    columns = map(payload -> _table_column(payload.values), values(published))
    isempty(columns) ||
        all(length(column) == length(first(columns)) for column in columns) ||
        throw(DimensionMismatch("published report columns must have equal lengths"))
    names = _observation_names(published)
    table = DataFrame(NamedTuple{names}(columns))
    return _observation_columns!(table, names, published)
end

function illustrate(definition::TableReportDefinition, source, published, table)
    definition.illustration === nothing && return nothing
    return PlotBuilder.make_render(
        definition.illustration,
        published;
        definition.plot_options...
    )
end

encode(::TableReportDefinition, source, published, table, illustration) = nothing
write(::TableReportDefinition, source, published, table, illustration, ::Nothing) = nothing

function finish(
        ::AbstractReportDefinition,
        source,
        published,
        table,
        illustration,
        encoded,
        written
)
    return ReportArtifact(table, illustration, written)
end

"""
$(TYPEDSIGNATURES)

Build a report through `entitle`, `select`, `tabulate`, `illustrate`, `encode`,
`write`, and `finish`, in that order.

# Arguments

- `definition`: Report policy and requested output.
- `source`: Completed scientific result or published-product owner.

# Returns

- A [`ReportArtifact`](@ref), or the value returned by a specialised `finish`
  method.

# Errors

- Throws when the definition does not accept the source or requests an
  unsupported observation.
"""
@orchestrator AbstractReportDefinition function report(
        definition::AbstractReportDefinition,
        source
)
    entitled = entitle(definition, source)
    published = select(definition, entitled)
    table = tabulate(definition, entitled, published)
    illustration = illustrate(definition, entitled, published, table)
    encoded = encode(definition, entitled, published, table, illustration)
    written = write(
        definition,
        entitled,
        published,
        table,
        illustration,
        encoded
    )
    return finish(
        definition,
        entitled,
        published,
        table,
        illustration,
        encoded,
        written
    )
end
