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

Request one generic table from a completed scientific source.

The request and unit tuples follow [`observables`](@ref). `illustration` may be
a PlotBuilder definition type; `plot_options` are passed to `make_render`.

$(TYPEDFIELDS)
"""
struct TableReportDefinition{R <: NamedTuple, U <: NamedTuple, P, O <: NamedTuple} <:
       AbstractReportDefinition
    "Named scientific observation requests."
    requests::R
    "Optional display-unit overrides keyed like `requests`."
    units::U
    "PlotBuilder definition type, or `nothing` for a table-only report."
    illustration::P
    "Keyword options passed to PlotBuilder."
    plot_options::O
end

function TableReportDefinition(
        requests::NamedTuple;
        units::NamedTuple = (;),
        illustration = nothing,
        plot_options::NamedTuple = (;)
)
    return TableReportDefinition(requests, units, illustration, plot_options)
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

"Return the completed report artifact."
function finish end

function entitle(definition::TableReportDefinition, source)
    validate_observables(source, definition.requests, definition.units)
    return source
end

function select(definition::TableReportDefinition, source)
    return observables(source, definition.requests; units = definition.units)
end

function _table_column(values)
    values isa Number && return [values]
    values isa AbstractVector && return collect(values)
    values isa AbstractArray && return collect(vec(values))
    return [values]
end

function tabulate(::TableReportDefinition, source, published::NamedTuple)
    columns = map(payload -> _table_column(payload.values), values(published))
    isempty(columns) ||
        all(length(column) == length(first(columns)) for column in columns) ||
        throw(DimensionMismatch("published report columns must have equal lengths"))
    table = DataFrame(; NamedTuple{keys(published)}(columns)...)
    unit_labels = Dict(
        key => Units.label(payload.unit) for (key, payload) in pairs(published)
    )
    headings = Dict(
        key => Units.label(payload.quantity, payload.unit)
    for (key, payload) in pairs(published)
    )
    metadata!(table, "units", unit_labels, style = :note)
    metadata!(table, "headings", headings, style = :note)
    return table
end

function illustrate(definition::TableReportDefinition, source, published, table)
    definition.illustration === nothing && return nothing
    return PlotBuilder.make_render(
        definition.illustration,
        source;
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
function report(definition::AbstractReportDefinition, source)
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
