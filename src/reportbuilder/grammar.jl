"""
$(TYPEDEF)

Supertype for policies consumed by [`report`](@ref).
"""
abstract type AbstractReportDefinition end

"""
$(TYPEDEF)

Hold the table and optional PlotBuilder artifact produced by [`report`](@ref).

$(TYPEDFIELDS)
"""
struct ReportArtifact{T, I}
    "Human-facing table or structured collection of tables."
    table::T
    "Optional backend-neutral plot artifact."
    illustration::I
end

"""
$(TYPEDEF)

Request one generic table from a completed scientific source.

The request and unit tuples follow [`observables`](@ref). `illustration` may be
a PlotBuilder definition type; `plot_options` are passed to `make_render`.

$(TYPEDFIELDS)
"""
struct TableReport{R <: NamedTuple, U <: NamedTuple, P, O <: NamedTuple} <:
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

function TableReport(
        requests::NamedTuple;
        units::NamedTuple = (;),
        illustration = nothing,
        plot_options::NamedTuple = (;)
)
    return TableReport(requests, units, illustration, plot_options)
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

function _request_identity(request, supported)
    request isa Function && return request
    request isa Tuple && !isempty(request) || throw(
        ArgumentError("report requests must be selector functions or nonempty tuples"),
    )
    pair = length(request) >= 2 ? (request[1], request[2]) : nothing
    return pair !== nothing && pair in supported ? pair : first(request)
end

function entitle(definition::TableReport, source)
    applicable(observables, typeof(source)) || throw(ArgumentError(
        "$(typeof(source)) does not declare reportable observations",
    ))
    supported = observables(typeof(source))
    for request in values(definition.requests)
        identity = _request_identity(request, supported)
        identity in supported || throw(ArgumentError(
            "$(typeof(source)) does not publish selector $(repr(identity))",
        ))
    end
    return source
end

function select(definition::TableReport, source)
    return observables(source, definition.requests; units = definition.units)
end

function _table_column(values)
    values isa Number && return [values]
    values isa AbstractVector && return collect(values)
    values isa AbstractArray && return collect(vec(values))
    return [values]
end

function tabulate(::TableReport, source, published::NamedTuple)
    columns = map(payload -> _table_column(payload.values), values(published))
    isempty(columns) || all(length(column) == length(first(columns)) for column in columns) ||
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

illustrate(::AbstractReportDefinition, source, published, table) = nothing
function illustrate(definition::TableReport, source, published, table)
    definition.illustration === nothing && return nothing
    return PlotBuilder.make_render(
        definition.illustration,
        source;
        definition.plot_options...
    )
end

encode(::AbstractReportDefinition, source, published, table, illustration) = nothing
write(::AbstractReportDefinition, source, published, table, illustration, encoded) = nothing

function finish(
        ::AbstractReportDefinition,
        source,
        published,
        table,
        illustration,
        encoded,
        written
)
    return ReportArtifact(table, illustration)
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
