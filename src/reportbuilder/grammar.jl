"""
$(TYPEDEF)

Supertype for policies consumed by [`report`](@ref).
"""
abstract type AbstractReportDefinition end

"""
$(TYPEDEF)

Hold the table, optional plot, and written output produced by
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

The request and unit tuples follow [`observables`](@ref). Set `illustration =
true` to plot the resulting [`ObservationPublication`](@ref), or provide a
callable that accepts it. `plot_options` are forwarded as keywords.

$(TYPEDFIELDS)
"""
struct TableReportDefinition{R <: Tuple, U <: Tuple, P, O <: NamedTuple} <:
       AbstractReportDefinition
    "Positional scientific observation requests."
    requests::R
    "Optional display-unit overrides aligned with `requests`."
    units::U
    "`true`, a publication-plotting callable, or `nothing` for table only."
    illustration::P
    "Keyword options passed to the illustration call."
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

"Publish the observations required by a report definition."
function select end

"Construct the report-owned table representation."
function tabulate end

"Construct an optional plot from the published observations."
function illustrate end

"Encode a supported report representation."
function encode end

"Write an encoded report when the definition requests an external artifact."
function write end

@required AbstractReportDefinition begin
    select(::AbstractReportDefinition, source)
    tabulate(::AbstractReportDefinition, source, published)
end

illustrate(::AbstractReportDefinition, source, published, table) = nothing
encode(::AbstractReportDefinition, source, published, table, illustration) = nothing
function write(
        ::AbstractReportDefinition, source, published, table, illustration, ::Nothing)
    nothing
end

function select(definition::TableReportDefinition, source)
    return observables(
        source,
        definition.requests;
        units = definition.units,
        clip = definition.clip
    )
end

function _observation_columns!(table::DataFrame, published::ObservationPublication)
    metadata!(
        table,
        "observation_columns",
        published.metadata.observation_columns,
        style = :note
    )
    metadata!(table, "basis", published.metadata.basis, style = :note)
    metadata!(table, "row_order", published.metadata.row_order, style = :note)
    return table
end

"""
$(TYPEDSIGNATURES)

Construct one DataFrame from an explicit detached observation publication.
"""
function DataFrame(published::ObservationPublication)
    return _observation_columns!(DataFrame(published.columns), published)
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

function tabulate(::TableReportDefinition, source, published::ObservationPublication)
    return DataFrame(published)
end

function illustrate(definition::TableReportDefinition, source, published, table)
    illustration = definition.illustration
    (illustration === nothing || illustration === false) && return nothing
    illustration === true && return PlotBuilder.plot(published; definition.plot_options...)
    applicable(illustration, published) || throw(ArgumentError(
        "illustration must be true, a callable accepting the publication, or nothing",
    ))
    return illustration(published; definition.plot_options...)
end

"""
$(TYPEDSIGNATURES)

Build a report through `select`, `tabulate`, optional `illustrate`, optional
`encode`, and optional `write`, in that order.

# Arguments

- `definition`: Report policy and requested output.
- `source`: Completed scientific result or published-product owner.

# Returns

- A [`ReportArtifact`](@ref).

# Errors

- Throws when the definition does not accept the source or requests an
  unsupported observation.
"""
function report(
        definition::AbstractReportDefinition,
        source
)
    published = select(definition, source)
    table = tabulate(definition, source, published)
    illustration = illustrate(definition, source, published, table)
    encoded = encode(definition, source, published, table, illustration)
    written = write(
        definition,
        source,
        published,
        table,
        illustration,
        encoded
    )
    return ReportArtifact(table, illustration, written)
end
