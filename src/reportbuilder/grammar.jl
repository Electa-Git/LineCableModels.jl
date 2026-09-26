"""
$(TYPEDEF)

Supertype for definitions consumed by [`report`](@ref).
"""
abstract type AbstractReportDefinition end

"""
$(TYPEDEF)

Retain the observed inputs, separate observed reference, tables, illustration,
and written destinations of one completed report.

Plain and HTML display show the completed quantity tables. `artifact[R]`
retrieves the reported resistance DataFrame, or an ordered vector of DataFrames
for a collection. `artifact[i, R]` retrieves the table for result position `i`.
Lookup and display use completed tables without further observation or calculation.

$(TYPEDFIELDS)
"""
struct ReportArtifact{T,I,O}
    "One atomic observation or an ordinary vector of observations."
    observed::Union{ObservedResult,Vector{ObservedResult}}
    "Separate observed reference, when supplied."
    reference::Union{Nothing,ObservedResult}
    "Quantity-wise tables and other retained scientific summaries."
    tables::T
    "Optional rendered illustration."
    illustration::I
    "Written destinations, or nothing for an in-memory report."
    output::O
end

# Only ordinary report containers are traversed; DataFrame cells are leaves.
_reported_leaves(table::DataFrame) = (table,)
_reported_leaves(tables::Union{NamedTuple,Tuple,AbstractVector}) =
    Iterators.flatten(_reported_leaves(child) for child in tables)
_reported_leaves(_) = ()

function _reported_tables(artifact::ReportArtifact)
    output=artifact.tables
    benchmark=output isa NamedTuple && haskey(output,:features)
    benchmark && (output=output.quantities)
    collection=artifact.observed isa AbstractVector
    # A specialized definition may return one aggregate table. Its display is
    # still native, but a collection does not associate it with one gridpoint.
    output isa DataFrame && return [(point_index=collection ? nothing : 1,tables=[output])]
    return map(eachindex(_observed_points(artifact.observed))) do index
        tables=benchmark || collection ? output[index] : output
        (point_index=index,tables=collect(DataFrame,_reported_leaves(tables)))
    end
end

function _reported_table(artifact,groups,request,index)
    identity=Grammar.normalize_observation_selector(request_identity(request))
    indices=request_indices(request)
    tables=collect(DataFrame,Iterators.flatten(group.tables for group in groups
        if group.point_index==index))
    matches=filter(tables) do table
        retained=metadata(table,"request",nothing)
        retained===nothing && return false
        isequal(Grammar.normalize_observation_selector(request_identity(retained)),identity) &&
            (isempty(indices) || isequal(request_indices(retained),indices))
    end
    length(matches)==1 && return only(matches)
    unassociated=any(group -> group.point_index===nothing,groups)
    available_tables=unassociated ? Iterators.flatten(group.tables for group in groups) : tables
    available=join((repr(metadata(table,"request",nothing)) for table in available_tables),", ")
    problem=isempty(matches) ? "absent" : "ambiguous"
    unassociated && (problem="not associated with a gridpoint")
    id=_observed_points(artifact.observed)[index].gridpoint.id
    location="reported result $index"*(id===nothing ? "" : " (gridpoint $(repr(id)))")
    throw(ArgumentError("reported request $(repr(request)) is $problem for $location; " *
        "available reported requests: [$available] (nothing means absent quantity descriptors). " *
        "Use report(...; values=...) for a different selection."))
end

"""
$(TYPEDSIGNATURES)

Retrieve an already-produced quantity DataFrame using an observation request.
An atomic report returns one DataFrame; a collection returns an ordered vector,
including for one result. `artifact[i, request]` returns the table for result
position `i`, independently of its recorded scientific gridpoint identifier.

An unindexed request returns the product with its reported selection intact.
An indexed request must match that original selection exactly. Complete
transformation and statistical identities remain distinct. Lookup never
acquires, slices, converts, or tabulates quantities.

# Arguments

- `artifact`: A completed report.
- `request`: A quantity selector or complete observation request.

# Returns

The stored DataFrame itself, or a vector of those DataFrames. Use `copy` for an
independent table; editing a returned table does not edit retained observations.

# Errors

Absent or ambiguous products, missing descriptors, and Boolean point indices
raise `ArgumentError`. Invalid result positions raise `BoundsError`.

# Examples

```julia
resistance = constants_report[R]
resistance_tables = collection_report[R]
second_resistance = collection_report[2, R]
```
"""
function Base.getindex(artifact::ReportArtifact,request)
    groups=_reported_tables(artifact)
    tables=map(eachindex(_observed_points(artifact.observed))) do index
        _reported_table(artifact,groups,request,index)
    end
    return artifact.observed isa ObservedResult ? only(tables) : tables
end

function Base.getindex(artifact::ReportArtifact,index::Integer,request)
    index isa Bool && throw(ArgumentError("result position must be an integer, not Bool"))
    points=_observed_points(artifact.observed)
    index in eachindex(points) || throw(BoundsError(artifact,(index,request)))
    return _reported_table(artifact,_reported_tables(artifact),request,index)
end

"""
$(TYPEDEF)

Select retained quantities for separate tables. Raw-input conveniences first
construct [`ObservedResult`](@ref). Illustrations consume those same observations.

$(TYPEDFIELDS)
"""
struct TableReportDefinition{R<:Tuple,U<:Tuple,P,O<:NamedTuple} <: AbstractReportDefinition
    "Quantity requests; an empty tuple selects all retained quantities."
    requests::R
    "Display units used only when constructing a raw-input observation."
    units::U
    "True, a plotting callable, or nothing."
    illustration::P
    "Options passed to the illustration call."
    plot_options::O
    "Engineering recentering used only when constructing a raw-input observation."
    clip::Bool
end
TableReportDefinition(requests::Tuple=();units::Tuple=(),illustration=nothing,
    plot_options::NamedTuple=(;),clip::Bool=true) =
    TableReportDefinition(requests,units,illustration,plot_options,clip)

"""Select retained products without extracting or recomputing numerical values."""
function select end
"""Build quantity-wise tables from detached observations."""
function tabulate end
"""Render an optional illustration from detached observations."""
function illustrate end
"""Encode observed tables for a requested output format."""
function encode end
"""Write already encoded report output."""
function write end

_observed_points(observed::ObservedResult) = (observed,)
_observed_points(observed::AbstractVector{<:ObservedResult}) = observed

select(definition::AbstractReportDefinition,observed::AbstractVector{<:ObservedResult};reference=nothing) =
    map(point -> select(definition,point;reference),observed)

function select(definition::TableReportDefinition,observed::ObservedResult;reference=nothing)
    isempty(definition.requests) && return observed.quantities
    return [Grammar.observation_product(observed,request)
        for request in Grammar.observation_requests(observed,definition.requests).retained]
end

tabulate(definition::AbstractReportDefinition,observed;reference=nothing) =
    tabulate(definition,observed,select(definition,observed;reference);reference)

illustrate(::AbstractReportDefinition,observed,tables;reference=nothing) = nothing
encode(::AbstractReportDefinition,observed,tables,illustration;reference=nothing) = nothing
write(::AbstractReportDefinition,::Nothing) = nothing

function illustrate(definition::TableReportDefinition,observed,tables;reference=nothing)
    illustration=definition.illustration
    (illustration===nothing || illustration===false) && return nothing
    options=merge(definition.plot_options,isempty(definition.requests) ? (;) : (ydata=definition.requests,))
    callable=illustration===true ? PlotBuilder.plot : illustration
    return reference===nothing ? callable(observed;options...) : callable(observed;reference,options...)
end

"""
$(TYPEDSIGNATURES)

Select retained products, build tables, render an optional illustration, encode
output, and write artifacts, in that order. Only observations enter these report
stages. `reference` is an atomic
observation kept outside the reported-result collection.
"""
function report(definition::AbstractReportDefinition,
        observed::Union{ObservedResult,AbstractVector{<:ObservedResult}};
        reference::Union{Nothing,ObservedResult}=nothing)
    selected=select(definition,observed;reference)
    tables=tabulate(definition,observed,selected;reference)
    illustration=illustrate(definition,observed,tables;reference)
    encoded=encode(definition,observed,tables,illustration;reference)
    written=write(definition,encoded)
    points=observed isa ObservedResult ? observed : collect(ObservedResult,observed)
    return ReportArtifact(points,reference,tables,illustration,written)
end

function report(definition::TableReportDefinition,source;kwargs...)
    return report(source;values=definition.requests,units=definition.units,clip=definition.clip,
        illustration=definition.illustration,plot_options=definition.plot_options,kwargs...)
end

"""Return the quantity and unit metadata of an observed table's columns."""
observation_columns(table::DataFrame) = metadata(table,"observation_columns")

# Resolve the Julia dispatch intersection between a generic raw convenience and
# the common observed workflow without adding a second execution path.
function report(definition::TableReportDefinition,
        observed::Union{ObservedResult,AbstractVector{<:ObservedResult}};reference=nothing)
    return invoke(report,Tuple{AbstractReportDefinition,typeof(observed)},definition,observed;reference)
end

# These are acquisition keywords, not a second set of observation defaults.
function _report_observation_options(options; retained = false)
    for key in keys(options)
        key in (:ydata, :rdata, :requests, :quantities) && throw(ArgumentError(
            "use values to select report quantities; $key is not a reporting keyword"))
        key in (:units, :length_unit, :quantity_units, :frequency_unit, :freq_unit,
            :clip, :atol, :frequencies) || throw(ArgumentError(
            "unknown reporting keyword $key; illustration options belong in plot_options"))
        retained && key in (:clip, :atol, :frequencies) &&
            throw(ArgumentError(
                "$key requires a raw result; retained reports only select or re-express recorded values"))
    end
    haskey(options, :frequency_unit) && haskey(options, :freq_unit) &&
        throw(ArgumentError(
            "use frequency_unit or freq_unit, not both"))
    return (;
        (key===:freq_unit ? :frequency_unit=>value : key=>value
    for (key, value) in options)...)
end

function _report_plot_options(source, requests, illustration, options::NamedTuple)
    isempty(options) || !(illustration===nothing || illustration===false) ||
        throw(ArgumentError(
            "plot_options requires an explicit illustration"))
    for key in (:values, :rdata, :requests, :quantities)
        haskey(options, key) &&
            throw(ArgumentError("select report and illustration quantities with values"))
    end
    if haskey(options, :ydata)
        selected=Grammar.observation_selection(source, options.ydata)
        expected=Grammar.observation_requests(source, requests; complete_pairs = true).displayed
        Grammar.observation_requests(source, selected; complete_pairs = true).displayed==expected ||
            throw(ArgumentError("plot_options.ydata conflicts with the report's values selection"))
    end
    return (; (key=>value for (key, value) in pairs(options) if key!==:ydata)...)
end

"""
$(TYPEDSIGNATURES)

Build separate in-memory quantity tables from a completed result or collection.
The same scientific selection used by plotting's `ydata` is named `values` here.
Raw results are observed first; retained observations preserve their recorded
units and numerical eligibility unless compatible display units are requested.

# Arguments

- `source`: A completed primary result, standalone Z/Y tensor, result space,
  ordinary collection, or retained observation.
- `selection`: Optional positional alternative to `values`; do not supply both.

# Keywords

- `values=nothing`: A selector, one `@observe` request, or a tuple of requests.
  `nothing` and `()` use the source owner's defaults, or all retained products.
- `units`, `length_unit`, `quantity_units`, `frequency_unit`: Observation unit
  options. Raw defaults belong to the observation owner; omitted retained
  options preserve recorded units. `freq_unit` is an alternative spelling of
  `frequency_unit`; supplying both is an error.
- `clip`, `atol`, `frequencies`: Raw observation options. Cutoffs use native
  units; standalone tensor frequency context is in \\[Hz\\]. These keywords
  cannot be supplied for retained inputs.
- `reference=nothing`: A separate atomic raw or observed reference. It does not
  trigger a numerical comparison or become another reported result.
- `illustration=nothing`: `true` or a plotting callable requests an illustration
  of the prepared observations with the matching `ydata` selection.
- `plot_options=(;)`: Options for an explicitly requested illustration.

# Returns

- A [`ReportArtifact`](@ref) containing observed inputs and separate quantity
  DataFrames. Default reporting creates no figure and writes no files.

# Errors

Unknown options, competing selection keywords, conflicting illustration
selections, and acquisition options on retained inputs raise `ArgumentError`.
Specialized definitions remain available through `report(definition, observed)`.

# Examples

```julia
report(constants)
report(constants; values=(R, L, G, C), length_unit=:kilo,
    quantity_units=(R=:base, L=:milli, G=:micro, C=:micro))
report(line_parameters; values=@observe(R[1, 1, 1:12]))
```
"""
function report(
        source::Union{Grammar.AbstractCoreResult, Grammar.AbstractResultSpace,
            Engine.SeriesImpedance, Engine.ShuntAdmittance, AbstractVector, Tuple};
        values = nothing, reference = nothing, illustration = nothing,
        plot_options::NamedTuple = (;), kwargs...)
    collection=source isa Union{AbstractVector, Tuple, Grammar.AbstractParametricResult}
    collection && isempty(source) &&
        throw(ArgumentError("report requires at least one result"))
    acquisition=_report_observation_options(kwargs;
        retained = collection && any(point -> point isa ObservedResult, source))
    point=collection ? first(source) : source
    requests=Grammar.observation_selection(point, values)
    displayed=isempty(requests) ? () :
              Grammar.observation_requests(point, requests; complete_pairs = true).displayed
    options=_report_plot_options(point, displayed, illustration, plot_options)
    reference===nothing ||
        reference isa Union{ObservedResult, Grammar.AbstractCoreResult,
            Engine.SeriesImpedance, Engine.ShuntAdmittance} ||
        throw(ArgumentError(
            "reference must be an atomic result or ObservedResult"))
    observed=observables(source, requests; complete_pairs = true, acquisition...)
    if reference isa ObservedResult
        display_units=(;
            (key=>value
        for (key, value) in pairs(acquisition)
        if key in (:units, :length_unit, :quantity_units, :frequency_unit))...)
        isempty(display_units) || (reference=ObservedResult(reference; display_units...))
    elseif reference!==nothing
        reference=ObservedResult(reference, requests; complete_pairs = true, acquisition...)
    end
    return report(
        observed; values = displayed, reference, illustration, plot_options = options)
end

function report(observed::Union{ObservedResult, AbstractVector{<:ObservedResult}};
        values = nothing, reference::Union{Nothing, ObservedResult} = nothing,
        illustration = nothing, plot_options::NamedTuple = (;), kwargs...)
    display_units=_report_observation_options(kwargs; retained = true)
    observed isa AbstractVector && isempty(observed) &&
        throw(ArgumentError("report requires at least one observation"))
    point=observed isa ObservedResult ? observed : first(observed)
    requests=Grammar.observation_selection(point, values)
    displayed=Grammar.observation_requests(point, requests).displayed
    options=_report_plot_options(point, displayed, illustration, plot_options)
    if !isempty(display_units)
        observed=observables(observed, requests; display_units...)
        reference===nothing ||
            (reference=ObservedResult(reference, requests; display_units...))
    end
    definition=TableReportDefinition(requests; illustration,
        plot_options = illustration===nothing || illustration===false ? options :
                       merge(options, (ydata = displayed,)))
    return report(definition, observed; reference)
end

function report(
        source::Union{Grammar.AbstractCoreResult, Grammar.AbstractResultSpace,
            Engine.SeriesImpedance, Engine.ShuntAdmittance, ObservedResult, AbstractVector, Tuple},
        selection; kwargs...)
    haskey(kwargs, :values) &&
        throw(ArgumentError("use positional selection or values, not both"))
    return report(source; values = selection, kwargs...)
end
