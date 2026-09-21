"""
$(TYPEDEF)

Supertype for definitions consumed by [`report`](@ref).
"""
abstract type AbstractReportDefinition end

"""
$(TYPEDEF)

Retain the observed inputs, separate observed reference, tables, illustration,
and written destinations of one completed report.

$(TYPEDFIELDS)
"""
struct ReportArtifact{T,I,O}
    "One atomic observation or an ordinary vector of observations."
    observed::Union{ObservedResult,Vector{ObservedResult}}
    "Separate observed reference, if this is a comparison report."
    reference::Union{Nothing,ObservedResult}
    "Quantity-wise tables and other retained scientific summaries."
    tables::T
    "Optional rendered illustration."
    illustration::I
    "Written destinations, or nothing for an in-memory report."
    output::O
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

function select(definition::TableReportDefinition,observed::ObservedResult)
    isempty(definition.requests) && return observed.quantities
    return [_selected_quantity(observed,request)
        for request in Grammar.observation_requests(observed,definition.requests).retained]
end

_selected_quantity(observed::ObservedResult,request) = Grammar.observation_product(observed,request)

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

Build tables, an optional illustration, encoded output, and written artifacts,
in that order. Only observations cross this boundary. `reference` is an atomic
observation kept outside the candidate collection.
"""
function report(definition::AbstractReportDefinition,
        observed::Union{ObservedResult,AbstractVector{<:ObservedResult}};
        reference::Union{Nothing,ObservedResult}=nothing)
    tables=tabulate(definition,observed;reference)
    illustration=illustrate(definition,observed,tables;reference)
    encoded=encode(definition,observed,tables,illustration;reference)
    written=write(definition,encoded)
    points=observed isa ObservedResult ? observed : collect(ObservedResult,observed)
    return ReportArtifact(points,reference,tables,illustration,written)
end

function report(definition::TableReportDefinition,source;kwargs...)
    observed=observables(source,definition.requests;units=definition.units,clip=definition.clip,complete_pairs=true,kwargs...)
    return report(definition,observed)
end

"""Return the quantity and unit metadata of an observed table's columns."""
observation_columns(table::DataFrame) = metadata(table,"observation_columns")

# Resolve the Julia dispatch intersection between a generic raw convenience and
# the common observed workflow without adding a second execution path.
function report(definition::TableReportDefinition,
        observed::Union{ObservedResult,AbstractVector{<:ObservedResult}};reference=nothing)
    return invoke(report,Tuple{AbstractReportDefinition,typeof(observed)},definition,observed;reference)
end
