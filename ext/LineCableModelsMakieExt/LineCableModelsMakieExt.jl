"""
    LineCableModelsMakieExt

Add compact high-level LineCableModels plotting methods to native Makie.
"""
module LineCableModelsMakieExt

import LineCableModels
using LineCableModels: EarthLayer, LineParameters, Material, RadialDielectric,
    SeriesImpedance, ShuntAdmittance, Y, Z, label, nominal,
    observables, outer_radius
import Makie
using Makie: Auto, Axis, Button, Colorbar, DataAspect, Figure,
    Fixed, GridLayout, Label, Legend, LineElement, Mixed,
    Observable, Outside, Rect2f, Relative, Theme, Toggle,
    colgap!, colsize!, content, errorbars!, height,
    hlines!, hspan!, lift, lines!,
    off, on, onany, poly!, reset_limits!, rowgap!,
    rowsize!, scatter!, stairs!, text!, to_value, translate!,
    update!, width, widths, with_theme
using LinearAlgebra: diag
using Printf: @sprintf
using Colors: HSV, Oklab, RGB, RGBA, blue, green, red
import Dates
using Statistics: mean
import Statistics

import LineCableModels.Units
import LineCableModels.Engine
import LineCableModels.DataModel
import LineCableModels.Grammar
import LineCableModels.ImportExport
import LineCableModels.UQ
import Makie.GridLayoutBase
import LineCableModels.Grammar:
                                request_identity,
                                request_indices,
                                request_quantity

function current_backend_symbol()
    backend = Makie.current_backend()
    backend isa Module || return :none
    name = nameof(backend)
    name === :CairoMakie && return :cairo
    name === :GLMakie && return :gl
    name === :WGLMakie && return :wgl
    return :unknown
end

include("recipes/line_data.jl")
include("recipes/comparison_data.jl")
include("material_colors.jl")
include("recipes/preview_types.jl")
include("recipes/preview_data.jl")
include("attributes.jl")
include("series_styles.jl")
include("shell.jl")
include("recipes/line_facets.jl")
include("recipes/preview_render.jl")
include("montecarlo.jl")
include("native_export.jl")

import LineCableModels.PlotBuilder: plot, preview, show_material_scale
include("recipes/formulation_comparisons.jl")

const _LineSource = Union{LineCableModels.LineParameters,LineCableModels.SeriesImpedance,LineCableModels.ShuntAdmittance}

function _plot_ydata(positional,keyword,default)
    keyword===nothing && return positional===nothing ? default : positional
    positional===nothing || throw(ArgumentError("use positional ydata or the ydata keyword"))
    return keyword
end
function _plot_requests(source,selection)
    selection===nothing && return ()
    selection isa Function && return (selection,)
    selection isa Tuple || throw(ArgumentError("ydata must be a selector or tuple of requests"))
    isempty(selection) && return ()
    if first(selection) isa Function
        identity=request_identity(selection)
        declared=source isa Grammar.ObservedResult ? Tuple(request_identity(q.request) for q in source.quantities) : observables(typeof(source))
        identity in declared && (identity isa Tuple || !isempty(request_indices(selection))) && return (selection,)
    end
    return selection
end

function plot(source::_LineSource,selection=nothing;ydata=nothing,frequencies=nothing,
        freq_unit=:base,length_unit=:kilo,quantity_units=nothing,clip::Bool=true,atol=nothing,kwargs...)
    requests=_plot_requests(source,_plot_ydata(selection,ydata,()))
    normalized=Grammar.observation_requests(source,requests;complete_pairs=true)
    observed=Grammar.ObservedResult(source,requests;complete_pairs=true,frequencies,
        frequency_unit=freq_unit,length_unit,quantity_units,clip,atol)
    return plot(observed;ydata=normalized.displayed,kwargs...)
end
function plot(source::Union{SeriesImpedance,ShuntAdmittance},f::AbstractVector,selection=nothing;kwargs...)
    return plot(source,selection;frequencies=f,kwargs...)
end
Makie.plot(source::_LineSource,args...;kwargs...) = plot(source,args...;kwargs...)

function plot(observed::Grammar.ObservedResult,selection=nothing;ydata=nothing,reference=nothing,kwargs...)
    return plot([observed],selection;ydata,reference,kwargs...)
end
function plot(observed::AbstractVector{<:Grammar.ObservedResult},selection=nothing;
        ydata=nothing,reference::Union{Nothing,Grammar.ObservedResult}=nothing,
        series_labels=nothing,series_attributes=nothing,problem=nothing,formulations=nothing,band=nothing,kwargs...)
    isempty(observed) && throw(ArgumentError("plot requires at least one observed point"))
    for key in (:clip,:atol,:freq_unit,:length_unit,:quantity_units,:frequencies)
        haskey(kwargs,key) && throw(ArgumentError("$key belongs to observation construction; this plot selects retained products"))
    end
    selected=findall(observed) do point
        id=get(point.gridpoint,:id,nothing)
        (problem===nothing || id!==nothing && id.problem_index in (problem isa Integer ? (problem,) : problem)) &&
            (formulations===nothing || id!==nothing && id.formulation_index in (formulations isa Integer ? (formulations,) : formulations))
    end
    isempty(selected) && throw(ArgumentError("no retained observations match the original point selection"))
    if formulations isa Union{AbstractVector,Tuple}
        sort!(selected;by=index -> findfirst(==(observed[index].gridpoint.id.formulation_index),formulations))
    end
    candidates=observed[selected]
    requests=_plot_requests(first(candidates),_plot_ydata(selection,ydata,()))
    requests=Grammar.observation_requests(first(candidates),requests).displayed
    if band!==nothing
        records=filter(error -> isequal(error.band,band),first(candidates).errors)
        isempty(records) && throw(ArgumentError("the requested comparison band was not retained"))
        samples=first(records).settings.indices
        all(error -> error.settings.indices==samples,records) || throw(ArgumentError("retained comparisons disagree on band coordinates"))
        requests=map(requests) do request
            product=Grammar.observation_product(first(candidates),request)
            c=product.coordinates
            selected_samples=intersect(c.samples,samples)
            isempty(selected_samples) && throw(ArgumentError("the requested band contains no retained samples"))
            identity=Grammar.request_identity(request)
            prefix=identity isa Tuple ? identity : (identity,)
            indices=c.kind===:matrix ? (c.rows,c.columns,selected_samples) : (c.rows,selected_samples)
            (prefix...,indices...)
        end
    end
    points=reference===nothing ? Tuple(candidates) : (Tuple(candidates)...,reference)
    count=length(observed)+(reference===nothing ? 0 : 1)
    displayed=reference===nothing ? selected : [selected;count]
    labels=series_labels===nothing ? nothing :
        Tuple(_comparison_labels(series_labels,count)[index] for index in displayed)
    attributes=series_attributes isa Union{Tuple,AbstractVector} ?
        Tuple(_series_attributes(series_attributes,count)[index] for index in displayed) : series_attributes
    roles=reference===nothing ? nothing : [fill(:candidate,length(candidates));:reference]
    indices=reference===nothing ? selected : [selected.+1;1]
    return _addon_line_pages(points;ydata=requests,series_labels=labels,formulation_roles=roles,
        series_indices=indices,series_count=count,series_attributes=attributes,kwargs...)
end
Makie.plot(source::Union{Grammar.ObservedResult,AbstractVector{<:Grammar.ObservedResult}},args...;kwargs...) = plot(source,args...;kwargs...)

function plot(first::LineParameters,second::LineParameters,rest...;ydata=nothing,freq_unit=:base,length_unit=:kilo,quantity_units=nothing,clip=true,atol=nothing,kwargs...)
    sources=LineParameters[first,second]
    selection=nothing
    for item in rest
        item isa LineParameters && selection===nothing ? push!(sources,item) :
            selection===nothing ? (selection=item) : throw(ArgumentError("one trailing ydata selection is accepted"))
    end
    requests=_plot_requests(first,_plot_ydata(selection,ydata,()))
    normalized=Grammar.observation_requests(first,requests;complete_pairs=true)
    observed=observables(sources,requests;complete_pairs=true,frequency_unit=freq_unit,length_unit,quantity_units,clip,atol)
    return plot(observed;ydata=normalized.displayed,kwargs...)
end
Makie.plot(first::LineParameters,second::LineParameters,rest...;kwargs...) = plot(first,second,rest...;kwargs...)

function plot(sources::NamedTuple,selection=nothing;ydata=nothing,series_labels=nothing,
        freq_unit=:base,length_unit=:kilo,quantity_units=nothing,clip=true,atol=nothing,kwargs...)
    raw=collect(values(sources))
    isempty(raw) && throw(ArgumentError("plot requires at least one result"))
    requests=_plot_requests(first(raw),_plot_ydata(selection,ydata,()))
    normalized=Grammar.observation_requests(first(raw),requests;complete_pairs=true)
    observed=observables(raw,requests;complete_pairs=true,frequency_unit=freq_unit,length_unit,quantity_units,clip,atol)
    labels=series_labels===nothing ? Tuple(string.(keys(sources))) : series_labels
    return plot(observed;ydata=normalized.displayed,series_labels=labels,kwargs...)
end
Makie.plot(sources::NamedTuple,args...;kwargs...) = plot(sources,args...;kwargs...)

function preview(
        design::DataModel.CableDesign;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        design;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

# The public extension method consumes backend/display choices and forwards the
# remaining preview options unchanged. DataModel retains only detached geometry
# and material attributes; Makie objects and backend state stay here.
function preview(
        designs::AbstractVector{<:DataModel.CableDesign};
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        designs;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function preview(
        system::DataModel.LineCableSystem;
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_preview(
        system;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

function show_material_scale(
        ; backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        kwargs...
)
    return _addon_material_scale(;
        backend,
        display_plot,
        controls,
        kwargs...
    )
end

end # module LineCableModelsMakieExt
