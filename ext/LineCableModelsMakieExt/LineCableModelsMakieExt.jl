"""
    LineCableModelsMakieExt

Add compact high-level LineCableModels plotting methods to native Makie.
"""
module LineCableModelsMakieExt

import LineCableModels
using LineCableModels: EarthLayer, Material, RadialDielectric,
                       SeriesImpedance, ShuntAdmittance, label, nominal,
                       observables, outer_radius
import Makie
using Makie: Auto, Axis, Button, Colorbar, DataAspect, Figure,
             Block, CategoricalConversion, Fixed, GridLayout, Label, Legend, LineElement,
             Mixed,
             Observable, Outside, Rect2f, Relative, RichText, Theme, Toggle,
             colgap!, colsize!, content, defaultlimits, errorbars!,
             fast_string_boundingboxes_obs,
             hlines!, hspan!, lift, lines!,
             off, on, onany, poly!, reset_limits!, rowgap!,
             rowsize!, scatter!, stairs!, text!, to_value, translate!, update!, widths,
             with_theme
using Printf: @sprintf
using Colors: HSV, Oklab, RGB, RGBA, blue, green, red
import Dates
import Base: resize!
using Statistics: mean

import LineCableModels.Units
import LineCableModels.Engine
import LineCableModels.DataModel
import LineCableModels.Grammar
import LineCableModels.ImportExport
import LineCableModels.UQ
import Makie.GridLayoutBase
using Makie.GridLayoutBase: HorizontalAlignment, VerticalAlignment, firstrow, lastrow
import LineCableModels.Grammar:
                                request_identity

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
include("controls.jl")
include("guides.jl")
include("layout.jl")
include("recipes/line_facets.jl")
include("recipes/preview_render.jl")
include("montecarlo.jl")
include("export_presentation.jl")
include("native_export.jl")

import LineCableModels.PlotBuilder: plot, preview, show_material_scale

const _PrimarySource = Union{
    LineCableModels.LineParameters, LineCableModels.SeriesImpedance,
    LineCableModels.ShuntAdmittance, LineCableModels.CableConstants}
const _PlotSource = Union{_PrimarySource, Grammar.ObservedResult}
struct _Omitted end
const _omitted = _Omitted()

function _plot_ydata(positional, keyword, default)
    keyword===nothing && return positional===nothing ? default : positional
    positional===nothing ||
        throw(ArgumentError("use positional ydata or the ydata keyword"))
    return keyword
end

# Split only public keyword ownership. Scientific normalization stays in Grammar.
function _plot_observation_options(kwargs; retained = false)
    haskey(kwargs, :complete_pairs) &&
        throw(ArgumentError("complete_pairs is owned by observation construction"))
    haskey(kwargs, :freq_unit) && haskey(kwargs, :frequency_unit) &&
        throw(ArgumentError("use frequency_unit or freq_unit, not both"))
    keys=(:clip, :atol, :units, :length_unit, :quantity_units,
        :frequency_unit, :frequencies, :freq_unit)
    if retained
        for key in (:clip, :atol, :frequencies)
            haskey(kwargs, key) &&
                throw(ArgumentError("$key requires primary observation construction; retained plotting only selects or re-expresses recorded values"))
        end
    end
    acquisition=(;
        (key===:freq_unit ? :frequency_unit=>value : key=>value
    for (key, value) in kwargs if key in keys)...)
    presentation=(; (key=>value for (key, value) in kwargs if key ∉ keys)...)
    return acquisition, presentation
end

function plot(source::_PrimarySource, selection = nothing;
        ydata = nothing, reference = nothing, kwargs...)
    acquisition, presentation=_plot_observation_options(kwargs)
    requests=Grammar.observation_selection(source, _plot_ydata(selection, ydata, ()))
    normalized=Grammar.observation_requests(source, requests; complete_pairs = true)
    observed=Grammar.ObservedResult(source, requests; complete_pairs = true, acquisition...)
    if reference!==nothing && !(reference isa Grammar.ObservedResult)
        reference isa _PrimarySource ||
            throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Grammar.ObservedResult(reference, requests; complete_pairs = true, acquisition...)
    end
    return plot(observed; ydata = normalized.displayed, reference, presentation...)
end
function plot(source::Union{SeriesImpedance, ShuntAdmittance},
        f::AbstractVector, selection = nothing; kwargs...)
    return plot(source, selection; frequencies = f, kwargs...)
end
Makie.plot(source::_PrimarySource, args...; kwargs...) = plot(source, args...; kwargs...)

function plot(observed::Grammar.ObservedResult, selection = nothing; kwargs...)
    return plot([observed], selection; kwargs...)
end

function plot(observed::AbstractVector{<:Grammar.ObservedResult}, selection = nothing;
        ydata = nothing, reference = nothing,
        series_labels = nothing, series_attributes = nothing, problem = nothing,
        formulations = nothing, band = nothing,
        errorbar_sampling = nothing, title = nothing, title_prefix = nothing, figure_title = nothing,
        title_attributes::NamedTuple = (;), panel_titles = nothing, fig_size = nothing, layout = nothing,
        xscale = nothing, yscale = :linear, legend_position = _omitted, legend_title = nothing,
        legend_attributes::NamedTuple = (;), legend_overflow::Symbol = :show_all, panel_legends = (),
        backend = nothing, display_plot::Bool = true, controls::Bool = true,
        export_theme::Symbol = :default, open_export::Bool = true, kwargs...)
    isempty(observed) && throw(ArgumentError("plot requires at least one observed point"))
    display_units, kwargs=_plot_observation_options(kwargs; retained = true)
    for key in (:blocks, :series_defaults, :series_count, :series_indices, :formulation_roles,
        :dependent_plots, :marker_coordinates, :signed_ylog, :anchor, :legend_anchor)
        haskey(kwargs, key) && throw(ArgumentError("unsupported plotting keyword $key"))
    end
    selected=findall(observed) do point
        id=get(point.gridpoint, :id, nothing)
        (problem===nothing ||
         id!==nothing &&
         id.problem_index in (problem isa Integer ? (problem,) : problem)) &&
            (formulations===nothing ||
             id!==nothing && id.formulation_index in
             (formulations isa Integer ? (formulations,) : formulations))
    end
    isempty(selected) &&
        throw(ArgumentError("no retained observations match the original point selection"))
    if formulations isa Union{AbstractVector, Tuple}
        sort!(selected;
            by = index -> findfirst(==(observed[index].gridpoint.id.formulation_index), formulations))
    end
    candidates=observed[selected]
    requests=Grammar.observation_selection(first(candidates), _plot_ydata(selection, ydata, ()))
    requests=Grammar.observation_requests(first(candidates), requests).displayed
    if reference!==nothing && !(reference isa Grammar.ObservedResult)
        reference isa _PrimarySource ||
            throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Grammar.ObservedResult(reference, requests; complete_pairs = true)
    end
    if !isempty(display_units)
        candidates=[Grammar.ObservedResult(point, requests; display_units...)
                    for point in candidates]
        reference===nothing ||
            (reference=Grammar.ObservedResult(reference, requests; display_units...))
    end
    sources=reference===nothing ? Tuple(candidates) : (Tuple(candidates)..., reference)
    count=length(observed)+(reference===nothing ? 0 : 1)
    displayed=reference===nothing ? selected : [selected; count]
    explicit_labels=series_labels!==nothing
    candidate_labels=explicit_labels && reference!==nothing &&
                     length(series_labels)==length(observed)
    labels=explicit_labels ?
           Tuple(_comparison_labels(series_labels, candidate_labels ? count-1 : count)[i]
    for i in (candidate_labels ? selected : displayed)) : nothing
    attributes=Tuple(_series_attributes(series_attributes, count)[i] for i in displayed)
    roles=reference===nothing ? fill(:candidate, length(candidates)) :
          [fill(:candidate, length(candidates)); :reference]
    slots=reference===nothing ? selected : [selected; 0]
    reference_id=reference===nothing ? nothing : get(reference.gridpoint, :id, nothing)
    products=map(request -> Grammar.observation_product(sources, request; band, reference_id), requests)
    published=map(eachindex(sources)) do index
        _prepare_line_observations(Tuple(records[index] for records in products))
    end
    facets=_semantic_line_facets(published, requests)
    matrix_positions=[[(f.row, f.column)
                       for f in facets if f.kind===:matrix && f.request_index==i]
                      for i in eachindex(requests)]
    filter!(!isempty, matrix_positions)
    flow_counts=[Base.count(f -> f.request_index==i, facets) for i in eachindex(requests)]
    capacity=_addon_capacity(layout, matrix_positions, flow_counts)
    reference_size=_addon_figure_size(fig_size, capacity)
    pages=_semantic_line_pages(facets, capacity; automatic = layout===nothing)
    if panel_titles isa Union{Tuple, AbstractVector}
        length(panel_titles)==length(facets) ||
            throw(DimensionMismatch("panel_titles must contain one title per selected panel before pagination"))
        panel_titles=Dict((f.identity, f.row, f.column)=>label
        for (f, label) in zip(facets, panel_titles))
    end
    identities=Set((f.row, f.column) for f in facets)
    all(pair -> _addon_panel_identity(first(pair)) in identities, _addon_panel_legend_pairs(panel_legends)) ||
        throw(ArgumentError("panel legend identity is absent from the selected panels"))
    # Resolve descriptions and scientific groups once per request, before pages.
    groups=map(request -> Grammar.observation_groups(candidates; request), requests)
    request_labels=map(requests) do request
        result=explicit_labels ? Any[labels...] :
               Grammar.observation_labels(sources; request)
        candidate_labels &&
            push!(result, last(Grammar.observation_labels(sources; request)))
        (!explicit_labels || candidate_labels) && reference!==nothing &&
            (result[end]*=" (reference)")
        result
    end
    legend_overflow in (:ellipsis, :show_all) ||
        throw(ArgumentError("legend_overflow must be :ellipsis or :show_all"))
    errorbar_sampling===nothing || errorbar_sampling in (:staggered, :all) ||
        throw(ArgumentError("errorbar_sampling must be :staggered or :all"))
    _addon_activate_backend(backend)
    built=LineCableModels.UIPlot[]
    for (page_index, page) in enumerate(pages)
        request_index=first(page.facets).request_index
        retained=[group.representative for group in groups[request_index]]
        reference===nothing || push!(retained, length(sources))
        styles=_addon_comparison_styles(Tuple(slots[i] for i in retained),
            Tuple(roles[i] for i in retained), length(observed))
        position=legend_position===_omitted ?
                 (length(retained)>1 || explicit_labels ? :bottom : nothing) :
                 legend_position
        sampling=errorbar_sampling===nothing ? (length(retained)>1 ? :staggered : :all) :
                 errorbar_sampling
        automatic_title=_semantic_page_title(first(sources), page)
        page_title=title===nothing ? automatic_title : title
        title!==nothing && length(pages)>1 && (page_title="$title — $automatic_title")
        title===nothing && title_prefix!==nothing && !ismissing(title_prefix) &&
            (page_title="$title_prefix — $automatic_title")
        visible_title=_semantic_page_option(figure_title, page_index, length(pages), "figure_title")
        push!(built,
            with_theme(_addon_theme(export_theme = export_theme)) do
                _addon_semantic_line_page(
                    first(sources), Tuple(published[i] for i in retained),
                    Tuple(request_labels[request_index][i] for i in retained), page;
                    series_indices = Tuple(slots[i] for i in retained), series_count = length(observed),
                    errorbar_sampling = sampling, series_defaults = styles,
                    series_attributes = Tuple(attributes[i] for i in retained), title = page_title,
                    figure_title = visible_title, title_attributes, panel_titles,
                    fig_size = reference_size, xscale, yscale,
                    legend_position = position, legend_title, legend_attributes, legend_overflow, panel_legends,
                    controls, display_plot = false, export_theme, open_export, kwargs...)
            end)
        last(built).addon_state=merge(last(built).addon_state,
            (observed = sources,
                display_groups = ((
                    request = requests[request_index], groups = groups[request_index]),),
                displayed_indices = retained, nominal_capacity = capacity))
    end
    _addon_calibrate_frames!(built, capacity)
    if layout===nothing
        for (p, page) in zip(built, pages)
            first(page.facets).kind===:diagonal && _addon_responsive_axis_grid!(p)
        end
    end
    if display_plot
        foreach(page -> _addon_display!(page.figure, page.export_name), built)
    end
    return length(built)==1 ? only(built) : built
end
function Makie.plot(source::Grammar.ObservedResult, args...; kwargs...)
    plot(source, args...; kwargs...)
end
function Makie.plot(source::AbstractVector{<:Grammar.ObservedResult}, args...; kwargs...)
    plot(source, args...; kwargs...)
end

function plot(sources::Union{AbstractVector{<:_PlotSource}, Tuple{Vararg{_PlotSource}}},
        selection = nothing;
        ydata = nothing, reference = nothing, kwargs...)
    isempty(sources) && throw(ArgumentError("plot requires at least one result"))
    if all(source -> source isa Grammar.ObservedResult, sources)
        return plot(
            Grammar.ObservedResult[sources...], selection; ydata, reference, kwargs...)
    end
    acquisition, presentation=_plot_observation_options(kwargs)
    requests=Grammar.observation_selection(first(sources), _plot_ydata(selection, ydata, ()))
    normalized=Grammar.observation_requests(first(sources), requests; complete_pairs = true)
    retained_units=(;
        (key=>value
    for (key, value) in pairs(acquisition)
    if key in
       (:units, :length_unit, :quantity_units, :frequency_unit))...)
    observed=Grammar.ObservedResult[source isa Grammar.ObservedResult ?
                                    (isempty(retained_units) ? source :
                                     Grammar.ObservedResult(source, normalized.displayed; retained_units...)) :
                                    Grammar.ObservedResult(source, requests;
                                        complete_pairs = true, acquisition...)
                                    for source in sources]
    if reference!==nothing && !(reference isa Grammar.ObservedResult)
        reference isa _PrimarySource ||
            throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Grammar.ObservedResult(reference, requests; complete_pairs = true, acquisition...)
    end
    return plot(observed; ydata = normalized.displayed, reference, presentation...)
end
function Makie.plot(
        sources::Union{AbstractVector{<:_PlotSource}, Tuple{Vararg{_PlotSource}}},
        args...; kwargs...)
    plot(sources, args...; kwargs...)
end

# Public package plotting can validate heterogeneous owned collections. Native
# Makie array methods remain untouched for ordinary arrays.
function plot(sources::AbstractVector, selection = nothing; kwargs...)
    all(source -> source isa _PlotSource, sources) || throw(ArgumentError(
        "result collections must contain supported primary or observed results"))
    return plot(_PlotSource[sources...], selection; kwargs...)
end

function plot(
        first::_PrimarySource, second::_PrimarySource, rest...; ydata = nothing, kwargs...)
    sources=_PrimarySource[first, second]
    selection=nothing
    for item in rest
        item isa _PrimarySource && selection===nothing ? push!(sources, item) :
        selection===nothing ? (selection=item) :
        throw(ArgumentError("one trailing ydata selection is accepted"))
    end
    return plot(sources; ydata = _plot_ydata(selection, ydata, ()), kwargs...)
end
function Makie.plot(first::_PrimarySource, second::_PrimarySource, rest...; kwargs...)
    plot(first, second, rest...; kwargs...)
end

function plot(sources::NamedTuple, selection = nothing;
        ydata = nothing, series_labels = nothing, kwargs...)
    all(source -> source isa _PlotSource, values(sources)) ||
        throw(ArgumentError("named result plotting requires supported primary or observed members"))
    labels=series_labels===nothing ? Tuple(string.(keys(sources))) : series_labels
    return plot(Tuple(values(sources)), selection; ydata, series_labels = labels, kwargs...)
end
function Makie.plot(sources::NamedTuple{K, <:Tuple{Vararg{_PlotSource}}}, args...; kwargs...) where {K}
    plot(sources, args...; kwargs...)
end

include("recipes/formulation_comparisons.jl")

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
