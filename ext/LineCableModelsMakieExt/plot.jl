const _PrimarySource = Union{
    LineCableModels.LineParameters, LineCableModels.SeriesImpedance,
    LineCableModels.ShuntAdmittance, LineCableModels.CableConstants,
    LineCableModels.PropagationParameters}
const _PlotSource = Union{_PrimarySource, Grammar.ObservedResult}

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
    return plot((source,), selection; ydata, reference, kwargs...)
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
        overlay::Symbol = :auto,
        xscale = nothing, yscale = :linear, legend_position = _omitted, legend_title = nothing,
        legend_attributes::NamedTuple = (;), legend_cap = 0.5, panel_legends = (),
        backend = nothing, display_plot::Bool = true, controls::Bool = true,
        export_theme::Symbol = :default, open_export::Bool = true, kwargs...)
    isempty(observed) && throw(ArgumentError("plot requires at least one observed point"))
    overlay in (:auto, :gridpoints, :coordinates, :rows) || throw(ArgumentError(
        "overlay must be :auto, :gridpoints, :coordinates, or :rows"))
    display_units, kwargs=_plot_observation_options(kwargs; retained = true)
    for key in (:legend_overflow, :blocks, :series_defaults,
        :series_count, :series_indices, :formulation_roles,
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
    input_point_count=length(observed)+(reference===nothing ? 0 : 1)
    displayed=reference===nothing ? selected : [selected; input_point_count]
    slots=reference===nothing ? selected : [selected; 0]
    reference_id=reference===nothing ? nothing : get(reference.gridpoint, :id, nothing)
    products=map(request -> Grammar.observation_product(sources, request; band, reference_id), requests)
    published=map(eachindex(sources)) do index
        _prepare_line_observations(Tuple(records[index] for records in products))
    end
    # The selected cardinality includes an explicit reference and precedes
    # equivalence grouping. The grouping only changes gridpoint curves.
    orientations=map(eachindex(requests)) do index
        kind=first(published).observations[index].coordinates.kind
        overlay===:rows && kind!==:matrix && throw(ArgumentError(
            "overlay=:rows requires matrix observations; select matrix quantities in a separate call"))
        overlay===:auto ?
        (kind===:vector && length(sources)==1 &&
         (layout===nothing || layout==(1, 1)) ? :coordinates : :gridpoints) : overlay
    end
    groups=map(request -> Grammar.observation_groups(candidates; request), requests)
    retained=map(eachindex(requests)) do index
        indices=orientations[index]!==:gridpoints ? collect(eachindex(sources)) :
                [group.representative for group in groups[index]]
        orientations[index]===:gridpoints && reference!==nothing &&
            push!(indices, length(sources))
        indices
    end
    coordinate_identities=map(eachindex(requests)) do index
        first(published).observations[index].coordinates.kind===:assemblies ?
        unique(vcat((source.observations[index].coordinates.labels[
                         source.observations[index].coordinates.assemblies]
        for source in published)...)) : nothing
    end
    coordinate_labels=map(eachindex(requests)) do index
        coordinate_identities[index]===nothing ?
        _coordinate_labels(first(published).observations[index],
            first(published).coordinates[index], orientations[index]) :
        unique(vcat((_coordinate_labels(source.observations[index],
                         source.coordinates[index], orientations[index]) for source in published)...))
    end
    # Positional overrides bind to the overlaid dimension. Preflight every
    # family before constructing the first native figure.
    if series_labels isa Union{Tuple, AbstractVector} ||
       series_attributes isa Union{Tuple, AbstractVector}
        counts=unique(orientations[index]!==:gridpoints ?
                      length(coordinate_labels[index]) :
                      input_point_count for index in eachindex(requests))
        length(counts)==1 || throw(ArgumentError(
            "positional series overrides span figure families with different overlaid dimensions; use separate plot calls"))
    end
    explicit_labels=series_labels!==nothing
    candidate_labels=explicit_labels && reference!==nothing &&
                     all(==(:gridpoints), orientations) &&
                     length(series_labels)==length(observed)
    labels=explicit_labels && any(==(:gridpoints), orientations) ?
           Tuple(_comparison_labels(series_labels, candidate_labels ? input_point_count-1 :
                                                   input_point_count)[i]
    for i in (candidate_labels ? selected : displayed)) : nothing
    attributes=any(==(:gridpoints), orientations) ?
               _series_attributes(series_attributes, input_point_count)[displayed] : nothing
    coordinate_attributes=map(eachindex(requests)) do index
        orientations[index]!==:gridpoints ?
        _series_attributes(series_attributes, length(coordinate_labels[index])) : nothing
    end
    request_labels=map(eachindex(requests)) do index
        request=requests[index]
        if orientations[index]!==:gridpoints
            result=Grammar.observation_labels(sources; request, fallback = "")
            reference!==nothing &&
                (result[end]*=isempty(result[end]) ? "Reference" : " (reference)")
            result
        else
            result=explicit_labels ? Any[labels...] :
                   Grammar.observation_labels(sources; request,
                fallback = length(sources)==1 ? "" : nothing)
            candidate_labels && push!(result,
                last(Grammar.observation_labels(sources; request)))
            (!explicit_labels || candidate_labels) && reference!==nothing &&
                (result[end]*=" (reference)")
            result
        end
    end
    if explicit_labels
        for index in eachindex(requests)
            orientations[index]===:gridpoints && continue
            coordinate_labels[index]=collect(_comparison_labels(series_labels,
                length(coordinate_labels[index])))
        end
    end
    facets=_addon_observation_facets(
        published, requests, orientations, retained, slots, input_point_count,
        request_labels, coordinate_labels, coordinate_identities; explicit_labels)
    capacities=map(eachindex(requests)) do index
        selected_facets=filter(f -> f.request_index==index, facets)
        matrix=first(selected_facets).kind===:matrix && orientations[index]===:gridpoints ?
               [[(f.row, f.column) for f in selected_facets]] : ()
        panel_count=orientations[index]===:rows ?
                    length(first(published).coordinates[index][2]) : length(selected_facets)
        _addon_capacity(layout, matrix, (panel_count,))
    end
    pages=_addon_observation_pages(facets, capacities; automatic = layout===nothing)
    if panel_titles isa Union{Tuple, AbstractVector}
        length(panel_titles)==length(facets) ||
            throw(DimensionMismatch("panel_titles must contain one title per selected panel before pagination"))
        panel_titles=Dict((f.identity, f.panel_identity)=>label
        for (f, label) in zip(facets, panel_titles))
    end
    identities=Set(f.panel_identity for f in facets)
    all(pair -> _addon_panel_identity(first(pair)) in identities, _addon_panel_legend_pairs(panel_legends)) ||
        throw(ArgumentError("panel legend identity is absent from the selected panels"))
    legend_cap=_addon_legend_fraction(legend_cap)
    errorbar_sampling===nothing || errorbar_sampling in (:staggered, :all) ||
        throw(ArgumentError("errorbar_sampling must be :staggered or :all"))
    if figure_title isa Union{Tuple, AbstractVector}
        length(figure_title)==length(pages) || throw(DimensionMismatch(
            "figure_title must contain one entry per generated figure"))
    end
    _addon_activate_backend(backend)
    built=LineCableModels.UIPlot[]
    for (page_index, page) in enumerate(pages)
        request_index=first(page.facets).request_index
        curves=unique(curve -> curve.group,
            [curve for facet in page.facets for curve in facet.curves])
        coordinate_overlay=orientations[request_index]!==:gridpoints
        series_count=coordinate_overlay ? length(coordinate_labels[request_index]) :
                     length(observed)
        series_indices=[curve.slot for curve in curves]
        styles=_addon_comparison_styles(series_indices,
            [curve.role for curve in curves], series_count)
        page_attributes=coordinate_overlay ?
                        [coordinate_attributes[request_index][curve.position]
                         for curve in curves] :
                        [attributes[curve.source_index] for curve in curves]
        position=legend_position===_omitted ?
                 (length(curves)>1 || explicit_labels ? :bottom : nothing) :
                 legend_position
        sampling=errorbar_sampling===nothing ? (length(curves)>1 ? :staggered : :all) :
                 errorbar_sampling
        automatic_title=_quantity_page_title(page)
        page_title=title===nothing ? automatic_title : title
        title!==nothing && length(pages)>1 && (page_title="$title — $automatic_title")
        title===nothing && title_prefix!==nothing && !ismissing(title_prefix) &&
            (page_title="$title_prefix — $automatic_title")
        visible_title=figure_title isa Union{Tuple, AbstractVector} ?
                      figure_title[page_index] : figure_title
        push!(built,
            with_theme(_addon_theme(export_theme = export_theme)) do
                _addon_line_page(
                    published, page;
                    series_indices,
                    errorbar_sampling = sampling, series_defaults = styles,
                    series_attributes = page_attributes, title = page_title,
                    figure_title = visible_title, title_attributes, panel_titles,
                    fig_size = _addon_figure_size(fig_size, page.capacity),
                    xscale, yscale,
                    legend_position = position, legend_title, legend_attributes, legend_cap, panel_legends,
                    controls, display_plot = false, export_theme, open_export, kwargs...)
            end)
        last(built).addon_state=merge(last(built).addon_state,
            (observed = sources,
                display_groups = ((
                    request = requests[request_index], groups = groups[request_index]),),
                displayed_indices = orientations[request_index]===:rows ?
                                    [first(page.facets).source_index] : retained[request_index],
                nominal_capacity = page.capacity))
        _addon_frame_budget(last(built), page.capacity)
    end
    if layout===nothing
        for index in eachindex(requests)
            family=findall(page -> first(page.facets).request_index==index, pages)
            _addon_calibrate_frames!(built[family], capacities[index])
        end
    else
        _addon_calibrate_frames!(built, first(capacities))
    end
    if layout===nothing
        for (p, page) in zip(built, pages)
            (first(page.facets).orientation===:rows ||
             first(page.facets).kind in (:diagonal, :vector) &&
             first(page.facets).orientation===:gridpoints) &&
                _addon_responsive_axis_grid!(p)
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
    # Source names describe gridpoint curves. Row overlays use physical row labels.
    labels=series_labels===nothing && get(kwargs,:overlay,:auto)!==:rows ?
           Tuple(string.(keys(sources))) : series_labels
    return plot(Tuple(values(sources)), selection; ydata, series_labels = labels, kwargs...)
end
function Makie.plot(sources::NamedTuple{K, <:Tuple{Vararg{_PlotSource}}}, args...; kwargs...) where {K}
    plot(sources, args...; kwargs...)
end
