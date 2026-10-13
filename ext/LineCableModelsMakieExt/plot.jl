const _PrimarySource = Union{
    LineCableModels.LineParameters, LineCableModels.SeriesImpedance,
    LineCableModels.ShuntAdmittance, LineCableModels.CableConstants,
    LineCableModels.PropagationParameters}
const _PlotSource = Union{_PrimarySource, Commons.ObservedResult}

# Split only public keyword ownership. Scientific normalization stays in Commons.
function _plot_observation_options(kwargs; retained = false)
    haskey(kwargs, :complete_pairs) &&
        throw(ArgumentError("complete_pairs is owned by observation construction"))
    keys=(:clip, :atol, :units, :length_unit, :quantity_units, :frequency_unit, :frequencies)
    if retained
        for key in (:clip, :atol, :frequencies)
            haskey(kwargs, key) &&
                throw(ArgumentError("$key requires primary observation construction; retained plotting only selects or re-expresses recorded values"))
        end
    end
    acquisition=(; (key=>value for (key, value) in kwargs if key in keys)...)
    presentation=(; (key=>value for (key, value) in kwargs if key ∉ keys)...)
    return acquisition, presentation
end

function plot(source::_PrimarySource, selection = nothing; reference = nothing, kwargs...)
    return plot((source,), selection; reference, kwargs...)
end
Makie.plot(source::_PrimarySource, args...; kwargs...) = plot(source, args...; kwargs...)

function plot(observed::Commons.ObservedResult, selection = nothing; kwargs...)
    return plot([observed], selection; kwargs...)
end

function plot(observed::AbstractVector{<:Commons.ObservedResult}, selection = nothing;
        reference = nothing,
        series_labels = nothing, series_attributes = nothing, problem = nothing,
        formulations = nothing, band = nothing,
        errorbar_sampling = nothing, title = nothing, title_prefix = nothing, figure_title = nothing,
        title_attributes::NamedTuple = (;), panel_titles = nothing, fig_size = nothing, layout = nothing,
        overlay::Symbol = :auto,
        xscale = nothing, yscale = :linear, legend_position = _omitted, legend_title = nothing,
        legend_attributes::NamedTuple = (;), legend_max_fraction = 0.5, panel_legends = nothing,
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
    # A number selects one original identity, as a collection of one. The order of
    # `formulations` orders the points.
    selected=findall(observed) do point
        id=get(point.gridpoint, :id, nothing)
        (problem===nothing || id!==nothing && id.problem_index in problem) &&
            (formulations===nothing || id!==nothing && id.formulation_index in formulations)
    end
    isempty(selected) &&
        throw(ArgumentError("no retained observations match the original point selection"))
    formulations===nothing || sort!(selected;
        by = index -> findfirst(==(observed[index].gridpoint.id.formulation_index), formulations))
    points=observed[selected]
    requests=Commons.observation_selection(first(points), selection)
    requests=Commons.observation_requests(first(points), requests).displayed
    if reference!==nothing && !(reference isa Commons.ObservedResult)
        reference isa _PrimarySource ||
            throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Commons.ObservedResult(reference, requests; complete_pairs = true)
    end
    if !isempty(display_units)
        points=[Commons.ObservedResult(point, requests; display_units...)
                    for point in points]
        reference===nothing ||
            (reference=Commons.ObservedResult(reference, requests; display_units...))
    end
    sources=reference===nothing ? Tuple(points) : (Tuple(points)..., reference)
    input_point_count=length(observed)+(reference===nothing ? 0 : 1)
    displayed=reference===nothing ? selected : [selected; input_point_count]
    slots=reference===nothing ? selected : [selected; 0]
    reference_id=reference===nothing ? nothing : get(reference.gridpoint, :id, nothing)
    products=map(request -> Commons.observation_product(sources, request; band, reference_id), requests)
    plot_data=map(eachindex(sources)) do index
        _line_plot_data(Tuple(records[index] for records in products))
    end
    # The selected cardinality includes an explicit reference and precedes
    # equivalence grouping. The grouping only changes gridpoint curves.
    orientations=map(eachindex(requests)) do index
        kind=first(plot_data).observations[index].coordinates.kind
        overlay===:rows && kind!==:matrix && throw(ArgumentError(
            "overlay=:rows requires matrix observations; select matrix quantities in a separate call"))
        overlay===:auto ?
        (kind===:vector && length(sources)==1 &&
         (layout===nothing || layout==(1, 1)) ? :coordinates : :gridpoints) : overlay
    end
    groups=map(request -> Commons.observation_groups(points; request), requests)
    retained=map(eachindex(requests)) do index
        indices=orientations[index]!==:gridpoints ? collect(eachindex(sources)) :
                [group.representative for group in groups[index]]
        orientations[index]===:gridpoints && reference!==nothing &&
            push!(indices, length(sources))
        indices
    end
    coordinate_identities=map(eachindex(requests)) do index
        first(plot_data).observations[index].coordinates.kind===:assemblies ?
        unique(vcat((source.observations[index].coordinates.labels[
                         source.observations[index].coordinates.assemblies]
        for source in plot_data)...)) : nothing
    end
    coordinate_labels=map(eachindex(requests)) do index
        coordinate_identities[index]===nothing ?
        _coordinate_labels(first(plot_data).observations[index],
            first(plot_data).coordinates[index], orientations[index]) :
        unique(vcat((_coordinate_labels(source.observations[index],
                         source.coordinates[index], orientations[index]) for source in plot_data)...))
    end
    # Positional overrides bind to the overlaid dimension. Preflight every
    # family before constructing the first native figure.
    if series_labels !== nothing || series_attributes !== nothing
        counts=unique(orientations[index]!==:gridpoints ?
                      length(coordinate_labels[index]) :
                      input_point_count for index in eachindex(requests))
        length(counts)==1 || throw(ArgumentError(
            "positional series overrides span figure families with different overlaid dimensions; use separate plot calls"))
    end
    explicit_labels=series_labels!==nothing
    result_labels=explicit_labels && reference!==nothing &&
                     all(==(:gridpoints), orientations) &&
                     length(series_labels)==length(observed)
    labels=explicit_labels && any(==(:gridpoints), orientations) ?
           Tuple(validate(series_labels, Legend, result_labels ? input_point_count-1 :
                                                   input_point_count)[i]
    for i in (result_labels ? selected : displayed)) : nothing
    attributes=any(==(:gridpoints), orientations) ?
               _series_attributes(series_attributes, input_point_count)[displayed] : nothing
    coordinate_attributes=map(eachindex(requests)) do index
        orientations[index]!==:gridpoints ?
        _series_attributes(series_attributes, length(coordinate_labels[index])) : nothing
    end
    request_labels=map(eachindex(requests)) do index
        request=requests[index]
        if orientations[index]!==:gridpoints
            result=Commons.observation_labels(sources; request, fallback = "")
            reference!==nothing &&
                (result[end]*=isempty(result[end]) ? "Reference" : " (reference)")
            result
        else
            result=explicit_labels ? Any[labels...] :
                   Commons.observation_labels(sources; request,
                fallback = length(sources)==1 ? "" : nothing)
            result_labels && push!(result,
                last(Commons.observation_labels(sources; request)))
            (!explicit_labels || result_labels) && reference!==nothing &&
                (result[end]*=" (reference)")
            result
        end
    end
    if explicit_labels
        for index in eachindex(requests)
            orientations[index]===:gridpoints && continue
            coordinate_labels[index]=collect(validate(series_labels, Legend,
                length(coordinate_labels[index])))
        end
    end
    facets=_observation_facets(
        plot_data, requests, orientations, retained, slots, input_point_count,
        request_labels, coordinate_labels, coordinate_identities; explicit_labels)
    capacities=map(eachindex(requests)) do index
        selected_facets=filter(f -> f.request_index==index, facets)
        matrix=first(selected_facets).kind===:matrix && orientations[index]===:gridpoints ?
               [[(f.row, f.column) for f in selected_facets]] : ()
        panel_count=orientations[index]===:rows ?
                    length(first(plot_data).coordinates[index][2]) : length(selected_facets)
        _page_capacity(layout, matrix, (panel_count,))
    end
    pages=_observation_pages(facets, capacities; automatic = layout===nothing)
    panel_titles===nothing || (panel_titles=Dict((f.identity, f.panel_identity)=>label
        for (f, label) in zip(facets, validate(panel_titles, Axis, length(facets)))))
    identities=Set(f.panel_identity for f in facets)
    all(pair -> _panel_identity(first(pair)) in identities, _panel_legend_pairs(panel_legends)) ||
        throw(ArgumentError("panel legend identity is absent from the selected panels"))
    legend_max_fraction=_legend_fraction(legend_max_fraction)
    errorbar_sampling===nothing || errorbar_sampling in (:staggered, :all) ||
        throw(ArgumentError("errorbar_sampling must be :staggered or :all"))
    _activate_plot_backend(backend)
    built=LineCableModels.UIPlot[]
    for page in pages
        request_index=first(page.facets).request_index
        curves=unique(curve -> curve.group,
            [curve for facet in page.facets for curve in facet.curves])
        coordinate_overlay=orientations[request_index]!==:gridpoints
        series_count=coordinate_overlay ? length(coordinate_labels[request_index]) :
                     length(observed)
        series_indices=[curve.slot for curve in curves]
        styles=_comparison_styles(series_indices,
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
        push!(built,
            with_theme(_plot_theme(export_theme = export_theme)) do
                _line_page(
                    plot_data, page;
                    series_indices,
                    errorbar_sampling = sampling, series_defaults = styles,
                    series_attributes = page_attributes, title = page_title,
                    figure_title, title_attributes, panel_titles,
                    fig_size = _figure_size(fig_size, page.capacity),
                    xscale, yscale,
                    legend_position = position, legend_title, legend_attributes, legend_max_fraction, panel_legends,
                    controls, display_plot = false, export_theme, open_export, kwargs...)
            end)
        last(built).plot_state=merge(last(built).plot_state,
            (observed = sources,
                display_groups = ((
                    request = requests[request_index], groups = groups[request_index]),),
                displayed_indices = orientations[request_index]===:rows ?
                                    [first(page.facets).source_index] : retained[request_index],
                nominal_capacity = page.capacity))
        _frame_budget(last(built), page.capacity)
    end
    if layout===nothing
        for index in eachindex(requests)
            family=findall(page -> first(page.facets).request_index==index, pages)
            _calibrate_frames!(built[family], capacities[index])
        end
    else
        _calibrate_frames!(built, first(capacities))
    end
    if layout===nothing
        for (p, page) in zip(built, pages)
            (first(page.facets).orientation===:rows ||
             first(page.facets).kind in (:diagonal, :vector) &&
             first(page.facets).orientation===:gridpoints) &&
                _responsive_axis_grid!(p)
        end
    end
    if display_plot
        foreach(page -> _display_figure!(page.figure, page.export_name), built)
    end
    return length(built)==1 ? only(built) : built
end
function Makie.plot(source::Commons.ObservedResult, args...; kwargs...)
    plot(source, args...; kwargs...)
end
function Makie.plot(source::AbstractVector{<:Commons.ObservedResult}, args...; kwargs...)
    plot(source, args...; kwargs...)
end

function plot(sources::Union{AbstractVector{<:_PlotSource}, Tuple{Vararg{_PlotSource}}},
        selection = nothing;
        reference = nothing, kwargs...)
    isempty(sources) && throw(ArgumentError("plot requires at least one result"))
    if all(source -> source isa Commons.ObservedResult, sources)
        return plot(Commons.ObservedResult[sources...], selection; reference, kwargs...)
    end
    acquisition, presentation=_plot_observation_options(kwargs)
    requests=Commons.observation_selection(first(sources), selection)
    normalized=Commons.observation_requests(first(sources), requests; complete_pairs = true)
    retained_units=(;
        (key=>value
    for (key, value) in pairs(acquisition)
    if key in
       (:units, :length_unit, :quantity_units, :frequency_unit))...)
    observed=Commons.ObservedResult[source isa Commons.ObservedResult ?
                                    (isempty(retained_units) ? source :
                                     Commons.ObservedResult(source, normalized.displayed; retained_units...)) :
                                    Commons.ObservedResult(source, requests;
                                        complete_pairs = true, acquisition...)
                                    for source in sources]
    if reference!==nothing && !(reference isa Commons.ObservedResult)
        reference isa _PrimarySource ||
            throw(ArgumentError("construct an atomic ObservedResult for a reference collection"))
        reference=Commons.ObservedResult(reference, requests; complete_pairs = true, acquisition...)
    end
    return plot(observed, normalized.displayed; reference, presentation...)
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

function plot(first::_PrimarySource, second::_PrimarySource, rest...; kwargs...)
    sources=_PrimarySource[first, second]
    selection=nothing
    for item in rest
        item isa _PrimarySource && selection===nothing ? push!(sources, item) :
        selection===nothing ? (selection=item) :
        throw(ArgumentError("one trailing selection is accepted"))
    end
    return plot(sources, selection; kwargs...)
end
function Makie.plot(first::_PrimarySource, second::_PrimarySource, rest...; kwargs...)
    plot(first, second, rest...; kwargs...)
end

function plot(sources::NamedTuple, selection = nothing; series_labels = nothing, kwargs...)
    all(source -> source isa _PlotSource, values(sources)) ||
        throw(ArgumentError("named result plotting requires supported primary or observed members"))
    # Source names describe gridpoint curves. Row overlays use physical row labels.
    labels=series_labels===nothing && get(kwargs,:overlay,:auto)!==:rows ?
           Tuple(string.(keys(sources))) : series_labels
    return plot(Tuple(values(sources)), selection; series_labels = labels, kwargs...)
end
function Makie.plot(sources::NamedTuple{K, <:Tuple{Vararg{_PlotSource}}}, args...; kwargs...) where {K}
    plot(sources, args...; kwargs...)
end
