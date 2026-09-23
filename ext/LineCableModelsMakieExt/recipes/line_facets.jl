# Semantic line faceting -----------------------------------------------------
#
# Matrix coordinates identify axes. Result containers identify series. These
# helpers deliberately normalize only the data needed to construct native Makie
# blocks; they are not a second plot-specification model.

function _semantic_line_facets(published, ydata)
    facets = NamedTuple[]
    for request_index in eachindex(ydata)
        product = first(published).observations[request_index]
        c = product.coordinates
        rows, columns, _ = first(published).coordinates[request_index]
        extent=c.kind===:matrix ?
               ntuple(
            d -> maximum(source.observations[request_index].coordinates.extent[d]
            for source in published),
            length(c.extent)) : c.extent
        for (local_row, row) in enumerate(rows),
            (local_column, column) in enumerate(columns)

            push!(facets,
                (; request_index, local_row, local_column, row,
                    column = c.kind===:diagonal ? row : column,
                    kind = c.kind, extent, domain = get(c, :domain, :unspecified),
                    quantity = product.quantity, identity = request_identity(ydata[request_index])))
        end
    end
    return facets
end

function _semantic_line_pages(facets, capacity; automatic = false)
    pages=NamedTuple[]
    for request_index in unique(facet.request_index for facet in facets)
        selected=filter(facet -> facet.request_index==request_index, facets)
        first_facet=first(selected)
        if first_facet.kind!==:matrix
            append!(pages, _addon_flow_pages(selected, capacity))
            continue
        end
        positions=[(f.row, f.column) for f in selected]
        origin=automatic ? _addon_panel_footprint(positions).origin : (1, 1)
        for page in _addon_matrix_pages(positions, first_facet.extent, capacity; origin)
            push!(pages,
                (; facets = selected[page.members], positions = page.positions,
                    dimensions = page.dimensions, origin = page.origin, index = page.index))
        end
    end
    return pages
end

function _semantic_quantity_title(object, facet)
    label=Units.label(facet.quantity)
    facet.identity isa Tuple && first(facet.identity)===LineCableModels.statistics &&
        (label *= " · " * string(last(facet.identity)))
    facet.kind in (:matrix, :diagonal) || return label
    is_modal = facet.domain === :ModalDomain
    coordinate = is_modal ? "mode" : "conductor"

    if facet.row == facet.column
        return "$label, $coordinate $(facet.row)"
    elseif is_modal
        return "$label, intermodal residual ($(facet.row), $(facet.column))"
    else
        return "$label, conductors $(facet.row) → $(facet.column)"
    end
end

function _semantic_page_title(object, page)
    facet=first(page.facets)
    label=Units.label(facet.quantity)
    facet.identity isa Tuple && first(facet.identity)===LineCableModels.statistics &&
        (label *= " · " * string(last(facet.identity)))
    return page.index==(1, 1) ? label : "$label ($(page.index[1]),$(page.index[2]))"
end

function _semantic_page_option(value, page_index::Int, page_count::Int, name::AbstractString)
    value === nothing && return nothing
    if value isa Tuple || value isa AbstractVector
        length(value) == page_count || throw(DimensionMismatch(
            "$name must contain one entry per generated figure",
        ))
        return value[page_index]
    end
    return value
end

function _semantic_panel_title(
        panel_titles, object, facet, panel_index::Int, panel_count::Int)
    panel_titles === nothing && return _semantic_quantity_title(object, facet)
    panel_titles isa Function && return panel_titles(facet)
    if panel_titles isa AbstractDict
        quantity_symbol = Symbol(Units.symbol(facet.quantity))
        candidates = (
            (facet.identity, facet.row, facet.column),
            (quantity_symbol, facet.row, facet.column),
            (facet.row, facet.column),
            facet.identity,
            quantity_symbol
        )
        for candidate in candidates
            haskey(panel_titles, candidate) && return panel_titles[candidate]
        end
        return _semantic_quantity_title(object, facet)
    end
    panel_titles isa Tuple || panel_titles isa AbstractVector ||
        throw(ArgumentError(
            "panel_titles must be a tuple, vector, dictionary, function, or nothing",
        ))
    length(panel_titles) == panel_count || throw(DimensionMismatch(
        "panel_titles must contain one entry per subplot on each generated figure",
    ))
    return panel_titles[panel_index]
end

function _addon_semantic_line_page(
        object,
        published,
        source_labels,
        page,
        ;
        series_indices,
        series_count,
        errorbar_sampling,
        series_defaults,
        series_attributes,
        title,
        figure_title,
        title_attributes,
        panel_titles,
        fig_size,
        xscale,
        yscale,
        legend_position,
        legend_title,
        legend_attributes,
        legend_overflow,
        panel_legends,
        controls,
        display_plot,
        export_theme,
        open_export,
        kwargs...
)
    shell = _addon_shell(; size = fig_size, controls, kwargs...)
    shell.canvas.default_rowgap = Fixed(24)
    shell.canvas.default_colgap = Fixed(48)
    rowgap!(shell.canvas, 24)
    colgap!(shell.canvas, 48)
    # The occupied rectangle preserves internal holes, not unselected perimeter.
    cells=Dict{Tuple{Int, Int}, Any}()
    for row in 1:page.dimensions[1], column in 1:page.dimensions[2]

        cell=_addon_panel!(shell, (row, column))
        rowsize!(shell.canvas, row, Auto(false, 1))
        colsize!(shell.canvas, column, Auto(false, 1))
        cells[(row, column)]=cell
    end
    axes = Any[]
    axis_series = Vector{NamedTuple}[]
    panels = Any[]
    resets = Function[]
    requested_scales = NamedTuple[]
    groups = Dict{Symbol, Vector{Any}}()
    dependent_plots = Pair{Makie.Plot, Makie.Plot}[]
    marker_coordinates = Dict{Makie.Plot, Any}()
    group_order = Symbol[]
    group_labels = Dict{Symbol, Any}()
    panel_group_labels = Any[]
    colors = Tuple(series_defaults === nothing ? _addon_comparison_color(index) :
                   series_defaults[i].attributes.color
    for (i, index) in enumerate(series_indices))

    for (panel_index, (facet, position)) in enumerate(zip(page.facets, page.positions))
        observation = first(published).observations[facet.request_index]
        xvalues = collect(Iterators.flatten(source.frequencies[facet.request_index].values
        for source in published))
        yvalues = collect(Iterators.flatten(
            view(source.observations[facet.request_index].values,
                facet.local_row, facet.local_column, :) for source in published
        ))
        xobservation = merge(first(published).frequencies[facet.request_index], (;
            values = xvalues))
        yobservation = merge(observation, (; values = yvalues))
        panel=merge(cells[position], (; logical_position = (facet.row, facet.column)))
        row, column = position
        bottom_row = maximum(first, page.positions)
        attributes = (;
            xlabelvisible = row == bottom_row,
            xticklabelsvisible = row == bottom_row,
            xticksvisible = row == bottom_row,
            # Matrix cells have independent y limits; each must expose its scale.
            ylabelvisible = true,
            yticklabelsvisible = true,
            yticksvisible = true
        )
        if all(
            source -> source.resolutions[facet.request_index].clip &&
                      source.resolutions[facet.request_index].kind === :declared_floor,
            published)
            if all(ismissing, yvalues)
                attributes = merge(attributes,
                    (subtitle = facet.identity isa Tuple && angle in facet.identity ?
                                "Undefined phase" : "Unavailable quantity",))
            end
        end
        axis, scales = _addon_axis!(
            panel.content,
            xobservation,
            yobservation;
            title = _semantic_panel_title(
                panel_titles, object, facet, panel_index, length(page.facets)),
            xscale = xscale===nothing ?
                     (facet.kind in (:matrix, :diagonal) ? :log10 : :linear) : xscale,
            yscale,
            xlabel = get(xobservation, :label, nothing),
            attributes = facet.kind===:assemblies ?
                         merge(attributes, (dim1_conversion = CategoricalConversion(),)) :
                         attributes,
            native_attributes = shell.axis_attributes
        )
        series = NamedTuple[]
        scoped_labels = Dict{Symbol, Any}()
        for source_index in eachindex(published)
            source = published[source_index]
            curve = collect(view(
                source.observations[facet.request_index].values,
                facet.local_row,
                facet.local_column,
                :
            ))
            errors = get(source.observations[facet.request_index], :errors, nothing)
            yerror = errors === nothing ? nothing :
                     collect(view(errors, facet.local_row, facet.local_column, :))
            group = Symbol("result_$source_index")
            source_label = source_labels[source_index]
            interval_support=Dict{Makie.Plot, Any}()
            draw! = facet.kind in (:assemblies, :array) ? _addon_points! : _addon_line!
            plots = draw!(
                axis,
                source.frequencies[facet.request_index].values,
                curve;
                dependent_plots,
                label = source_label,
                color = colors[source_index],
                phase = series_defaults === nothing ?
                        (series_indices[source_index], series_count) :
                        series_defaults[source_index].phase,
                endpoints = series_defaults !== nothing &&
                            series_defaults[source_index].endpoints,
                marker_coordinates = series_defaults === nothing ? nothing :
                                     marker_coordinates,
                errorbar_sampling, yerror, interval_support
            )
            if !haskey(groups, group)
                groups[group] = Any[]
                push!(group_order, group)
                group_labels[group] = source_label
            end
            append!(groups[group], plots)
            scoped_labels[group] = source_label
            binding=(;
                xdata = facet.kind===:assemblies ? nothing :
                        source.frequencies[facet.request_index].values,
                ydata = curve,
                yerror,
                interval_support,
                sampled_intervals = facet.kind in (:matrix, :diagonal) &&
                                    errorbar_sampling === :staggered,
                plots
            )
            support(dimension) = begin
                values=_addon_visible_values((binding,), dimension; include_uncertainty = true)
                isempty(values) ? () : extrema(values)
            end
            push!(series, merge(binding, (full_support = (
                x = support(:x), y = support(:y)),)))
        end
        reset! = _addon_reset!(axis, series)
        push!(axes, axis)
        push!(axis_series, series)
        push!(panels, panel)
        push!(panel_group_labels, scoped_labels)
        push!(resets, reset!)
        push!(requested_scales, scales)
    end
    selected=Set((facet.row, facet.column) for facet in page.facets)
    local_panel_legends=Tuple(pair
    for pair in _addon_panel_legend_pairs(panel_legends) if first(pair) in selected)
    built = _addon_finish!(
        shell, axes, resets, groups, group_order, group_labels;
        requested_scales, axis_series,
        dependent_plots,
        series_attributes,
        series_defaults,
        marker_coordinates,
        title,
        figure_title,
        title_attributes,
        legend_position,
        legend_attributes,
        legend_overflow,
        legend_title,
        panels, frame_cells = Tuple(values(cells)),
        panel_legends = local_panel_legends,
        panel_group_labels,
        controls,
        display_plot,
        export_name = title,
        export_theme,
        open_export
    )
    built.addon_state=merge(built.addon_state,
        (
            panel_page = (;
                index = page.index, dimensions = page.dimensions, origin = get(page, :origin, nothing),
                coordinates = Tuple((f.row, f.column) for f in page.facets)),
            page_cells = cells))
    return built
end
