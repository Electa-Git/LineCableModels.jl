# Semantic line faceting -----------------------------------------------------
#
# Matrix coordinates identify axes. Result containers identify series. These
# helpers deliberately normalize only the data needed to construct native Makie
# blocks; they are not a second plot-specification model.

function _semantic_line_facets(published, requests)
    facets = NamedTuple[]
    for request_index in eachindex(requests)
        observation = published.observations[request_index]
        rows, columns, _ = published.coordinates[request_index]
        diagonal = _diagonal_request(requests[request_index])
        for (local_row, source_row) in enumerate(rows),
            (local_column, source_column) in enumerate(columns)
            push!(facets, (;
                request_index,
                local_row,
                local_column,
                row = source_row,
                column = diagonal ? source_row : source_column,
                diagonal,
                family = _line_request_family(requests[request_index]),
                quantity = observation.quantity,
                identity = request_identity(requests[request_index])
            ))
        end
    end
    return facets
end

function _semantic_line_layout_mode(object, facets, layout)
    layout === nothing || begin
        layout isa Tuple && length(layout) == 2 &&
        all(value -> value isa Integer && !(value isa Bool) && value > 0, layout) ||
            throw(ArgumentError(
                "layout must be a tuple of two positive integers or nothing",
            ))
        dimensions = Tuple(Int.(layout))
        dimensions == (1, 1) && return :individual
        dimensions in ((1, 2), (2, 1)) && return :paired
        matrix_size = size(object isa LineParameters ? Z(object) : object, 1)
        dimensions == (matrix_size, matrix_size) && return :matrix
        throw(ArgumentError(
            "line layout must be (1, 1), (1, 2), (2, 1), " *
            "($matrix_size, $matrix_size), or nothing",
        ))
    end
    coordinates = unique((facet.row, facet.column) for facet in facets)
    families = unique(facet.family for facet in facets)
    return length(coordinates) == 1 && length(families) == 1 &&
           length(facets) > 1 ? :paired : :matrix
end

function _semantic_line_pages(object, facets, layout)
    mode = _semantic_line_layout_mode(object, facets, layout)
    if mode === :individual
        return mode, [(; facets = Any[facet], positions = ((1, 1),),
            dimensions = (1, 1)) for facet in facets]
    end

    keys = Any[]
    grouped = Vector{Vector{Any}}()
    key = mode === :paired ?
          (facet -> (facet.family, facet.row, facet.column)) :
          (facet -> facet.quantity)
    for facet in facets
        facet_key = key(facet)
        index = findfirst(==(facet_key), keys)
        if index === nothing
            push!(keys, facet_key)
            push!(grouped, Any[facet])
        else
            push!(grouped[index], facet)
        end
    end

    matrix_size = size(object isa LineParameters ? Z(object) : object, 1)
    pages = NamedTuple[]
    for page_facets in grouped
        positions, dimensions = if mode === :paired
            requested = layout === nothing ? (1, length(page_facets)) : layout
            _addon_positions(length(page_facets), requested)
        else
            matrix_positions = Tuple((facet.row, facet.column) for facet in page_facets)
            length(unique(matrix_positions)) == length(matrix_positions) || throw(
                ArgumentError(
                    "the observation selection maps more than one curve to the same physical subplot",
                ),
            )
            matrix_positions, (matrix_size, matrix_size)
        end
        push!(pages, (; facets = page_facets, positions, dimensions))
    end
    return mode, pages
end

_semantic_coordinate_name(::LineCableModels.LineParameters{T, U, D}) where {
    T, U, D <: LineCableModels.ModalDomain} = "mode"
_semantic_coordinate_name(_) = "conductor"

function _semantic_quantity_title(object, facet)
    quantity_label = LineCableModels.Units.label(facet.quantity)
    description = lowercasefirst(quantity_label)
    coordinate = _semantic_coordinate_name(object)
    if coordinate == "mode"
        return facet.row == facet.column ?
               "$quantity_label — mode $(facet.row)" :
               "$quantity_label — mode $(facet.row) → mode $(facet.column)"
    elseif facet.row == facet.column
        return "Self $description — conductor $(facet.row)"
    end
    return "Mutual $description — conductor $(facet.row) → conductor $(facet.column)"
end

function _semantic_relation_title(object, quantity, row, column)
    quantity_label = LineCableModels.Units.label(quantity)
    description = lowercasefirst(quantity_label)
    coordinate = _semantic_coordinate_name(object)
    if coordinate == "mode"
        return row == column ?
               "$quantity_label — mode $row" :
               "$quantity_label — mode $row → mode $column"
    elseif row == column
        return "Self $description — conductor $row"
    end
    return "Mutual $description — conductor $row → conductor $column"
end

function _semantic_page_title(object, page, mode)
    first_facet = first(page.facets)
    mode === :individual && return _semantic_quantity_title(object, first_facet)
    mode === :paired && return _semantic_relation_title(
        object,
        LineCableModels.Units.quantity(_family_parent(first_facet.family)),
        first_facet.row,
        first_facet.column
    )
    return LineCableModels.Units.label(first_facet.quantity)
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

function _semantic_panel_title(panel_titles, object, facet, panel_index::Int, panel_count::Int)
    panel_titles === nothing && return _semantic_quantity_title(object, facet)
    panel_titles isa Function && return String(panel_titles(facet))
    if panel_titles isa AbstractDict
        quantity_symbol = Symbol(LineCableModels.Units.symbol(facet.quantity))
        candidates = (
            (facet.identity, facet.row, facet.column),
            (quantity_symbol, facet.row, facet.column),
            (facet.row, facet.column),
            facet.identity,
            quantity_symbol
        )
        for candidate in candidates
            haskey(panel_titles, candidate) && return String(panel_titles[candidate])
        end
        return _semantic_quantity_title(object, facet)
    end
    panel_titles isa Tuple || panel_titles isa AbstractVector || throw(ArgumentError(
        "panel_titles must be a tuple, vector, dictionary, function, or nothing",
    ))
    length(panel_titles) == panel_count || throw(DimensionMismatch(
        "panel_titles must contain one entry per subplot on each generated figure",
    ))
    return String(panel_titles[panel_index])
end

function _semantic_figure_size(fig_size, dimensions)
    if fig_size !== nothing
        fig_size isa Tuple{Int, Int} && all(>(0), fig_size) || throw(ArgumentError(
            "fig_size must be a tuple of two positive integers or nothing",
        ))
        return fig_size
    end
    rows, columns = dimensions
    return (max(680, 390columns + 180), max(440, 290rows + 100))
end

function _addon_semantic_line_page(
        object,
        published,
        source_labels,
        page,
        mode;
        title,
        figure_title,
        title_attributes,
        panel_titles,
        fig_size,
        xscale,
        yscale,
        legend_position,
        legend_anchor,
        legend_title,
        legend_attributes,
        legend_overflow,
        panel_legends,
        controls,
        display_plot,
        export_theme,
        open_export
)
    shell = _addon_shell(; size = _semantic_figure_size(fig_size, page.dimensions), controls)
    shell.canvas.default_rowgap = Fixed(24)
    shell.canvas.default_colgap = Fixed(48)
    rowgap!(shell.canvas, 24)
    colgap!(shell.canvas, 48)
    axes = Any[]
    panels = Any[]
    resets = Function[]
    xsetters = Function[]
    ysetters = Function[]
    groups = Dict{Symbol, Vector{Any}}()
    group_order = Symbol[]
    group_labels = Dict{Symbol, String}()
    panel_group_labels = Any[]
    colors = Tuple(_addon_comparison_color(index) for index in eachindex(published))

    for (panel_index, (facet, position)) in enumerate(zip(page.facets, page.positions))
        observation = first(published).observations[facet.request_index]
        xvalues = collect(Iterators.flatten(source.frequency.values for source in published))
        yvalues = collect(Iterators.flatten(
            view(source.observations[facet.request_index].values,
                facet.local_row, facet.local_column, :) for source in published
        ))
        xobservation = merge(first(published).frequency, (; values = xvalues))
        yobservation = merge(observation, (; values = yvalues))
        xscales = _axis_scales(xvalues)
        yscales = _axis_scales(yvalues)
        xscale in xscales || throw(DomainError(
            xvalues,
            "logarithmic frequency axes require positive finite data and uncertainty bounds",
        ))
        yscale in yscales || throw(DomainError(
            yvalues,
            "logarithmic ordinate axes require positive finite data and uncertainty bounds",
        ))
        panel = _addon_panel!(shell, position)
        row, column = position
        attributes = (;
            xlabelvisible = row == page.dimensions[1],
            xticklabelsvisible = row == page.dimensions[1],
            xticksvisible = row == page.dimensions[1],
            ylabelvisible = mode !== :matrix || column == 1,
            yticklabelsvisible = mode !== :matrix || column == 1,
            yticksvisible = mode !== :matrix || column == 1
        )
        axis, axis_labels = _addon_axis!(
            panel.content,
            xobservation,
            yobservation;
            title = _semantic_panel_title(
                panel_titles, object, facet, panel_index, length(page.facets)),
            xscale,
            yscale,
            xscales,
            yscales,
            attributes
        )
        series = NamedTuple[]
        scoped_labels = Dict{Symbol, String}()
        for source_index in eachindex(published)
            source = published[source_index]
            curve = collect(view(
                source.observations[facet.request_index].values,
                facet.local_row,
                facet.local_column,
                :
            ))
            group = Symbol("result_$source_index")
            source_label = source_labels[source_index]
            plots = _addon_line!(
                axis,
                source.frequency.values,
                curve;
                label = source_label,
                color = colors[source_index]
            )
            if !haskey(groups, group)
                groups[group] = Any[]
                push!(group_order, group)
                group_labels[group] = source_label
            end
            append!(groups[group], plots)
            scoped_labels[group] = source_label
            push!(series, (;
                xdata = source.frequency.values,
                ydata = curve,
                plots
            ))
        end
        reset! = () -> _addon_reset!(axis, series, axis_labels)
        xsetter = scale -> _addon_set_axis!(
            axis, :x, axis_labels.x, xscales,
            something(_addon_scientific_exponent(xvalues), 0), scale)
        ysetter = scale -> _addon_set_axis!(
            axis, :y, axis_labels.y, yscales,
            something(_addon_scientific_exponent(yvalues), 0), scale)
        push!(axes, axis)
        push!(panels, panel)
        push!(panel_group_labels, scoped_labels)
        push!(resets, reset!)
        :log10 in xscales && push!(xsetters, xsetter)
        :log10 in yscales && push!(ysetters, ysetter)
        on(shell.figure.scene, axis.finallimits) do _
            _addon_refresh_format!(axis, series, axis_labels)
            return nothing
        end
        reset!()
    end
    length(xsetters) == length(axes) || empty!(xsetters)
    length(ysetters) == length(axes) || empty!(ysetters)
    mode === :paired && length(axes) > 1 && _addon_responsive_axis_grid!(
        shell.figure, shell.canvas, panels, axes, page.dimensions)
    return _addon_finish!(
        shell, axes, resets, xsetters, ysetters, groups, group_order, group_labels;
        title,
        figure_title,
        title_attributes,
        legend_position,
        legend_anchor,
        legend_attributes,
        legend_overflow,
        legend_title,
        panels,
        panel_legends,
        panel_group_labels,
        controls,
        display_plot,
        export_name = title,
        export_theme,
        open_export
    )
end

function _addon_line_pages(
        sources::Tuple;
        frequencies = nothing,
        requests,
        series_labels = nothing,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        panel_titles = nothing,
        labels = nothing,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        clip::Bool = true,
        fig_size = nothing,
        layout = nothing,
        xscale::Symbol = :linear,
        yscale::Symbol = :linear,
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_attributes::NamedTuple = (;),
        legend_overflow::Symbol = :ellipsis,
        panel_legends = (),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    _addon_activate_backend(backend)
    legend_overflow in (:ellipsis, :show_all) || throw(ArgumentError(
        "legend_overflow must be :ellipsis or :show_all",
    ))
    legend_attributes isa NamedTuple || throw(ArgumentError(
        "legend_attributes must be a NamedTuple",
    ))
    panel_titles === nothing || labels === nothing || throw(ArgumentError(
        "use panel_titles; labels is only a compatibility alias",
    ))
    resolved_panel_titles = panel_titles === nothing ? labels : panel_titles
    explicit_source_labels = series_labels !== nothing
    source_labels = explicit_source_labels ?
                    _comparison_labels(series_labels, length(sources)) :
                    Tuple("Result $index" for index in eachindex(sources))
    published = if length(sources) == 1
        (_prepare_line_observations(
            only(sources);
            frequencies,
            requests,
            freq_unit,
            length_unit,
            quantity_units,
            clip
        ),)
    else
        _prepare_line_comparison(
            sources;
            requests,
            series_labels = source_labels,
            freq_unit,
            length_unit,
            quantity_units,
            clip
        ).published
    end
    any(source -> length(source.frequency.values) <= 1, published) &&
        return LineCableModels.UIPlot[]

    facets = _semantic_line_facets(first(published), requests)
    mode, pages = _semantic_line_pages(first(sources), facets, layout)
    effective_legend_position = length(sources) > 1 || explicit_source_labels ?
                                legend_position : nothing
    built = LineCableModels.UIPlot[]
    for (page_index, page) in enumerate(pages)
        automatic_title = _semantic_page_title(first(sources), page, mode)
        page_title = title === nothing ? automatic_title : String(title)
        length(pages) > 1 && title !== nothing &&
            (page_title = "$page_title — $automatic_title")
        visible_title = _semantic_page_option(
            figure_title, page_index, length(pages), "figure_title")
        push!(built, with_theme(_addon_theme(export_theme = export_theme)) do
            _addon_semantic_line_page(
                first(sources), published, source_labels, page, mode;
                title = page_title,
                figure_title = visible_title,
                title_attributes,
                panel_titles = resolved_panel_titles,
                fig_size,
                xscale,
                yscale,
                legend_position = effective_legend_position,
                legend_anchor,
                legend_title,
                legend_attributes,
                legend_overflow,
                panel_legends,
                controls,
                display_plot,
                export_theme,
                open_export
            )
        end)
    end
    return length(built) == 1 ? only(built) : built
end

function _addon_semantic_line_plots(
        object;
        frequencies = nothing,
        requests,
        series_labels = nothing,
        legend = nothing,
        legend_labels = nothing,
        kwargs...
)
    supplied_labels = count(!isnothing, (series_labels, legend_labels, legend))
    supplied_labels <= 1 || throw(ArgumentError(
        "use series_labels; legend_labels and legend are compatibility aliases",
    ))
    resolved_labels = series_labels !== nothing ? series_labels :
                      legend_labels !== nothing ? legend_labels : legend
    return _addon_line_pages(
        (object,);
        frequencies,
        requests,
        series_labels = resolved_labels,
        kwargs...
    )
end

function _addon_comparison_color(index::Int)
    hue = mod(210.0 + 137.50776405003785 * (index - 1), 360.0)
    return RGB(HSV(hue, 0.72, 0.78))
end
