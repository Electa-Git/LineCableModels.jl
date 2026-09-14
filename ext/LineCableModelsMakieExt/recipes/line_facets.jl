# Semantic line faceting -----------------------------------------------------
#
# Matrix coordinates identify axes. Result containers identify series. These
# helpers deliberately normalize only the data needed to construct native Makie
# blocks; they are not a second plot-specification model.

function _semantic_line_facets(published, ydata)
    facets = NamedTuple[]
    for request_index in eachindex(ydata)
        observation = published.observations[request_index]
        rows, columns, _ = published.coordinates[request_index]
        diagonal = _diagonal_request(ydata[request_index])
        for (local_row, source_row) in enumerate(rows),
            (local_column, source_column) in enumerate(columns)

            push!(facets,
                (;
                    request_index,
                    local_row,
                    local_column,
                    row = source_row,
                    column = diagonal ? source_row : source_column,
                    diagonal,
                    family = _line_request_family(ydata[request_index]),
                    quantity = observation.quantity,
                    identity = request_identity(ydata[request_index])
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
        matrix_size = size(object isa LineParameters ? Z(object) :
            object isa ObservationPublication ? first(object).values : object, 1)
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

function _semantic_line_pages(object, facets, layout, blocks=nothing)
    if blocks !== nothing
        blocks isa Tuple && length(blocks) == 2 &&
            all(value -> value isa Integer && !(value isa Bool) && value > 0, blocks) ||
            throw(ArgumentError("blocks must be a tuple of two positive integers or nothing"))
        matrix_size = size(object isa LineParameters ? Z(object) :
            object isa ObservationPublication ? first(object).values : object, 1)
        layout === nothing || layout == (matrix_size, matrix_size) ||
            throw(ArgumentError("blocks requires matrix layout; omit layout or use the full matrix dimensions"))
    end
    mode = blocks === nothing ? _semantic_line_layout_mode(object, facets, layout) : :matrix
    if mode === :individual
        return mode,
        [(; facets = Any[facet], positions = ((1, 1),),
             dimensions = (1, 1)) for facet in facets]
    end

    keys = Any[]
    grouped = Vector{Vector{Any}}()
    key = mode === :paired ?
          (facet -> (facet.family, facet.row, facet.column)) :
          (facet -> (facet.quantity,facet.identity))
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

    matrix_size = size(object isa LineParameters ? Z(object) :
        object isa ObservationPublication ? first(object).values : object, 1)
    pages = NamedTuple[]
    for page_facets in grouped
        positions,
        dimensions = if mode === :paired
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
        if blocks === nothing
            push!(pages, (; facets = page_facets, positions, dimensions))
        else
            block_rows, block_columns = Int.(blocks)
            for block_row in 1:cld(matrix_size, block_rows),
                    block_column in 1:cld(matrix_size, block_columns)
                retained = filter(page_facets) do facet
                    cld(facet.row, block_rows) == block_row &&
                        cld(facet.column, block_columns) == block_column
                end
                isempty(retained) && continue
                local_positions = Tuple((mod1(facet.row, block_rows),
                    mod1(facet.column, block_columns)) for facet in retained)
                push!(pages, (; facets=retained, positions=local_positions,
                    dimensions=(block_rows, block_columns), block=(block_row, block_column)))
            end
        end
    end
    return mode, pages
end

function _semantic_coordinate_name(::LineCableModels.LineParameters{
        T, U, D}) where {
        T, U, D <: LineCableModels.ModalDomain}
    "mode"
end
_semantic_coordinate_name(_) = "conductor"

function _semantic_quantity_title(object, facet)
    quantity_label = LineCableModels.Units.label(facet.quantity)
    facet.identity isa Tuple && first(facet.identity) === LineCableModels.statistics &&
        (quantity_label *= " · " * string(last(facet.identity)))
    description = lowercasefirst(quantity_label)
    coordinate = _semantic_coordinate_name(object)
    if coordinate == "mode"
        return facet.row == facet.column ?
               "$quantity_label, mode $(facet.row)" :
               "$quantity_label, modes $(facet.row) → $(facet.column)"
    elseif facet.row == facet.column
        return "Self- $description, conductor $(facet.row)"
    end
    return "Mutual- $description, conductors $(facet.row) → $(facet.column)"
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
        return "Self- $description, conductors $row"
    end
    return "Mutual- $description, conductors $row → $column"
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
    label = LineCableModels.Units.label(first_facet.quantity)
    first_facet.identity isa Tuple && first(first_facet.identity) === LineCableModels.statistics &&
        (label *= " · " * string(last(first_facet.identity)))
    return haskey(page, :block) ? "$label ($(page.block[1]),$(page.block[2]))" : label
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
    panel_titles isa Tuple || panel_titles isa AbstractVector ||
        throw(ArgumentError(
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
        series_indices,
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
        legend_anchor,
        legend_title,
        legend_attributes,
        legend_overflow,
        panel_legends,
        signed_ylog,
        controls,
        display_plot,
        export_theme,
        open_export,
        kwargs...
)
    shell = _addon_shell(; size = _semantic_figure_size(fig_size, page.dimensions), controls, kwargs...)
    shell.canvas.default_rowgap = Fixed(24)
    shell.canvas.default_colgap = Fixed(48)
    rowgap!(shell.canvas, 24)
    colgap!(shell.canvas, 48)
    blocked = haskey(page, :block)
    # Reserve the complete footprint, including residual/empty cells. These are
    # layout tracks, not dummy axes or numerical observations.
    cells = Dict{Tuple{Int,Int},Any}()
    if blocked
        for row in 1:page.dimensions[1], column in 1:page.dimensions[2]
            cell = _addon_panel!(shell, (row, column))
            rowsize!(shell.canvas, row, Auto(false, 1))
            colsize!(shell.canvas, column, Auto(false, 1))
            cell.layout.alignmode = Outside()
            cells[(row, column)] = cell
        end
    end
    axes = Any[]
    panels = Any[]
    resets = Function[]
    requested_scales = NamedTuple[]
    groups = Dict{Symbol, Vector{Any}}()
    dependent_plots = Pair{Makie.Plot,Makie.Plot}[]
    group_order = Symbol[]
    group_labels = Dict{Symbol, String}()
    panel_group_labels = Any[]
    colors = Tuple(_addon_comparison_color(index) for index in series_indices)

    for (panel_index, (facet, position)) in enumerate(zip(page.facets, page.positions))
        observation = first(published).observations[facet.request_index]
        xvalues = collect(Iterators.flatten(source.frequency.values
        for source in published))
        yvalues = collect(Iterators.flatten(
            view(source.observations[facet.request_index].values,
                facet.local_row, facet.local_column, :) for source in published
        ))
        xobservation = merge(first(published).frequency, (; values = xvalues))
        yobservation = merge(observation, (; values = yvalues))
        panel = if blocked
            merge(cells[position], (; logical_position=(facet.row, facet.column)))
        else
            _addon_panel!(shell, position)
        end
        row, column = position
        bottom_row = blocked ? maximum(first, page.positions) : page.dimensions[1]
        attributes = (;
            xlabelvisible = row == bottom_row,
            xticklabelsvisible = row == bottom_row,
            xticksvisible = row == bottom_row,
            # Matrix cells have independent y limits; each must expose its scale.
            ylabelvisible = true,
            yticklabelsvisible = true,
            yticksvisible = true
        )
        if all(source -> source.resolutions[facet.request_index].clip &&
                source.resolutions[facet.request_index].kind === :declared_floor, published)
            if all(ismissing, yvalues)
                attributes = merge(attributes, (subtitle="Undefined phase",))
            end
        end
        axis,scales = _addon_axis!(
            panel.content,
            xobservation,
            yobservation;
            title = _semantic_panel_title(
                panel_titles, object, facet, panel_index, length(page.facets)),
            xscale,
            yscale,
            attributes,
            native_attributes=shell.axis_attributes
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
                dependent_plots,
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
        reset! = _addon_reset!(axis, series)
        push!(axes, axis)
        push!(panels, panel)
        push!(panel_group_labels, scoped_labels)
        push!(resets, reset!)
        push!(requested_scales,scales)
    end
    mode === :paired && length(axes) > 1 &&
        _addon_responsive_axis_grid!(
            shell.figure, shell.canvas, panels, axes, page.dimensions)
    local_panel_legends = if blocked
        selected = Set((facet.row, facet.column) for facet in page.facets)
        Tuple(pair for pair in _addon_panel_legend_pairs(panel_legends) if first(pair) in selected)
    else
        panel_legends
    end
    built = _addon_finish!(
        shell, axes, resets, groups, group_order, group_labels;
        requested_scales, signed_ylog,
        dependent_plots,
        series_attributes,
        series_defaults,
        title,
        figure_title,
        title_attributes,
        legend_position,
        legend_anchor,
        legend_attributes,
        legend_overflow,
        legend_title,
        panels,
        panel_legends=local_panel_legends,
        panel_group_labels,
        controls,
        display_plot,
        export_name = title,
        export_theme,
        open_export
    )
    if blocked
        built.addon_state = merge(built.addon_state, (matrix_block=(;
            index=page.block, dimensions=page.dimensions,
            cells=Tuple(values(cells)),
            coordinates=Tuple((facet.row, facet.column) for facet in page.facets)),))
    end
    return built
end

function _addon_line_pages(
        sources::Tuple;
        publications = nothing,
        frequencies = nothing,
        ydata,
        series_labels = nothing,
        series_indices = collect(eachindex(sources)),
        series_family_labels = nothing,
        formulation_sources = nothing,
        formulation_roles = nothing,
        series_attributes = nothing,
        series_defaults = nothing,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        panel_titles = nothing,
        freq_unit = :base,
        length_unit = :kilo,
        quantity_units = nothing,
        clip::Bool = true,
        atol = nothing,
        fig_size = nothing,
        layout = nothing,
        blocks = nothing,
        xscale = :linear,
        yscale = :linear,
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_attributes::NamedTuple = (;),
        legend_overflow::Symbol = :ellipsis,
        panel_legends = (),
        signed_ylog::Bool = false,
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true,
        kwargs...
)
    _addon_activate_backend(backend)
    legend_overflow in (:ellipsis, :show_all) || throw(ArgumentError(
        "legend_overflow must be :ellipsis or :show_all",
    ))
    legend_attributes isa NamedTuple || throw(ArgumentError(
        "legend_attributes must be a NamedTuple",
    ))
    explicit_source_labels = series_labels !== nothing
    source_labels = explicit_source_labels ?
                    _comparison_labels(series_labels, length(sources)) :
                    Tuple("Result $index" for index in eachindex(sources))
    published = if publications !== nothing
        length(publications)==length(sources) || throw(DimensionMismatch("one publication is required per source"))
        publications
    elseif length(sources) == 1
        (_prepare_line_observations(
            only(sources);
            frequencies,
            ydata,
            freq_unit,
            length_unit,
            quantity_units,
            clip,
            atol
        ),)
    else
        _prepare_line_comparison(
            sources;
            ydata,
            series_labels = source_labels,
            freq_unit,
            length_unit,
            quantity_units,
            clip,
            atol,
            frequencies
        ).published
    end
    any(source -> length(source.frequency.values) <= 1, published) &&
        return LineCableModels.UIPlot[]

    facets = _semantic_line_facets(first(published), ydata)
    mode, pages = _semantic_line_pages(first(sources), facets, layout, blocks)
    if blocks !== nothing
        selected = Set((facet.row, facet.column) for facet in facets)
        for (position, _) in _addon_panel_legend_pairs(panel_legends)
            position in selected || throw(ArgumentError("panel legend $position is not a selected matrix coordinate"))
        end
        if panel_titles isa Union{Tuple,AbstractVector}
            length(panel_titles) == length(facets) || throw(DimensionMismatch(
                "blocked panel_titles must contain one title per selected matrix facet"))
            panel_titles = Dict((facet.identity, facet.row, facet.column) => String(label)
                for (facet, label) in zip(facets, panel_titles))
        end
    end
    effective_legend_position = length(sources) > 1 || explicit_source_labels ?
                                legend_position : nothing
    built = LineCableModels.UIPlot[]
    styles = formulation_sources === nothing ? nothing :
        LineCableModels.PlotBuilder._series_attributes(series_attributes,length(sources))
    for (page_index, page) in enumerate(pages)
        page_labels = series_family_labels === nothing ? source_labels :
            _comparison_labels(series_family_labels[first(page.facets).family],length(sources))
        retained = collect(eachindex(sources))
        page_defaults = series_defaults
        page_attributes = series_attributes
        if formulation_sources !== nothing
            quantities = unique(facet.quantity for facet in page.facets)
            identities = [Tuple(LineCableModels.formula_id(source,quantity) for quantity in quantities)
                for source in formulation_sources]
            retained = unique(eachindex(sources)) do index
                formulation_roles[index] === :reference || any(ismissing,identities[index]) ?
                    index : identities[index]
            end
            for index in setdiff(eachindex(sources),retained)
                representative = only(filter(other -> formulation_roles[other] !== :reference &&
                    isequal(identities[other],identities[index]),retained))
                left,right = published[representative],published[index]
                same = isequal(left.frequency.values,right.frequency.values)
                for request in unique(facet.request_index for facet in page.facets)
                    a,b = left.observations[request].values,right.observations[request].values
                    same &= isequal(left.coordinates[request],right.coordinates[request]) && size(a)==size(b) &&
                        all(zip(a,b)) do (a,b)
                            isequal(a,b) || a isa Number && b isa Number &&
                                isapprox(LineCableModels.nominal(a),LineCableModels.nominal(b)) &&
                                isapprox(LineCableModels.uncertainty(a),LineCableModels.uncertainty(b))
                        end
                end
                same || throw(ArgumentError(
                    "repeated formulation for $(first(page.facets).quantity) has conflicting saved observations; inspect the calculations separately"))
            end
            defaults = _addon_default_formulations(formulation_sources;quantity=first(quantities))
            roles = Tuple(formulation_roles[index] === :reference ? :reference :
                defaults[index] ? :default : :alternative for index in retained)
            page_defaults = _addon_comparison_styles(Tuple(series_indices[index] for index in retained),
                roles,maximum(series_indices))
            page_attributes = Tuple(styles[index] for index in retained)
        end
        automatic_title = _semantic_page_title(first(sources), page, mode)
        page_title = title === nothing ? automatic_title : String(title)
        (length(pages) > 1 || blocks !== nothing) && title !== nothing &&
            (page_title = "$page_title — $automatic_title")
        visible_title = _semantic_page_option(
            figure_title, page_index, length(pages), "figure_title")
        push!(built,
            with_theme(_addon_theme(export_theme = export_theme)) do
                _addon_semantic_line_page(
                    first(sources), Tuple(published[index] for index in retained),
                    Tuple(page_labels[index] for index in retained), page, mode;
                    series_indices=Tuple(series_indices[index] for index in retained),
                    series_defaults=page_defaults,
                    series_attributes=page_attributes,
                    title = page_title,
                    figure_title = visible_title,
                    title_attributes,
                    panel_titles,
                    fig_size,
                    xscale,
                    yscale,
                    legend_position = effective_legend_position,
                    legend_anchor,
                    legend_title,
                    legend_attributes,
                    legend_overflow,
                    panel_legends,
                    signed_ylog,
                    controls,
                    display_plot=false,
                    export_theme,
                    open_export,
                    kwargs...
                )
            end)
    end
    blocks === nothing || _addon_equal_matrix_cells!(built)
    if display_plot
        for page in built
            _addon_display!(page.figure, page.export_name)
        end
    end
    return length(built) == 1 ? only(built) : built
end

function _addon_semantic_line_plots(
        object;
        frequencies = nothing,
        ydata,
        series_labels = nothing,
        kwargs...
)
    return _addon_line_pages(
        (object,);
        frequencies,
        ydata,
        series_labels,
        kwargs...
    )
end
