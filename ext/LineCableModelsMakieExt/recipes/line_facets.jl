# Line faceting --------------------------------------------------------------
#
# The recipe supplies physical coordinate descriptions. Shared plotting chooses
# panel and curve identities before native blocks are constructed.
function _coordinate_labels(product, coordinates, orientation)
    c=product.coordinates
    rows, columns, _=coordinates
    if orientation===:rows
        name=get(c,:domain,nothing)===:ModalDomain ? "Mode" :
             get(c,:domain,nothing)===:PhaseDomain ? "Conductor" : "Row"
        return ["$name $(c.labels[row])" for row in rows]
    elseif c.kind===:vector
        return ["$(c.axis_label) $(c.labels[row])" for row in rows]
    elseif c.kind===:matrix && get(c, :column_domain, nothing)===:ModalDomain
        return ["Conductor $(c.labels[row]), Mode $(c.column_labels[column])"
                for row in rows for column in columns]
    elseif c.kind===:matrix
        return ["$(c.labels[row]) → $(c.labels[column])"
                for row in rows for column in columns]
    elseif c.kind===:diagonal
        return ["Conductor $(c.labels[row])" for row in rows]
    elseif c.kind===:assemblies
        return string.(c.labels[c.assemblies])
    end
    return string.(1:(length(rows) * length(columns)))
end

function _quantity_panel_title(facet)
    label=Units.label(facet.quantity)
    facet.identity isa Tuple && first(facet.identity)===LineCableModels.statistics &&
        (label *= " · " * string(last(facet.identity)))
    if facet.orientation===:rows
        domain=facet.column_domain===nothing ? facet.domain : facet.column_domain
        coordinate=domain===:ModalDomain ? "Mode" : domain===:PhaseDomain ? "Conductor" : "Column"
        panel_label=facet.column_domain===:ModalDomain ? Units.symbol(facet.quantity) : label
        return "$panel_label, $coordinate $(facet.column_label)"
    end
    if facet.orientation===:coordinates
        panel_label=facet.column_domain===:ModalDomain ? Units.symbol(facet.quantity) :
                    label
        return isempty(facet.source_label) ? panel_label :
               "$panel_label · $(facet.source_label)"
    end
    facet.kind in (:matrix, :diagonal, :vector) || return label
    facet.column_domain===:ModalDomain &&
        return "$(Units.symbol(facet.quantity)), mode $(facet.column_label) → $(facet.row_label)"
    is_modal = facet.domain === :ModalDomain
    coordinate = is_modal ? "mode" : "conductor"

    facet.kind===:vector && return "$label, $coordinate $(facet.row_label)"

    if facet.row == facet.column
        return "$label, $coordinate $(facet.row)"
    elseif is_modal
        return "$label, intermodal residual ($(facet.row), $(facet.column))"
    else
        return "$label, conductors $(facet.row) → $(facet.column)"
    end
end

function _quantity_page_title(page)
    facet=first(page.facets)
    label=Units.label(facet.quantity)
    facet.identity isa Tuple && first(facet.identity)===LineCableModels.statistics &&
        (label *= " · " * string(last(facet.identity)))
    facet.orientation===:rows && !isempty(facet.source_label) &&
        (label *= " · " * facet.source_label)
    return page.index==(1, 1) ? label : "$label ($(page.index[1]),$(page.index[2]))"
end

function _panel_title(panel_titles, facet)
    panel_titles === nothing && return _quantity_panel_title(facet)
    panel_titles isa Function && return panel_titles(facet)
    if panel_titles isa AbstractDict
        quantity_symbol = Symbol(Units.symbol(facet.quantity))
        candidates = (
            (facet.identity, facet.panel_identity),
            facet.panel_identity,
            (facet.identity, facet.row, facet.column),
            (quantity_symbol, facet.row, facet.column),
            (facet.row, facet.column),
            facet.identity,
            quantity_symbol
        )
        for candidate in candidates
            haskey(panel_titles, candidate) && return panel_titles[candidate]
        end
        return _quantity_panel_title(facet)
    end
    throw(ArgumentError("panel_titles must be a dictionary, function, or nothing"))
end

function _addon_line_page(
        published,
        page,
        ;
        series_indices,
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
        legend_cap,
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

    for (facet, position) in zip(page.facets, page.positions)
        observation = first(published).observations[facet.request_index]
        curve_records=map(facet.curves) do curve
            source=published[curve.source_index]
            selected=curve.local_sample===nothing ? Colon() :
                     curve.local_sample:curve.local_sample
            x=source.xdata[facet.request_index].values[selected]
            y=collect(view(source.observations[facet.request_index].values,
                curve.local_row, curve.local_column, selected))
            errors=get(source.observations[facet.request_index], :errors, nothing)
            yerror=errors===nothing ? nothing :
                   collect(view(errors, curve.local_row, curve.local_column, selected))
            (; x, y, yerror)
        end
        xvalues = collect(Iterators.flatten(
            record.x for record in curve_records))
        yvalues = collect(Iterators.flatten(record.y for record in curve_records))
        xobservation = merge(first(published).xdata[facet.request_index], (;
            values = xvalues))
        yobservation = merge(observation, (; values = yvalues))
        panel=merge(cells[position], (; logical_position = facet.panel_identity))
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
            curve -> begin
                resolution=published[curve.source_index].resolutions[facet.request_index]
                resolution.clip && resolution.kind===:declared_floor
            end, facet.curves)
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
            title = _panel_title(panel_titles, facet),
            xscale = xscale===nothing ?
                     (facet.kind in (:matrix, :diagonal, :vector) ? :log10 : :linear) :
                     xscale,
            yscale,
            xlabel = get(xobservation, :label, nothing),
            ylabel = facet.column_domain===:ModalDomain ?
                     string(Units.symbol(facet.quantity),
                isempty(Units.label(yobservation.unit)) ? "" :
                " [$(Units.label(yobservation.unit))]") : nothing,
            attributes = facet.kind===:assemblies ?
                         merge(attributes, (dim1_conversion = Makie.CategoricalConversion(),)) :
                         attributes,
            native_attributes = shell.axis_attributes
        )
        series = NamedTuple[]
        scoped_labels = Dict{Symbol, Any}()
        for (curve_identity, record) in zip(facet.curves, curve_records)
            style_index=findfirst(==(curve_identity.slot), series_indices)
            curve = record.y
            yerror = record.yerror
            group = curve_identity.group
            curve_label = curve_identity.label
            interval_support=Dict{Makie.Plot, Any}()
            draw! = facet.kind in (:assemblies, :array) ? _addon_points! : _addon_line!
            plots = draw!(
                axis,
                record.x,
                curve;
                dependent_plots,
                label = curve_label,
                color = series_defaults[style_index].attributes.color,
                phase = series_defaults[style_index].phase,
                endpoints = series_defaults[style_index].endpoints,
                marker_coordinates,
                errorbar_sampling, yerror, interval_support
            )
            if !haskey(groups, group)
                groups[group] = Any[]
                push!(group_order, group)
                group_labels[group] = curve_label
            end
            append!(groups[group], plots)
            scoped_labels[group] = curve_label
            binding=(;
                xdata = facet.kind===:assemblies ? nothing :
                        record.x,
                ydata = curve,
                yerror,
                interval_support,
                sampled_intervals = facet.kind in (:matrix, :diagonal, :vector) &&
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
    selected=Set(facet.panel_identity for facet in page.facets)
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
        legend_cap,
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
                coordinates = Tuple(f.panel_identity for f in page.facets)),
            page_cells = cells))
    return built
end
