function _publication_observation_name(publication, observation)
    for name in keys(publication.metadata.observation_columns)
        contract = getproperty(publication.metadata.observation_columns, name)
        contract.quantity == observation.quantity && contract.unit == observation.unit &&
            return name
    end
    return nothing
end

function _publication_frequency_name(publication, observation_name)
    for name in publication.metadata.row_order
        name == observation_name && continue
        hasproperty(publication.columns, name) || continue
        hasproperty(publication.metadata.observation_columns, name) || continue
        contract = getproperty(publication.metadata.observation_columns, name)
        contract.quantity isa LineCableModels.Units.Quantity{:frequency} && return name
    end
    return nothing
end

function _publication_series_data(publication, observation)
    observation_name = _publication_observation_name(publication, observation)
    yvalues = observation_name === nothing ? collect(vec(observation.values)) :
              collect(getproperty(publication.columns, observation_name))
    frequency_name = _publication_frequency_name(publication, observation_name)
    xvalues = frequency_name === nothing ? collect(eachindex(yvalues)) :
              collect(getproperty(publication.columns, frequency_name))
    length(xvalues) == length(yvalues) || throw(DimensionMismatch(
        "publication plot coordinates must align with observation values",
    ))
    coordinate_names = frequency_name === nothing ? () :
                       Tuple(name
    for name in publication.metadata.row_order
    if name != frequency_name && name != observation_name &&
       hasproperty(publication.columns, name))
    grouped = Dict{Tuple, Vector{Int}}()
    order = Tuple[]
    for index in eachindex(yvalues)
        coordinate = Tuple(getproperty(publication.columns, name)[index]
        for name in coordinate_names)
        haskey(grouped, coordinate) || begin
            grouped[coordinate] = Int[]
            push!(order, coordinate)
        end
        push!(grouped[coordinate], index)
    end
    isempty(order) && (grouped[()] = collect(eachindex(yvalues)); push!(order, ()))
    xobservation = if frequency_name === nothing
        (; values = xvalues, quantity = nothing, unit = nothing)
    else
        contract = getproperty(publication.metadata.observation_columns, frequency_name)
        (; values = xvalues, quantity = contract.quantity, unit = contract.unit)
    end
    return (;
        xobservation,
        yobservation = merge(observation, (; values = yvalues)),
        coordinate_names,
        series = Tuple((;
                           coordinate,
                           x = xvalues[grouped[coordinate]],
                           y = yvalues[grouped[coordinate]]
                       ) for coordinate in order),
        xlabel = frequency_name === nothing ? "Sample" : nothing
    )
end

function _publication_coordinate_label(names, coordinate)
    isempty(coordinate) && return ""
    names == (:row, :column) && return "[$(coordinate[1]),$(coordinate[2])]"
    names == (:row,) && return "[$(only(coordinate))]"
    return join(("$(name)=$(value)" for (name, value) in zip(names, coordinate)), ", ")
end

function _publication_group(observation, coordinate, index::Int)
    family = LineCableModels.Units.family(observation.quantity)
    prefix = if family === Val(:series) || family === Val(:shunt)
        LineCableModels.Units.symbol(LineCableModels.Units.quantity(_family_parent(family)))
    else
        LineCableModels.Units.symbol(observation.quantity)
    end
    isempty(coordinate) && return Symbol("$(prefix)_observation_$index"), prefix
    suffix = join(string.(coordinate), "_")
    return Symbol("$(prefix)_$suffix"), prefix
end

function _addon_publication_plot(
        publication::LineCableModels.Grammar.ObservationPublication;
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        panel_titles = nothing,
        fig_size::Tuple{Int, Int} = (800, 400),
        layout = nothing,
        xscale::Symbol = :linear,
        yscale::Symbol = :linear,
        display_legend::Bool = false,
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_labels = nothing,
        legend_attributes::NamedTuple = (;),
        legend_overflow::Symbol = :ellipsis,
        panel_legends = (),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    isempty(publication.observations) && throw(ArgumentError(
        "an observation publication plot requires at least one observation",
    ))
    resolved_panel_titles = _addon_panel_titles(panel_titles, length(publication))
    display_title = title === nothing ? "Observation publication" : String(title)
    _addon_activate_backend(backend)
    positions, dimensions = _addon_positions(length(publication), layout)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size = fig_size, controls)
        shell.canvas.default_rowgap = Fixed(24)
        shell.canvas.default_colgap = Fixed(64)
        rowgap!(shell.canvas, 24)
        colgap!(shell.canvas, 64)
        axes = Any[]
        panels = Any[]
        resets = Function[]
        xsetters = Function[]
        ysetters = Function[]
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        group_labels = Dict{Symbol, String}()
        group_colors = Dict{Symbol, Any}()
        panel_group_labels = Any[]
        for (index, observation) in enumerate(publication)
            data = _publication_series_data(publication, observation)
            xvalues = data.xobservation.values
            yvalues = data.yobservation.values
            xscales = _axis_scales(xvalues)
            yscales = _axis_scales(yvalues)
            xscale in xscales || throw(DomainError(
                xvalues,
                "logarithmic sample axes require positive finite values"
            ))
            yscale in yscales || throw(DomainError(
                yvalues,
                "logarithmic ordinate axes require positive finite data and uncertainty bounds"
            ))
            panel_title = resolved_panel_titles === nothing ?
                          LineCableModels.Units.label(observation.quantity) :
                          resolved_panel_titles[index]
            row, column = positions[index]
            panel = _addon_panel!(shell, (row, column))
            attributes = (;
                xlabelvisible = row == dimensions[1],
                xticklabelsvisible = row == dimensions[1],
                xticksvisible = row == dimensions[1]
            )
            axis,
            axis_labels = _addon_axis!(
                panel.content,
                data.xobservation,
                data.yobservation;
                title = panel_title,
                xlabel = data.xlabel,
                xscale,
                yscale,
                xscales,
                yscales,
                attributes
            )
            series = NamedTuple[]
            scoped_labels = Dict{Symbol, String}()
            quantity_symbol = LineCableModels.Units.symbol(observation.quantity)
            for item in data.series
                coordinate_label = _publication_coordinate_label(
                    data.coordinate_names, item.coordinate)
                group,
                parent_symbol = _publication_group(
                    observation, item.coordinate, index)
                if !haskey(groups, group)
                    groups[group] = Any[]
                    push!(order, group)
                    group_labels[group] = isempty(coordinate_label) ?
                                          String(parent_symbol) :
                                          "$(parent_symbol)$(coordinate_label)"
                    group_colors[group] = _addon_comparison_color(length(order))
                end
                local_label = isempty(coordinate_label) ? panel_title :
                              "$(quantity_symbol)$(coordinate_label)"
                plots = _addon_line!(
                    axis,
                    item.x,
                    item.y;
                    label = local_label,
                    color = group_colors[group]
                )
                append!(groups[group], plots)
                scoped_labels[group] = local_label
                push!(series, (; xdata = item.x, ydata = item.y, plots))
            end
            push!(panel_group_labels, scoped_labels)
            reset! = () -> _addon_reset!(axis, series, axis_labels)
            xsetter = scale -> _addon_set_axis!(
                axis, :x, axis_labels.x, xscales,
                something(_addon_scientific_exponent(xvalues), 0), scale)
            ysetter = scale -> _addon_set_axis!(
                axis, :y, axis_labels.y, yscales,
                something(_addon_scientific_exponent(yvalues), 0), scale)
            push!(axes, axis)
            push!(panels, panel)
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
        _addon_relabel_legend!(group_labels, groups, order, legend_labels)
        if legend_labels !== nothing
            for scoped_labels in panel_group_labels, group in keys(scoped_labels)

                haskey(group_labels, group) && (scoped_labels[group] = group_labels[group])
            end
        end
        layout === nothing && _addon_responsive_axis_grid!(
            shell.figure,
            shell.canvas,
            panels,
            axes,
            dimensions
        )
        _addon_finish!(
            shell,
            axes,
            resets,
            xsetters,
            ysetters,
            groups,
            order,
            group_labels;
            title = display_title,
            figure_title,
            title_attributes,
            legend_position = display_legend ? legend_position : nothing,
            legend_anchor,
            legend_title,
            legend_attributes,
            legend_overflow,
            panels,
            panel_legends,
            panel_group_labels,
            controls,
            display_plot,
            export_name = display_title,
            export_theme,
            open_export
        )
    end
end
