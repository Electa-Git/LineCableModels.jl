function _addon_preview_axis!(
        shell,
        position,
        title,
        polygons,
        references,
        limits,
        groups,
        group_order,
        group_labels
)
    unit = LineCableModels.Units.units(:base, :meter)
    unit_label = LineCableModels.Units.label(unit)
    panel = _addon_panel!(shell, position)
    axis = Axis(
        panel.content;
        title,
        xlabel = "y [$unit_label]",
        ylabel = "z [$unit_label]",
        aspect = DataAspect(),
        tellwidth = false,
        tellheight = false
    )
    for reference in references
        plot = hlines!(
            axis,
            reference.values;
            color = reference.color,
            linewidth = reference.width
        )
        if !haskey(groups, reference.group)
            groups[reference.group] = Any[]
            push!(group_order, reference.group)
        end
        push!(groups[reference.group], plot)
    end
    for polygon in polygons
        plot = poly!(
            axis,
            polygon.geometry;
            label = polygon.label,
            color = polygon.color,
            strokecolor = polygon.stroke,
            strokewidth = polygon.width
        )
        if !haskey(groups, polygon.group)
            groups[polygon.group] = Any[]
            push!(group_order, polygon.group)
        end
        push!(groups[polygon.group], plot)
        polygon.label === nothing || (group_labels[polygon.group] = polygon.label)
    end
    reset! = if limits === nothing
        () -> autolimits!(axis)
    else
        () -> begin
            xlims!(axis, limits[1]...)
            ylims!(axis, limits[2]...)
            axis
        end
    end
    reset!()
    return axis, reset!, panel
end

function _addon_preview_finish!(
        shell,
        axes,
        resets,
        groups,
        group_order,
        group_labels;
        title,
        figure_title = nothing,
        title_attributes = (;),
        panels,
        panel_legends,
        panel_legend_titles = nothing,
        display_legend,
        legend_position,
        legend_anchor,
        legend_title,
        legend_attributes,
        legend_overflow,
        color_scales,
        colorbar_position,
        colorbar_attributes,
        controls,
        display_plot,
        export_name,
        export_theme,
        open_export
)
    return _addon_finish!(
        shell,
        axes,
        resets,
        Function[],
        Function[],
        groups,
        group_order,
        group_labels;
        title,
        figure_title,
        title_attributes,
        legend_position = display_legend ? legend_position : nothing,
        legend_anchor,
        legend_title,
        legend_attributes,
        legend_overflow,
        panels,
        panel_legends,
        panel_legend_titles,
        color_scales,
        colorbar_position,
        colorbar_attributes,
        controls,
        display_plot,
        export_name,
        export_theme,
        open_export
    )
end

function _addon_preview(
        design::LineCableModels.DataModel.CableDesign;
        x_offset::Real = 0.0,
        y_offset::Real = 0.0,
        display_legend::Bool = true,
        display_id::Bool = false,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        panel_titles = nothing,
        display_colorbars::Bool = true,
        size::Tuple{Int, Int} = (900, 700),
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_attributes::NamedTuple = (; nbanks = 2),
        legend_overflow::Symbol = :ellipsis,
        panel_legends = (),
        legend_group = nothing,
        legend_labels = nothing,
        colorbar_position = :right,
        colorbar_attributes::NamedTuple = (; vertical = false),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    _addon_activate_backend(backend)
    isfinite(x_offset) && isfinite(y_offset) || throw(ArgumentError(
        "preview offsets must be finite",
    ))
    display_title = title === nothing ?
                    _native_cable_title(display_id, design) :
                    String(title)
    resolved_panel_titles = _addon_panel_titles(panel_titles, nothing, 1)
    panel_title = resolved_panel_titles === nothing ?
                  display_title : only(resolved_panel_titles)
    polygons = _native_design_shapes(
        design,
        x_offset,
        y_offset;
        # Retain presentation metadata even when the initial legend is hidden;
        # `figurelegend!` may place it later without rebuilding geometry.
        display_legend = true,
        legend_group,
        legend_labels
    )
    color_scales = display_colorbars ?
                   _material_schemes(
        LineCableModels.DataModel.material_property_ranges(design)
    ) : ()
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls)
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, String}()
        axis, reset!, panel = _addon_preview_axis!(
            shell,
            (1, 1),
            panel_title,
            polygons,
            (),
            nothing,
            groups,
            order,
            labels
        )
        _addon_center_aspect_canvas!(shell)
        _addon_preview_finish!(
            shell,
            Any[axis],
            Function[reset!],
            groups,
            order,
            labels;
            title = display_title,
            figure_title,
            title_attributes,
            panels = (panel,),
            panel_legends,
            display_legend,
            legend_position,
            legend_anchor,
            legend_title,
            legend_attributes,
            legend_overflow,
            color_scales,
            colorbar_position = display_colorbars ? colorbar_position : nothing,
            colorbar_attributes,
            controls,
            display_plot,
            export_name = design.cable_id,
            export_theme,
            open_export
        )
    end
end

function _addon_preview(
        designs::AbstractVector{<:LineCableModels.DataModel.CableDesign};
        layout = nothing,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        panel_titles = nothing,
        display_colorbars::Bool = true,
        size::Tuple{Int, Int} = (1200, 900),
        colorbar_position = :right,
        colorbar_attributes::NamedTuple = (; vertical = false),
        panel_legends = (),
        legend_group = nothing,
        legend_labels = nothing,
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    _addon_activate_backend(backend)
    rows, columns = _native_preview_layout(
        length(designs), layout)
    color_scales = display_colorbars ?
                   _material_schemes(
        LineCableModels.DataModel.material_property_ranges(designs)
    ) : ()
    display_title = title === nothing ? "Cable design previews" : String(title)
    resolved_panel_titles = _addon_panel_titles(panel_titles, nothing, length(designs))
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls)
        axes = Any[]
        panels = Any[]
        resets = Function[]
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, String}()
        for (index, design) in enumerate(designs)
            polygons = _native_design_shapes(
                design,
                0.0,
                0.0;
                display_legend = true,
                legend_group,
                legend_labels
            )
            axis, reset!, panel = _addon_preview_axis!(
                shell,
                (cld(index, columns), mod1(index, columns)),
                resolved_panel_titles === nothing ?
                design.cable_id : resolved_panel_titles[index],
                polygons,
                (),
                nothing,
                groups,
                order,
                labels
            )
            push!(axes, axis)
            push!(panels, panel)
            push!(resets, reset!)
        end
        for row in 1:rows
            rowsize!(shell.canvas, row, Relative(1 / rows))
        end
        for column in 1:columns
            colsize!(shell.canvas, column, Relative(1 / columns))
        end
        _addon_center_aspect_canvas!(shell)
        _addon_preview_finish!(
            shell,
            axes,
            resets,
            groups,
            order,
            labels;
            title = display_title,
            figure_title,
            title_attributes,
            panels,
            panel_legends,
            display_legend = false,
            legend_position = nothing,
            legend_anchor = :rt,
            legend_title = nothing,
            legend_attributes = (;),
            legend_overflow = :show_all,
            color_scales,
            colorbar_position = display_colorbars ? colorbar_position : nothing,
            colorbar_attributes,
            controls,
            display_plot,
            export_name = "cable_design_previews",
            export_theme,
            open_export
        )
    end
end

function _addon_preview(
        system::LineCableModels.DataModel.LineCableSystem;
        earth_model = nothing,
        zoom_factor = nothing,
        display_legend::Bool = true,
        display_id::Bool = false,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        panel_titles = nothing,
        display_colorbars::Bool = true,
        size::Tuple{Int, Int} = (900, 700),
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_attributes::NamedTuple = (;),
        legend_overflow::Symbol = :ellipsis,
        panel_legends = (),
        legend_group = nothing,
        legend_labels = nothing,
        colorbar_position = :right,
        colorbar_attributes::NamedTuple = (; vertical = false),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    _addon_activate_backend(backend)
    limits = _native_system_limits(system, zoom_factor)
    polygons, references = _native_system_shapes(
        system,
        earth_model,
        limits,
        display_legend || !isempty(_addon_panel_legend_pairs(panel_legends));
        legend_group,
        legend_labels
    )
    color_scales = display_colorbars ?
                   _native_earth_colorbars(earth_model) : ()
    display_title = title === nothing ?
                    _native_system_title(display_id, system) :
                    String(title)
    resolved_panel_titles = _addon_panel_titles(panel_titles, nothing, 1)
    panel_title = resolved_panel_titles === nothing ?
                  display_title : only(resolved_panel_titles)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls)
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, String}()
        axis, reset!, panel = _addon_preview_axis!(
            shell,
            (1, 1),
            panel_title,
            polygons,
            references,
            limits,
            groups,
            order,
            labels
        )
        _addon_center_aspect_canvas!(shell)
        _addon_preview_finish!(
            shell,
            Any[axis],
            Function[reset!],
            groups,
            order,
            labels;
            title = display_title,
            figure_title,
            title_attributes,
            panels = (panel,),
            panel_legends,
            display_legend,
            legend_position,
            legend_anchor,
            legend_title,
            legend_attributes,
            legend_overflow,
            color_scales,
            colorbar_position = display_colorbars ? colorbar_position : nothing,
            colorbar_attributes,
            controls,
            display_plot,
            export_name = system.system_id,
            export_theme,
            open_export
        )
    end
end

function _addon_material_scale(;
        size::Tuple{Int, Int} = (800, 400),
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        colorbar_position = nothing,
        colorbar_attributes::NamedTuple = (; vertical = false),
        backend = nothing,
        display_plot::Bool = true,
        controls::Bool = true,
        export_theme::Symbol = :default,
        open_export::Bool = true
)
    _addon_activate_backend(backend)
    title = "Material property colour scale"
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls)
        use_canvas = colorbar_position === nothing
        scale_canvas = if use_canvas
            grid = GridLayout(
                ; width = Relative(1), height = Relative(1),
                tellwidth = false, tellheight = false
            )
            shell.canvas[1, 1] = grid
            grid
        else
            nothing
        end
        _addon_finish!(
            shell,
            Any[],
            Function[],
            Function[],
            Function[],
            Dict{Symbol, Vector{Any}}(),
            Symbol[],
            Dict{Symbol, String}();
            title,
            figure_title,
            title_attributes,
            legend_position = nothing,
            legend_attributes = (;),
            legend_overflow = :show_all,
            color_scales = _material_schemes(
                LineCableModels.DataModel.material_property_ranges()
            ),
            colorbar_position = use_canvas ? :right : colorbar_position,
            colorbar_attributes,
            colorbar_target = use_canvas ? scale_canvas[1, 1] : nothing,
            colorbar_target_orientation = use_canvas ? :horizontal : nothing,
            controls,
            display_plot,
            export_name = "material_scale",
            export_theme,
            open_export
        )
    end
end
