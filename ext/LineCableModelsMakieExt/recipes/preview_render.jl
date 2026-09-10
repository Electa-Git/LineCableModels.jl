function _addon_preview_axis!(
        shell,
        position,
        title,
        polygons,
        references,
        limits,
        groups,
        group_order,
        group_labels;
        earth_model = nothing,
        display_surface_gradient::Bool = true
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
    earth_spans = NamedTuple[]
    surface_gradient = nothing
    if earth_model !== nothing && !earth_model.vertical_layers
        top = 0.0
        for (index, layer) in enumerate(earth_model.layers[2:end])
            bottom = top - nominal(layer.thickness)
            span = hspan!(axis, 0.0, 0.0;
                color = _material_color(layer; alpha = 0.25),
                strokewidth = 0, xautolimits = false, yautolimits = false)
            translate!(span, 0, 0, -100)
            group = Symbol("earth_$index")
            groups[group] = Any[span]
            push!(group_order, group)
            group_labels[group] = "Earth layer $index"

            if isfinite(bottom)
                boundary = hlines!(axis, [bottom];
                    color = RGBA(0.25, 0.27, 0.30, 0.25), linewidth = 0.6,
                    visible = span.visible, xautolimits = false, yautolimits = false)
                translate!(boundary, 0, 0, -80)
            end
            push!(earth_spans, (; top, bottom, span))
            top = bottom
        end
        if display_surface_gradient
            # One sky decoration, anchored at the surface and independent of
            # the soil layers. Transparency reveals the background near z=0.
            surface_gradient = hspan!(axis, zeros(64), zeros(64);
                color = [RGBA(0.55, 0.76, 0.90, 0.45 * ((i - 1) / 63)^1.3)
                         for i in 1:64],
                strokewidth = 0, xautolimits = false, yautolimits = false)
            translate!(surface_gradient, 0, 0, -90)
        end
    end
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
    if !isempty(earth_spans)
        # HSpan owns full-width coverage. This scene-owned callback clips only
        # the physical vertical extents, keeping infinite coordinates out of
        # rendering and all background geometry out of automatic limits.
        on(axis.scene, axis.finallimits; update = true) do view
            lower = view.origin[2]
            upper = lower + view.widths[2]
            for entry in earth_spans
                bottom = clamp(entry.bottom, lower, upper)
                top = clamp(entry.top, lower, upper)
                Makie.update!(entry.span, bottom, top)
            end
            if surface_gradient !== nothing
                # Increase the blue from the transparent physical surface to
                # the top of the view. Underground-only views contain no sky.
                sky_height = max(upper, 0)
                bands = length(surface_gradient[1][])
                lows = [clamp((i - 1) * sky_height / bands, lower, upper) for i in 1:bands]
                highs = [clamp(i * sky_height / bands, lower, upper) for i in 1:bands]
                Makie.update!(surface_gradient, lows, highs)
                surface_gradient.visible[] = upper > 0
            end
        end
    end
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
        series_attributes = nothing,
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
        series_attributes,
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
        display_dielectric_pattern::Bool = true,
        display_legend::Bool = true,
        display_id::Bool = false,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        series_attributes = nothing,
        panel_titles = nothing,
        display_colorbars::Bool = true,
        size::Tuple{Int, Int} = (900, 700),
        legend_position = :right,
        legend_anchor = :rt,
        legend_title = nothing,
        legend_attributes::NamedTuple = (; nbanks = 1),
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
    resolved_panel_titles = _addon_panel_titles(panel_titles, 1)
    panel_title = resolved_panel_titles === nothing ?
                  display_title : only(resolved_panel_titles)
    polygons = _native_design_shapes(
        design,
        x_offset,
        y_offset;
        # Retain presentation metadata even when the initial legend is hidden;
        # `figurelegend!` may place it later without rebuilding geometry.
        display_legend = true,
        display_dielectric_pattern,
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
        axis, reset!,
        panel = _addon_preview_axis!(
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
            series_attributes,
            export_theme,
            open_export
        )
    end
end

function _addon_preview(
        designs::AbstractVector{<:LineCableModels.DataModel.CableDesign};
        layout = nothing,
        display_dielectric_pattern::Bool = true,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        series_attributes = nothing,
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
    resolved_panel_titles = _addon_panel_titles(panel_titles, length(designs))
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
                display_dielectric_pattern,
                legend_group,
                legend_labels
            )
            axis, reset!,
            panel = _addon_preview_axis!(
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
            series_attributes,
            export_theme,
            open_export
        )
    end
end

function _addon_preview(
        system::LineCableModels.DataModel.LineCableSystem;
        earth_model = nothing,
        zoom_factor = nothing,
        display_dielectric_pattern::Bool = true,
        display_surface_gradient::Bool = true,
        display_legend::Bool = true,
        display_id::Bool = false,
        title = nothing,
        figure_title = nothing,
        title_attributes::NamedTuple = (;),
        series_attributes = nothing,
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
    polygons,
    references = _native_system_shapes(
        system,
        display_legend || !isempty(_addon_panel_legend_pairs(panel_legends));
        display_dielectric_pattern,
        legend_group,
        legend_labels
    )
    color_scales = display_colorbars ?
                   _native_earth_colorbars(earth_model) : ()
    display_title = title === nothing ?
                    _native_system_title(display_id, system) :
                    String(title)
    resolved_panel_titles = _addon_panel_titles(panel_titles, 1)
    panel_title = resolved_panel_titles === nothing ?
                  display_title : only(resolved_panel_titles)
    return with_theme(_addon_theme(export_theme = export_theme)) do
        shell = _addon_shell(; size, controls)
        groups = Dict{Symbol, Vector{Any}}()
        order = Symbol[]
        labels = Dict{Symbol, String}()
        axis, reset!,
        panel = _addon_preview_axis!(
            shell,
            (1, 1),
            panel_title,
            polygons,
            references,
            limits,
            groups,
            order,
            labels;
            earth_model,
            display_surface_gradient
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
            series_attributes,
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
