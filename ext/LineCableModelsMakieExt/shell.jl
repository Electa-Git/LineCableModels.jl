function build_context(
        page::PageSpec;
        backend,
        display::Bool
)
    active = PlotBuilder.ensure_backend!(backend)
    interactive = display && active in (:gl, :wgl)
    window = interactive && active === :gl ?
             PlotBuilder.make_screen(
        "Fig. $(PlotBuilder.next_fignum()) – $(page.title)";
        backend = :gl
    ) : nothing
    return UIContext(
        active,
        interactive,
        window,
        nothing,
        nothing,
        nothing,
        Observable(page.status.initial),
        Any[],
        Dict{Symbol, Any}(),
        nothing,
        nothing,
        nothing,
        Any[],
        Ref{Any}(nothing)
    )
end

function build_shell(context::UIContext, page::PageSpec)
    root = only(filter(grid -> grid.parent === nothing, page.layout.grids))
    context.figure = Figure(
        size = page.size,
        figure_padding = _window_padding(root.padding)
    )
    context.materialized = _materialize_layout(context.figure, page.layout)
    if haskey(context.materialized.slot_specs, :canvas)
        context.canvas = _slot_grid(
            context.materialized,
            :canvas;
            tellwidth = false,
            tellheight = false
        )
        context.canvas.default_rowgap = Fixed(GRID_ROW_GAP)
        context.canvas.default_colgap = Fixed(GRID_COLUMN_GAP)
    end
    return context
end

function draw!(context::UIContext, recipe::PlotRecipe, page::PageSpec)
    return draw!(context, recipe.spec, recipe, page)
end

function draw!(
        context::UIContext,
        ::Type{D},
        ::PlotRecipe,
        page::PageSpec
) where {D <: PlotBuilder.AbstractPlotDefinition}
    context.panels = _build_panels(page, context.materialized)
    return context
end

function format_axes!(context::UIContext, ::PageSpec)
    return context
end

function place_legend!(
        context::UIContext,
        page::PageSpec;
        export_mode::Bool
)
    materialized = context.materialized
    if page.legend.enabled
        overflow = export_mode ? :show_all : page.legend.overflow
        dock_width = _legend_dock_width(page)
        responsive = overflow === :ellipsis
        context.legend_slot_grid = _slot_grid(
            materialized,
            page.legend.slot;
            width = dock_width === nothing ? Auto() : dock_width,
            height = responsive ? nothing : Auto(),
            tellheight = !responsive
        )
        built_legend = _legend!(
            context.legend_slot_grid[1, 1],
            context.panels,
            materialized.slot_specs[page.legend.slot];
            width = dock_width,
            overflow
        )
        if built_legend === nothing
            _collapse_slot!(page.layout, materialized, page.legend.slot)
        else
            context.legend = built_legend.legend
            context.responsive_legend = built_legend.responsive
        end
    else
        _collapse_slot!(page.layout, materialized, page.legend.slot)
    end
    return context
end

function place_colorbars!(context::UIContext, page::PageSpec)
    if isempty(page.colorbars)
        for slot in page.layout.slots
            slot.name === :colorbars &&
                _collapse_slot!(page.layout, context.materialized, slot.name)
        end
    else
        _build_colorbars!(page, context.materialized)
    end
    return context
end

function build_widgets!(
        context::UIContext,
        page::PageSpec;
        controls::Bool
)
    panels = context.panels
    widgets = context.widgets
    materialized = context.materialized
    if controls
        definitions = page.controls
        xlog_available = _page_supports_log(panels, :x)
        ylog_available = _page_supports_log(panels, :y)
        toolbar_enabled = definitions.reset || definitions.export_svg ||
                          xlog_available || ylog_available
        if toolbar_enabled
            toolbar = _slot_grid(materialized, definitions.slot)
            toolbar.default_colgap = Fixed(4)
            column = 1
            if definitions.reset
                reset = Button(
                    toolbar[1, column];
                    label = _icon_label(MI_REFRESH),
                    width = BUTTON_SIZE,
                    height = BUTTON_SIZE,
                    buttoncolor = BUTTON_BACKGROUND
                )
                column += 1
                widgets[:reset] = reset
                on(reset.clicks) do _
                    foreach(_reset_panel_limits!, panels)
                    context.status[] = "Axis limits reset"
                end
            end
            if definitions.export_svg
                save_button = Button(
                    toolbar[1, column];
                    label = _icon_label(MI_SAVE),
                    width = BUTTON_SIZE,
                    height = BUTTON_SIZE,
                    buttoncolor = BUTTON_BACKGROUND
                )
                column += 1
                widgets[:export_svg] = save_button
                on(save_button.clicks) do _
                    try
                        PlotBuilder.export_svg(context.plot_reference[])
                    catch error
                        context.status[] = sprint(showerror, error)
                    end
                end
            end
            if xlog_available
                active = !isempty(panels) && all(
                    panel -> panel.axis.xscale[] === Makie.log10,
                    panels
                )
                xlog = Toggle(toolbar[1, column], active = active)
                column += 1
                Label(toolbar[1, column], "log x")
                column += 1
                widgets[:xlog] = xlog
                on(xlog.active) do enabled
                    scale = enabled ? :log10 : :linear
                    foreach(
                        panel -> _set_axis_scale!(
                            panel.axis,
                            panel.view.xaxis,
                            :x,
                            panel.view.xaxis.exponent,
                            scale
                        ),
                        panels
                    )
                    foreach(_reset_panel_limits!, panels)
                    context.status[] = enabled ?
                                       "x-axis scale set to log" :
                                       "x-axis scale set to linear"
                end
            end
            if ylog_available
                active = !isempty(panels) && all(
                    panel -> panel.axis.yscale[] === Makie.log10,
                    panels
                )
                ylog = Toggle(toolbar[1, column], active = active)
                column += 1
                Label(toolbar[1, column], "log y")
                widgets[:ylog] = ylog
                on(ylog.active) do enabled
                    scale = enabled ? :log10 : :linear
                    foreach(
                        panel -> _set_axis_scale!(
                            panel.axis,
                            panel.view.yaxis,
                            :y,
                            panel.view.yaxis.exponent,
                            scale
                        ),
                        panels
                    )
                    foreach(_reset_panel_limits!, panels)
                    context.status[] = enabled ?
                                       "y-axis scale set to log" :
                                       "y-axis scale set to linear"
                end
            end
        else
            _collapse_slot!(page.layout, materialized, page.controls.slot)
        end
        page.legend.interactive && context.legend !== nothing &&
            (widgets[:legend] = context.legend)
        if page.status.enabled
            Label(
                _slot_position(materialized, page.status.slot),
                context.status;
                halign = :left,
                fontsize = 11
            )
        else
            _collapse_slot!(page.layout, materialized, page.status.slot)
        end
    else
        _collapse_slot!(page.layout, materialized, page.controls.slot)
        _collapse_slot!(page.layout, materialized, page.status.slot)
    end
    return context
end

function assemble(
        context::UIContext,
        recipe::PlotRecipe,
        page::PageSpec
)
    _collapse_empty_dock!(page, context.materialized, context.legend)
    _apply_layout_specs!(page.layout, context.materialized)
    if context.responsive_legend !== nothing
        _observe_legend!(
            context.responsive_legend,
            context.legend_slot_grid,
            context
        )
        Makie.update_state_before_display!(context.figure)
        bounding_box = context.legend_slot_grid.layoutobservables.computedbbox[]
        _fit_legend!(context.responsive_legend, bounding_box.widths[2])
    end
    page.legend.interactive && context.legend !== nothing &&
        _observe_visibility_limits!(context.panels, context)
    built = UIPlot(
        recipe,
        page,
        context.figure,
        context.panels,
        context.widgets,
        context
    )
    context.plot_reference[] = built
    return built
end

function display!(context::UIContext, plot::UIPlot, display_plot::Bool)
    display_plot || return plot
    if context.interactive && context.window !== nothing
        Base.display(context.window, plot.figure)
    else
        PlotBuilder.renderfig(plot.figure)
    end
    return plot
end

function build_page(
        recipe::PlotRecipe,
        page::PageSpec;
        backend,
        display::Bool,
        controls::Bool,
        export_mode::Bool
)
    context = build_context(page; backend, display)
    build_shell(context, page)
    draw!(context, recipe, page)
    format_axes!(context, page)
    place_legend!(context, page; export_mode)
    place_colorbars!(context, page)
    build_widgets!(context, page; controls)
    plot = assemble(context, recipe, page)
    return display!(context, plot, display)
end

function build(
        recipe::PlotRecipe;
        backend = nothing,
        display::Bool = true,
        controls::Bool = true,
        export_mode::Bool = false,
        export_theme::Union{Nothing, Symbol} = nothing
)
    PlotBuilder.validate(recipe)
    built = UIPlot[]
    for page in recipe.figures
        page_export_theme = export_theme === nothing ?
                            page.export_spec.theme : export_theme
        with_theme(_theme(; export_mode, export_theme = page_export_theme)) do
            push!(
                built,
                build_page(
                    recipe,
                    page;
                    backend,
                    display,
                    controls,
                    export_mode
                )
            )
        end
    end
    return built
end
