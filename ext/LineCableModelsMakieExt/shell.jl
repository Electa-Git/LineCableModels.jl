function draw! end

_shell_kind(::Type{<:PlotBuilder.AbstractPlotDefinition}) = Val(:standard)

_page_legend(page::PlotPage) = page.payload.legend
function _page_colorbars(page::PlotPage)
    return hasproperty(page.payload, :colorbars) ? page.payload.colorbars : ()
end
_page_export(page::PlotPage) = page.payload.export_definition

function build_context(
        page::PlotPage;
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
        Observable("Ready."),
        Any[],
        Dict{Symbol, Any}(),
        nothing,
        Any[],
        nothing,
        nothing,
        Any[],
        Ref{Any}(nothing)
    )
end

function _standard_shell!(context::UIContext, page::PlotPage)
    figure = Figure(
        size = page.size,
        figure_padding = _window_padding((20, 20, 28, 28))
    )
    root = figure.layout
    root.default_rowgap = Fixed(GRID_ROW_GAP)
    root.default_colgap = Fixed(12)

    canvas = GridLayout(
        width = Auto(),
        height = Auto(),
        tellwidth = false,
        tellheight = false,
        halign = :center,
        valign = :center
    )
    canvas.default_rowgap = Fixed(GRID_ROW_GAP)
    canvas.default_colgap = Fixed(GRID_COLUMN_GAP)
    root[2, 1] = canvas

    side = GridLayout(2, 1)
    side.default_rowgap = Fixed(4)
    side.default_colgap = Fixed(0)
    rowsize!(side, 1, Relative(1))
    rowsize!(side, 2, Auto(true))
    colsize!(side, 1, Auto(true))
    root[2, 2] = side

    context.figure = figure
    context.canvas = canvas
    context.shell = (;
        kind = :standard,
        root,
        side,
        toolbar = root[1, 1:2],
        legend = side[1, 1],
        colorbars = side[2, 1],
        status = root[3, 1:2],
        collapsed = Set{Symbol}()
    )
    _apply_shell_tracks!(context)
    return context
end

function _colorbar_shell!(context::UIContext, page::PlotPage)
    figure = Figure(
        size = page.size,
        figure_padding = _window_padding((20, 20, 28, 28))
    )
    root = figure.layout
    root.default_rowgap = Fixed(GRID_ROW_GAP)
    root.default_colgap = Fixed(0)

    side = GridLayout(1, 1)
    side.default_rowgap = Fixed(0)
    side.default_colgap = Fixed(0)
    rowsize!(side, 1, Relative(1))
    colsize!(side, 1, Auto(true))
    root[2, 1] = side

    context.figure = figure
    context.canvas = nothing
    context.shell = (;
        kind = :colorbar,
        root,
        side,
        toolbar = root[1, 1],
        legend = nothing,
        colorbars = side[1, 1],
        status = root[3, 1],
        collapsed = Set{Symbol}()
    )
    _apply_shell_tracks!(context)
    return context
end

function build_shell(
        context::UIContext,
        definition::Type{<:PlotBuilder.AbstractPlotDefinition},
        page::PlotPage
)
    return build_shell(context, _shell_kind(definition), page)
end

function build_shell(context::UIContext, ::Val{:standard}, page::PlotPage)
    _standard_shell!(context, page)
end
function build_shell(context::UIContext, ::Val{:colorbar}, page::PlotPage)
    _colorbar_shell!(context, page)
end

function _set_rowsize!(grid, index::Int, size)
    index <= length(grid.rowsizes) && rowsize!(grid, index, size)
    return nothing
end

function _set_colsize!(grid, index::Int, size)
    index <= length(grid.colsizes) && colsize!(grid, index, size)
    return nothing
end

function _apply_shell_tracks!(context::UIContext)
    shell = context.shell
    collapsed = shell.collapsed
    _set_rowsize!(shell.root, 1, :toolbar in collapsed ? Fixed(0) : Fixed(36))
    _set_rowsize!(shell.root, 2, Relative(1))
    _set_rowsize!(shell.root, 3, :status in collapsed ? Fixed(0) : Fixed(20))
    if shell.kind === :standard
        _set_colsize!(shell.root, 1, Auto(false, 1))
        _set_colsize!(shell.root, 2, Auto(true))
        _set_rowsize!(shell.side, 1, :legend in collapsed ? Fixed(0) : Relative(1))
        _set_rowsize!(shell.side, 2, :colorbars in collapsed ? Fixed(0) : Auto(true))
        _set_colsize!(shell.side, 1, Auto(true))
    else
        _set_colsize!(shell.root, 1, Auto(true))
        _set_rowsize!(shell.side, 1, :colorbars in collapsed ? Fixed(0) : Relative(1))
        _set_colsize!(shell.side, 1, Auto(true))
    end
    return context
end

function draw!(context::UIContext, recipe::PlotRecipe, page::PlotPage)
    return draw!(context, recipe.definition, page)
end

format_axes!(context::UIContext, ::PlotPage) = context

function _collapse_toolbar!(context::UIContext)
    push!(context.shell.collapsed, :toolbar)
    _apply_shell_tracks!(context)
    return nothing
end

function _collapse_status!(context::UIContext)
    push!(context.shell.collapsed, :status)
    _apply_shell_tracks!(context)
    return nothing
end

function _collapse_legend!(context::UIContext)
    context.shell.kind === :standard || return nothing
    push!(context.shell.collapsed, :legend)
    _apply_shell_tracks!(context)
    return nothing
end

function _collapse_colorbars!(context::UIContext)
    push!(context.shell.collapsed, :colorbars)
    _apply_shell_tracks!(context)
    return nothing
end

function _collapse_empty_dock!(context::UIContext, page::PlotPage)
    context.shell.kind === :standard || return nothing
    context.legend === nothing && isempty(_page_colorbars(page)) || return nothing
    context.shell.side.width[] = 0
    context.shell.side.tellwidth[] = true
    return nothing
end

function place_legend!(
        context::UIContext,
        page::PlotPage;
        export_mode::Bool
)
    definition = _page_legend(page)
    if definition.enabled && context.shell.legend !== nothing
        overflow = export_mode ? :show_all : definition.overflow
        dock_width = isempty(_page_colorbars(page)) ? nothing : LEGEND_DOCK_WIDTH
        responsive = overflow === :ellipsis
        context.legend_slot_grid = GridLayout(
            width = dock_width === nothing ? Auto() : dock_width,
            height = responsive ? nothing : Auto(),
            tellwidth = true,
            tellheight = !responsive,
            halign = :left,
            valign = :top
        )
        context.shell.legend[] = context.legend_slot_grid
        built = _legend!(
            context.legend_slot_grid[1, 1],
            context.panels;
            width = dock_width,
            overflow
        )
        if built === nothing
            _collapse_legend!(context)
        else
            context.legend = built.legend
            context.responsive_legend = built.responsive
        end
    else
        _collapse_legend!(context)
    end
    return context
end

function place_colorbars!(context::UIContext, page::PlotPage)
    definitions = _page_colorbars(page)
    if isempty(definitions)
        _collapse_colorbars!(context)
    else
        alignment = context.shell.kind === :colorbar ? (:left, :center) : (:left, :top)
        built = _colorbars!(
            context.shell.colorbars,
            definitions;
            halign = first(alignment),
            valign = last(alignment)
        )
        append!(context.colorbars, built.colorbars)
    end
    return context
end

function _toolbar!(context::UIContext)
    toolbar = GridLayout(
        width = Auto(),
        height = Auto(),
        tellwidth = true,
        tellheight = true,
        halign = :left,
        valign = :bottom
    )
    toolbar.default_colgap = Fixed(4)
    context.shell.toolbar[] = toolbar
    return toolbar
end

function build_widgets!(
        context::UIContext,
        page::PlotPage;
        controls::Bool
)
    controls || begin
        _collapse_toolbar!(context)
        _collapse_status!(context)
        return context
    end

    panels = context.panels
    widgets = context.widgets
    xlog_available = _page_supports_log(panels, :x)
    ylog_available = _page_supports_log(panels, :y)
    reset_available = !isempty(panels)
    toolbar = _toolbar!(context)
    column = 1
    if reset_available
        reset = Button(
            toolbar[1, column];
            label = _icon_label(MI_REFRESH),
            width = BUTTON_SIZE,
            height = BUTTON_SIZE,
            buttoncolor = BUTTON_BACKGROUND
        )
        column += 1
        widgets[:reset] = reset
        push!(context.observers, on(reset.clicks) do _
            foreach(_reset_panel_limits!, panels)
            context.status[] = "Axis limits reset"
        end)
    end

    save_button = Button(
        toolbar[1, column];
        label = _icon_label(MI_SAVE),
        width = BUTTON_SIZE,
        height = BUTTON_SIZE,
        buttoncolor = BUTTON_BACKGROUND
    )
    column += 1
    widgets[:export_svg] = save_button
    push!(context.observers, on(save_button.clicks) do _
        try
            PlotBuilder.export_svg(context.plot_reference[])
        catch error
            context.status[] = sprint(showerror, error)
        end
    end)

    if xlog_available
        active = all(panel -> panel.axis.xscale[] === Makie.log10, panels)
        xlog = Toggle(toolbar[1, column], active = active)
        column += 1
        Label(toolbar[1, column], "log x")
        column += 1
        widgets[:xlog] = xlog
        push!(context.observers,
            on(xlog.active) do enabled
                scale = enabled ? :log10 : :linear
                foreach(panels) do panel
                    metadata = panel.metadata.xaxis
                    _set_axis_scale!(panel.axis, metadata, :x, metadata.exponent, scale)
                end
                foreach(_reset_panel_limits!, panels)
                context.status[] = enabled ?
                                   "x-axis scale set to log" :
                                   "x-axis scale set to linear"
            end)
    end

    if ylog_available
        active = all(panel -> panel.axis.yscale[] === Makie.log10, panels)
        ylog = Toggle(toolbar[1, column], active = active)
        column += 1
        Label(toolbar[1, column], "log y")
        widgets[:ylog] = ylog
        push!(context.observers,
            on(ylog.active) do enabled
                scale = enabled ? :log10 : :linear
                foreach(panels) do panel
                    metadata = panel.metadata.yaxis
                    _set_axis_scale!(panel.axis, metadata, :y, metadata.exponent, scale)
                end
                foreach(_reset_panel_limits!, panels)
                context.status[] = enabled ?
                                   "y-axis scale set to log" :
                                   "y-axis scale set to linear"
            end)
    end

    definition = _page_legend(page)
    definition.interactive && context.legend !== nothing &&
        (widgets[:legend] = context.legend)
    Label(
        context.shell.status,
        context.status;
        halign = :left,
        fontsize = 11
    )
    return context
end

function assemble(
        context::UIContext,
        recipe::PlotRecipe,
        page::PlotPage
)
    _apply_shell_tracks!(context)
    _collapse_empty_dock!(context, page)
    if context.responsive_legend !== nothing
        _observe_legend!(
            context.responsive_legend,
            context.legend_slot_grid,
            context
        )
        Makie.update_state_before_display!(context.figure)
        box = context.legend_slot_grid.layoutobservables.computedbbox[]
        _fit_legend!(context.responsive_legend, box.widths[2])
    end
    definition = _page_legend(page)
    definition.interactive && context.legend !== nothing &&
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
        page::PlotPage;
        backend,
        display::Bool,
        controls::Bool,
        export_mode::Bool
)
    context = build_context(page; backend, display)
    build_shell(context, recipe.definition, page)
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
    built = UIPlot[]
    for page in recipe.pages
        page_export_theme = export_theme === nothing ?
                            _page_export(page).theme : export_theme
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
