function _addon_observable_snapshot!(snapshot, object, names)
    for name in names
        hasproperty(object, name) || continue
        value = getproperty(object, name)
        value isa Observable || continue
        push!(snapshot, value => value[])
    end
    return snapshot
end

function _addon_set_observable!(snapshot, observable, value)
    push!(snapshot, observable => observable[])
    observable[] = value
    return observable
end

function _addon_hide_layout_content!(snapshot, content)
    if content isa GridLayout
        for entry in content.content
            _addon_hide_layout_content!(snapshot, entry.content)
        end
    elseif hasproperty(content, :blockscene)
        _addon_set_observable!(snapshot, content.blockscene.visible, false)
    end
    return snapshot
end

function _addon_publication_snapshot!(snapshot, plot::LineCableModels.UIPlot, theme::Symbol)
    theme in (:default, :publication) || throw(ArgumentError(
        "theme must be :default or :publication",
    ))
    _addon_observable_snapshot!(snapshot, plot.figure.scene, (:backgroundcolor,))
    plot.figure.scene.backgroundcolor[] = Makie.to_color(:white)
    theme === :default && return nothing

    latex_fonts = Makie.theme_latexfonts().attributes[:fonts][]
    figure_fonts = plot.figure.scene.theme[:fonts]
    for (role, latex_role) in ((:regular, :regular), (:italic, :italic),
        (:bold, :bold), (:bold_italic, :bolditalic))
        _addon_set_observable!(snapshot, figure_fonts[role], latex_fonts[latex_role][])
    end
    regular = :regular
    bold = :bold
    for block in plot.figure.content
        if block isa Axis
            _addon_observable_snapshot!(snapshot, block,
                (
                    :titlefont,
                    :xlabelfont,
                    :ylabelfont,
                    :xticklabelfont,
                    :yticklabelfont
                ))
            block.titlefont[] = bold
            block.xlabelfont[] = regular
            block.ylabelfont[] = regular
            block.xticklabelfont[] = regular
            block.yticklabelfont[] = regular
        elseif block isa Legend
            _addon_observable_snapshot!(snapshot, block, (:labelfont, :titlefont))
            block.labelfont[] = regular
            block.titlefont[] = bold
        elseif block isa Colorbar
            _addon_observable_snapshot!(snapshot, block, (:labelfont, :ticklabelfont))
            block.labelfont[] = regular
            block.ticklabelfont[] = regular
        elseif block isa Label
            _addon_observable_snapshot!(snapshot, block, (:font,))
            block.font[] = block === plot.title ? bold : regular
        end
    end
    return nothing
end

function _addon_restore_snapshot!(snapshot)
    for (observable, value) in Iterators.reverse(snapshot)
        observable[] = value
    end
    return nothing
end

function _addon_hide_chrome!(snapshot, p)
    shell=p.addon_state.shell
    root=shell.root
    first_row=1+offsets(root)[1]
    last_row=nrows(root)+offsets(root)[1]
    saved=(
        rows = copy(root.rowsizes), gaps = copy(root.addedrowgaps), offset = first_row)
    occupied=Set{Int}()
    for object in shell.chrome
        gc=GridLayoutBase.gridcontent(object)
        gc===nothing && continue
        _addon_hide_layout_content!(snapshot, object)
        union!(occupied, gc.span.rows)
    end
    with_updates_suspended(root) do
        for row in occupied
            rowsize!(root, row, Fixed(0))
            # GridLayoutBase 0.11.2 and 0.11.3 assign different meanings to
            # indexed rowgap! calls for offset grids. The stored gaps retain
            # one stable order, so change only the entries adjacent to chrome.
            gap_after=row-first_row+1
            row>first_row && (root.addedrowgaps[gap_after-1]=Fixed(0))
            row<last_row && (root.addedrowgaps[gap_after]=Fixed(0))
        end
    end
    return saved
end

function _addon_export_presentation!(write, p, theme)
    state=p.addon_state
    state===nothing && throw(ArgumentError("SVG presentation requires a managed UIPlot"))
    size=Tuple(p.figure.scene.viewport[].widths)
    (; frames, canvas, views)=_addon_frame_snapshot(p)
    grids=Any[state.shell.root, state.shell.body, state.shell.canvas]
    append!(grids, [data.panel.layout
                    for data in values(state.panel_data) if data.panel.layout!==nothing])
    tracks=[(; grid, rows = copy(grid.rowsizes), columns = copy(grid.colsizes),
                row_offset = 1+offsets(grid)[1], column_offset = 1+offsets(grid)[2])
            for grid in grids]
    padding=copy(state.frame_padding)
    alignments=[axis.alignmode[] for axis in p.axes]
    snapshot=Pair{Any, Any}[]
    for grid in grids
        _addon_observable_snapshot!(snapshot, grid, (
            :width, :height, :tellwidth, :tellheight))
    end
    chrome=nothing
    previous_guard=state.fitting_geometry[]
    state.fitting_geometry[]=true
    try
        chrome=_addon_hide_chrome!(snapshot, p)
        _addon_publication_snapshot!(snapshot, p, theme)
        _addon_compose_guides!(p)
        all(axis -> axis.aspect[]===nothing, p.axes) && state.panel_page!==nothing &&
            _addon_panel_padding!(p)
        _addon_fit_frames!(p, frames; canvas)
        return write()
    finally
        _addon_restore_snapshot!(snapshot)
        if chrome!==nothing
            with_updates_suspended(state.shell.root) do
                for (row, value) in enumerate(chrome.rows)
                    rowsize!(state.shell.root, row+chrome.offset-1, value)
                end
                state.shell.root.addedrowgaps .= chrome.gaps
            end
        end
        for (axis, alignment, view) in zip(p.axes, alignments, views)
            axis.alignmode[]=alignment
            axis.targetlimits[]==view || (axis.targetlimits[]=view)
        end
        empty!(state.frame_padding)
        merge!(state.frame_padding, padding)
        for record in tracks
            for (i, value) in enumerate(record.rows)
                rowsize!(record.grid, i+record.row_offset-1, value)
            end
            for (i, value) in enumerate(record.columns)
                colsize!(record.grid, i+record.column_offset-1, value)
            end
        end
        try
            _addon_resize_preserving_views!(p, size)
        finally
            state.fitting_geometry[]=previous_guard
        end
    end
end
