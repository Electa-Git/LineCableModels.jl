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
    for (role,latex_role) in ((:regular,:regular),(:italic,:italic),(:bold,:bold),(:bold_italic,:bolditalic))
        _addon_set_observable!(snapshot,figure_fonts[role],latex_fonts[latex_role][])
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


function _addon_hide_chrome!(snapshot,p)
    shell=p.addon_state.shell
    root=shell.root
    saved=(rows=copy(root.rowsizes),gaps=copy(root.addedrowgaps),offset=GridLayoutBase.firstrow(root))
    occupied=Set{Int}()
    for object in shell.chrome
        gc=GridLayoutBase.gridcontent(object)
        gc===nothing && continue
        _addon_hide_layout_content!(snapshot,object)
        union!(occupied,gc.span.rows)
    end
    for row in occupied
        rowsize!(root,row,Fixed(0))
        # Release only gaps adjacent to registered chrome, wherever it lives.
        row>GridLayoutBase.firstrow(root) && rowgap!(root,row-1,Fixed(0))
        row<GridLayoutBase.lastrow(root) && rowgap!(root,row,Fixed(0))
    end
    return saved
end

function _addon_export_presentation!(write,p,theme)
    state=p.addon_state
    state===nothing && throw(ArgumentError("SVG presentation requires a managed UIPlot"))
    size=Tuple(p.figure.scene.viewport[].widths)
    views=[axis.targetlimits[] for axis in p.axes]
    frames=isempty(p.axes) ? nothing : ntuple(i ->
        minimum(axis -> Float64(axis.layoutobservables.computedbbox[].widths[i]),p.axes),2)
    canvas=Tuple(state.shell.canvas.layoutobservables.computedbbox[].widths)
    canvas_size=(state.shell.canvas.width[],state.shell.canvas.height[])
    padding=copy(state.frame_padding)
    alignments=[axis.alignmode[] for axis in p.axes]
    snapshot=Pair{Any,Any}[]
    chrome=nothing
    previous_guard=state.fitting_geometry[]
    state.fitting_geometry[]=true
    try
        chrome=_addon_hide_chrome!(snapshot,p)
        _addon_publication_snapshot!(snapshot,p,theme)
        _addon_compose_guides!(p)
        if frames!==nothing && state.panel_page!==nothing
            _addon_panel_padding!(p)
            _addon_fit_frames!(p,frames)
        else
            body=state.shell.body.layoutobservables
            outside=Float64.(p.figure.scene.viewport[].widths-body.suggestedbbox[].widths+
                body.computedbbox[].widths-state.shell.canvas.layoutobservables.computedbbox[].widths)
            _addon_resize_preserving_views!(p,canvas.+outside)
        end
        return write()
    finally
        _addon_restore_snapshot!(snapshot)
        if chrome!==nothing
            for (row,value) in enumerate(chrome.rows)
                rowsize!(state.shell.root,row+chrome.offset-1,value)
            end
            for (row,value) in enumerate(chrome.gaps)
                rowgap!(state.shell.root,row+chrome.offset-1,value)
            end
        end
        for (axis,alignment,view) in zip(p.axes,alignments,views)
            axis.alignmode[]=alignment
            axis.targetlimits[]==view || (axis.targetlimits[]=view)
        end
        empty!(state.frame_padding);merge!(state.frame_padding,padding)
        state.shell.canvas.width[]=canvas_size[1]
        state.shell.canvas.height[]=canvas_size[2]
        try
            _addon_resize_preserving_views!(p,size)
        finally
            state.fitting_geometry[]=previous_guard
        end
    end
end
