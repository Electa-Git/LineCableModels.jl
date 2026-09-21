# One nominal capacity, independent of quantity and actual residual extent.
function _addon_capacity(layout, matrix_extents, flow_counts)
    if layout !== nothing
        layout isa Tuple && length(layout)==2 &&
            all(v -> v isa Integer && !(v isa Bool) && v>0,layout) ||
            throw(ArgumentError("layout must be a tuple of two positive integers or nothing"))
        return Int.(layout)
    end
    if !isempty(matrix_extents)
        return (maximum(first,matrix_extents),maximum(last,matrix_extents))
    end
    count=maximum(flow_counts;init=1)
    columns=ceil(Int,sqrt(count))
    return (cld(count,columns),columns)
end

function _addon_flow_pages(items, capacity)
    rows,columns=capacity
    pages=NamedTuple[]
    for start in 1:prod(capacity):length(items)
        selected=items[start:min(start+prod(capacity)-1,length(items))]
        positions=Tuple((cld(i,columns),mod1(i,columns)) for i in eachindex(selected))
        dimensions=(maximum(first,positions),maximum(last,positions))
        push!(pages,(;facets=selected,positions,dimensions,index=(cld(start,prod(capacity)),1)))
    end
    return pages
end


# Equalize actual native panel frames at construction. Native protrusions keep
# their scientific labels/ticks; only the space allocated around them changes.
# The original alignmode remains available for fresh local measurements after a
# title or guide edit. A caller's subsequent alignmode edit releases this padding.
function _addon_panel_padding!(p; minimum_decoration=(0.,0.,0.,0.))
    managed=p.addon_state.frame_padding
    for (axis,record) in managed
        axis.alignmode[]==record.applied && (axis.alignmode[]=record.original)
    end
    measurements=map(collect(values(p.addon_state.panel_data))) do data
        axis=data.axis
        data.panel.layout===nothing && return nothing
        panel=data.panel.layout.layoutobservables.computedbbox[]
        frame=axis.layoutobservables.computedbbox[]
        lower=frame.origin-panel.origin
        upper=panel.origin+panel.widths-frame.origin-frame.widths
        (;axis,decoration=(lower[1],upper[1],lower[2],upper[2]))
    end
    filter!(!isnothing,measurements)
    isempty(measurements) && return (0.,0.,0.,0.)
    maximum_decoration=ntuple(i -> max(minimum_decoration[i],maximum(m -> m.decoration[i],measurements)),4)
    for m in measurements
        axis=m.axis
        original=haskey(managed,axis) ? managed[axis].original : axis.alignmode[]
        # Respect constructor/native alignmode choices beyond ordinary Inside.
        axis.alignmode[] isa GridLayoutBase.Inside || continue
        native=axis.layoutobservables.protrusions[]
        values=ntuple(i -> Float32(getfield(native,i)+max(0.,maximum_decoration[i]-m.decoration[i])),4)
        applied=Mixed(; (key=>GridLayoutBase.Protrusion(value) for (key,value) in
            zip((:left,:right,:bottom,:top),values))...)
        axis.alignmode[]=applied
        managed[axis]=(;original,applied)
    end
    return maximum_decoration
end

function _addon_frame_measurements(p)
    state=p.addon_state
    canvas=state.shell.canvas
    viewport=p.figure.scene.viewport[]
    bounds=canvas.layoutobservables.computedbbox[]
    data=collect(values(state.panel_data))
    decorations=map(data) do panel
        box=panel.panel.layout.layoutobservables.computedbbox[]
        frame=panel.axis.layoutobservables.computedbbox[]
        Tuple(Float64.(box.widths-frame.widths))
    end
    decoration=ntuple(i -> maximum(x -> x[i],decorations;init=0.),2)
    gaps=(Float64(canvas.default_colgap.x),Float64(canvas.default_rowgap.x))
    body=state.shell.body.layoutobservables
    outer=Tuple(Float64.(viewport.widths-body.suggestedbbox[].widths+
        body.computedbbox[].widths-bounds.widths))
    return (;decoration,gaps,outer)
end

function _addon_resize_preserving_views!(p,size)
    views=[axis.targetlimits[] for axis in p.axes]
    try
        Makie.resize!(p.figure,Int.(ceil.(size))...)
        _addon_compose_guides!(p)
    finally
        for (axis,view) in zip(p.axes,views)
            axis.targetlimits[]==view || (axis.targetlimits[]=view)
        end
    end
    return p
end

function _addon_chrome_width(p)
    visible(object)=object isa Makie.Block ? object.blockscene.visible[] :
        any(block -> _addon_belongs_to_slot(block,object) && block.blockscene.visible[],p.figure.content)
    objects=Any[p.addon_state.shell.chrome...]
    p.title===nothing || push!(objects,p.title)
    required=maximum(objects;init=0.) do object
        visible(object) || return 0.
        Float64(something(object.layoutobservables.autosize[][1],0.))
    end
    padding=max(0.,Float64(p.figure.scene.viewport[].widths[1]-
        p.figure.layout.layoutobservables.computedbbox[].widths[1]))
    return required+padding
end

function _addon_fit_frames!(p,frames)
    state=p.addon_state
    haskey(state,:panel_page) || return p
    previous=state.fitting_geometry[]
    state.fitting_geometry[]=true
    try
        rows,columns=state.panel_page.dimensions
        measurements=_addon_frame_measurements(p)
        cell=frames .+ measurements.decoration
        state.shell.canvas.width[] isa Real &&
            (state.shell.canvas.width[]=columns*cell[1]+(columns-1)*measurements.gaps[1])
        state.shell.canvas.height[] isa Real &&
            (state.shell.canvas.height[]=rows*cell[2]+(rows-1)*measurements.gaps[2])
        width=measurements.outer[1]+columns*cell[1]+(columns-1)*measurements.gaps[1]
        guide_width=maximum(state.guide_docks;init=0.) do dock
            dock.width[] isa Real ? Float64(dock.width[]) : 0.
        end
        overflow_width=guide_width-Float64(state.inside_bbox[].widths[1])
        overflow_width=overflow_width>1 ? overflow_width : 0.
        minimum_width=max(_addon_chrome_width(p),width+overflow_width)
        if minimum_width>width
            # Long interactive content can enlarge the outer window while the
            # existing data frames remain at their current dimensions.
            state.shell.canvas.width[]=columns*cell[1]+(columns-1)*measurements.gaps[1]
            width=minimum_width
        end
        # Width determines measured legend wrapping before the height is fitted.
        _addon_resize_preserving_views!(p,(width,p.figure.scene.viewport[].widths[2]))
        measurements=_addon_frame_measurements(p)
        height=measurements.outer[2]+rows*cell[2]+(rows-1)*measurements.gaps[2]
        guide_height=maximum(state.guide_docks;init=0.) do dock
            dock.height[] isa Real ? Float64(dock.height[]) : 0.
        end
        overflow_height=guide_height-Float64(state.inside_bbox[].widths[2])
        overflow_height=overflow_height>1 ? overflow_height : 0.
        if overflow_height>0
            state.shell.canvas.height[]=rows*cell[2]+(rows-1)*measurements.gaps[2]
            height+=overflow_height
        end
        _addon_resize_preserving_views!(p,(width,height))
    finally
        state.fitting_geometry[]=previous
    end
    return p
end

function _addon_calibrate_frames!(pages,capacity)
    isempty(pages) && return pages
    measured=map(_addon_panel_padding!,pages)
    padding=ntuple(i -> maximum(x -> x[i],measured),4)
    for p in pages
        _addon_panel_padding!(p;minimum_decoration=padding)
    end
    br,bc=capacity
    candidates=map(pages) do p
        m=_addon_frame_measurements(p)
        nominal=p.figure.scene.viewport[].widths
        ((nominal[1]-m.outer[1]-(bc-1)*m.gaps[1])/bc-m.decoration[1],
         (nominal[2]-m.outer[2]-(br-1)*m.gaps[2])/br-m.decoration[2])
    end
    frames=ntuple(i -> minimum(x -> x[i],candidates),2)
    all(x -> isfinite(x) && x>1,frames) || throw(ArgumentError(
        "the nominal figure size cannot fit layout=$capacity and its measured decorations; increase the figure size or reduce layout"))
    for p in pages
        _addon_fit_frames!(p,frames)
    end
    return pages
end

# Physical aspect is a panel requirement. Fit native panel cells into the
# available canvas using each panel's data ratio, including wide systems and
# one-row collections. No figure or window aspect ratio is imposed.
function _addon_fit_panel_aspects!(p)
    isempty(p.axes) && return p
    all(axis -> axis.aspect[] isa DataAspect,p.axes) || return p
    state=p.addon_state
    shell=state.shell
    panels=collect(values(state.panel_data))
    all(data -> data.panel.layout!==nothing,panels) || return p
    rows,columns=size(shell.canvas)
    column_ratios=zeros(columns)
    column_decoration=zeros(columns)
    row_decoration=zeros(rows)
    for data in panels
        gc=GridLayoutBase.gridcontent(data.panel.layout)
        row=first(gc.span.rows);column=first(gc.span.cols)
        view=data.axis.targetlimits[]
        ratio=Float64(view.widths[1]/view.widths[2])
        isfinite(ratio) && ratio>0 || continue
        column_ratios[column]=max(column_ratios[column],ratio)
        panel=data.panel.layout.layoutobservables.computedbbox[]
        frame=data.axis.layoutobservables.computedbbox[]
        decoration=panel.widths-frame.widths
        column_decoration[column]=max(column_decoration[column],decoration[1])
        row_decoration[row]=max(row_decoration[row],decoration[2])
    end
    any(iszero,column_ratios) && return p
    body=shell.body.layoutobservables
    canvas=shell.canvas.layoutobservables.computedbbox[]
    outside=max.(0.,Float64.(body.computedbbox[].widths-canvas.widths))
    available=Float64.(body.suggestedbbox[].widths).-outside
    gaps=(Float64(shell.canvas.default_colgap.x),Float64(shell.canvas.default_rowgap.x))
    frame_height=min((available[1]-sum(column_decoration)-(columns-1)*gaps[1])/sum(column_ratios),
        (available[2]-sum(row_decoration)-(rows-1)*gaps[2])/rows)
    isfinite(frame_height) && frame_height>1 || return p # native zero-sized initialization
    column_widths=frame_height.*column_ratios.+column_decoration
    row_heights=frame_height.+row_decoration
    shell.canvas.width[]=sum(column_widths)+(columns-1)*gaps[1]
    shell.canvas.height[]=sum(row_heights)+(rows-1)*gaps[2]
    shell.canvas.tellwidth[]=true;shell.canvas.tellheight[]=true
    shell.body.width[]=Auto();shell.body.height[]=Auto()
    shell.body.tellwidth[]=true;shell.body.tellheight[]=true
    shell.body.halign[]=:center;shell.body.valign[]=:center
    colsize!(shell.body,2,Auto(true));rowsize!(shell.body,2,Auto(true))
    for column in 1:columns
        colsize!(shell.canvas,column,Fixed(column_widths[column]))
    end
    for row in 1:rows
        rowsize!(shell.canvas,row,Fixed(row_heights[row]))
    end
    return p
end

function _addon_frame_snapshot(p)
    frames=isempty(p.axes) ? nothing : ntuple(i ->
        minimum(axis -> Float64(axis.layoutobservables.computedbbox[].widths[i]),p.axes),2)
    canvas=Tuple(p.addon_state.shell.canvas.layoutobservables.computedbbox[].widths)
    views=[axis.targetlimits[] for axis in p.axes]
    return (;frames,canvas,views)
end

function _addon_edit_presentation!(action,p; before=nothing)
    state=p.addon_state
    state.presentation_ready[] || return action()
    state.fitting_geometry[] && return action()
    (;frames,canvas,views)=before===nothing ? _addon_frame_snapshot(p) : before
    state.fitting_geometry[]=true
    try
        result=action()
        _addon_compose_guides!(p)
        if frames!==nothing && state.panel_page!==nothing
            _addon_panel_padding!(p)
            _addon_fit_frames!(p,frames)
        else
            # Caller-owned nested canvases retain their topology and dimensions.
            # Only their surrounding chrome and guides change the outer window.
            body=state.shell.body.layoutobservables
            outside=Float64.(p.figure.scene.viewport[].widths-body.suggestedbbox[].widths+
                body.computedbbox[].widths-state.shell.canvas.layoutobservables.computedbbox[].widths)
            _addon_resize_preserving_views!(p,canvas.+outside)
        end
        return result
    finally
        for (axis,view) in zip(p.axes,views)
            axis.targetlimits[]==view || (axis.targetlimits[]=view)
        end
        state.fitting_geometry[]=false
    end
end

# Automatic flow pages may rearrange only their existing members on a manual
# window resize. Explicit layouts and matrix pages never install this callback.
function _addon_responsive_axis_grid!(p)
    state=p.addon_state
    identities=state.panel_page.coordinates
    length(identities)<=1 && return p
    panels=[state.panel_data[id] for id in identities]
    current_columns=Ref(state.panel_page.dimensions[2])
    on(p.figure.scene,p.figure.scene.viewport) do _
        p.addon_state.fitting_geometry[] && return nothing
        grid=state.shell.canvas
        box=state.shell.body.layoutobservables.suggestedbbox[]
        all(>(0),box.widths) || return nothing
        columns=clamp(round(Int,sqrt(length(panels)*box.widths[1]/box.widths[2])),1,length(panels))
        columns==current_columns[] && return nothing
        views=[data.axis.targetlimits[] for data in panels]
        state.fitting_geometry[]=true
        try
            # Empty matrix-selection cells are not part of flow membership.
            if haskey(state,:page_cells)
                used=Set(data.panel.layout for data in panels)
                for cell in values(state.page_cells)
                    cell.layout in used || _addon_delete_subtree!(cell.layout)
                end
                empty!(state.page_cells)
            end
            rows=cld(length(panels),columns)
            for (index,data) in enumerate(panels)
                row,column=cld(index,columns),mod1(index,columns)
                grid[row,column]=data.panel.layout
                for key in (:xlabelvisible,:xticklabelsvisible,:xticksvisible)
                    haskey(state.shell.axis_attributes,key) || (getproperty(data.axis,key)[]=row==rows)
                end
            end
            GridLayoutBase.trim!(grid)
            for row in 1:rows
                rowsize!(grid,row,Auto(false,1))
            end
            for column in 1:columns
                colsize!(grid,column,Auto(false,1))
            end
            current_columns[]=columns
            p.addon_state=merge(p.addon_state,(panel_page=merge(p.addon_state.panel_page,(dimensions=(rows,columns),)),))
            _addon_fit_panel_aspects!(p)
            _addon_compose_guides!(p)
        finally
            for (data,view) in zip(panels,views)
                data.axis.targetlimits[]==view || (data.axis.targetlimits[]=view)
            end
            state.fitting_geometry[]=false
        end
        nothing
    end
    return p
end

# Native attribute changes use the same local fitting operation as public guide
# mutators. Capture frames before native layout listeners run, then measure the
# completed native edit. Internal composition and resize notifications do not
# enter this path, and subscriptions belong to the actual edited block.
function _addon_watch_presentation!(p,block,names)
    subscriptions=Any[]
    for name in names
        name in propertynames(typeof(block)) || continue
        attribute=getproperty(block,name)
        attribute isa Observable || continue
        pending=Ref{Any}(nothing)
        push!(subscriptions,on(block.blockscene,attribute;priority=1) do _
            state=p.addon_state
            pending[]=state.presentation_ready[] && !state.fitting_geometry[] &&
                !state.composing_guides[] ? _addon_frame_snapshot(p) : nothing
            nothing
        end)
        push!(subscriptions,on(block.blockscene,attribute;priority=-100) do _
            before=pending[];pending[]=nothing
            before===nothing || _addon_edit_presentation!(() -> nothing,p;before)
            nothing
        end)
    end
    return subscriptions
end
