# One current placement for each native guide. Construction and mutations enter
# the same composition operation; detached guides retain their native edits.
function _addon_guide_position(position;legend=false,main=false)
    position===nothing && return nothing
    legend && position===:inside && return position
    main && position===Val(:content) && return position
    position in (:left,:right,:top,:bottom) && return position
    position isa Tuple && length(position)==2 &&
        all(v -> v isa Integer && !(v isa Bool) && v>0,position) && position!=(2,2) && return Int.(position)
    throw(ArgumentError("guide position must be a side or positive side-grid slot; (2,2) is reserved for content"))
end

function _addon_guide_gap(value)
    values=value isa Real ? ntuple(_ -> value,4) : value
    values isa Tuple && length(values)==4 && all(x -> x isa Real && isfinite(x) && x>=0,values) ||
        throw(ArgumentError("guide_gap must be nonnegative or a nonnegative (left,right,bottom,top) tuple"))
    return Float64.(values)
end

function _addon_guide_attributes(attributes)
    attributes isa NamedTuple || throw(ArgumentError("guide attributes must be a NamedTuple"))
    if haskey(attributes,:margin)
        _addon_guide_gap(to_value(attributes.margin))
    end
    for (key,symbols) in ((:halign,(:left,:center,:right)),(:valign,(:bottom,:center,:top)))
        haskey(attributes,key) || continue
        value=to_value(getproperty(attributes,key))
        value in symbols || value isa Real && 0<=value<=1 ||
            throw(ArgumentError("$key must be a native alignment symbol or a fraction in [0,1]"))
    end
    return attributes
end

function _addon_legend_owner(p,panel)
    panel===nothing && return (body=p.addon_state.shell.body,groups=p.addon_state.groups,
        order=p.addon_state.order,labels=p.addon_state.labels,bounds=p.addon_state.inside_bbox)
    panel=_addon_panel_identity(panel)
    haskey(p.addon_state.panel_data,panel) || throw(ArgumentError("panel $(repr(panel)) is absent from this figure"))
    data=p.addon_state.panel_data[panel]
    return (body=something(data.panel.layout,p.addon_state.shell.body),groups=data.groups,order=data.order,labels=data.labels,bounds=data.axis.scene.viewport)
end

function _addon_guide_state(kind,scope,position,attributes;overflow=:show_all,title=nothing)
    _addon_guide_position(position;legend=kind===:legend,main=kind===:colorbars)
    _addon_guide_attributes(attributes)
    return (;kind,scope,position=Ref{Any}(position),attributes=Ref{Any}(attributes),
        overflow=Ref(overflow),title=Ref{Any}(title),placed=Ref{Any}(_omitted),object=Ref{Any}(nothing),
        layout=Ref{Any}(nothing),wrapping=Ref{Any}(nothing),subscriptions=Any[],hidden=Ref(false),visibility=IdDict{Any,Bool}())
end

function _addon_native_guide_attributes!(object,attributes)
    for (key,value) in pairs(attributes)
        key in propertynames(typeof(object)) || throw(ArgumentError("unsupported native guide attribute $key"))
        getproperty(object,key)[]=value
    end
    return object
end

function _addon_hide_guide!(p,state)
    state.hidden[] && return nothing
    empty!(state.visibility)
    if state.layout[]!==nothing
        for object in p.figure.content
            if _addon_belongs_to_slot(object,state.layout[])
                state.visibility[object]=object.blockscene.visible[]
                object.blockscene.visible[]=false
            end
        end
    end
    state.layout[]===nothing || _addon_detach!(state.layout[])
    objects=state.kind===:legend ? (state.object[],) : something(state.object[],())
    for object in objects
        object===nothing && continue
        get!(state.visibility,object,object.blockscene.visible[])
        object.blockscene.visible[]=false
    end
    state.hidden[]=true
    return nothing
end

function _addon_restore_guide!(state)
    state.hidden[] || return nothing
    for (object,visible) in state.visibility
        object.blockscene.visible[]=visible
    end
    empty!(state.visibility)
    state.hidden[]=false
    return nothing
end

function _addon_guide_extent(state,dimension)
    object=state.layout[]
    object===nothing && return 0.0
    dimensions=object.layoutobservables.reporteddimensions[]
    value=dimensions.inner[dimension]
    protrusions=dimensions.outer
    extra=dimension==1 ? protrusions.left+protrusions.right : protrusions.bottom+protrusions.top
    return Float64(something(value,object.layoutobservables.computedbbox[].widths[dimension])+extra)
end

function _addon_alignment(value,dimension)
    value isa Union{GridLayoutBase.HorizontalAlignment,GridLayoutBase.VerticalAlignment} && return Float64(value.x)
    value isa Real && return Float64(value)
    value in (:left,:bottom) && return 0.0
    value in (:right,:top) && return 1.0
    return 0.5
end

function _addon_compose_guides!(p)
    state=p.addon_state
    state.composing_guides[] && return p
    state.composing_guides[]=true
    views=[axis.targetlimits[] for axis in p.axes]
    try
        for guide in values(state.guides)
            guide.layout[]===nothing || _addon_detach!(guide.layout[])
            if guide.kind===:legend && guide.object[]!==nothing
                _addon_detach!(guide.object[])
            end
        end
        for dock in state.guide_docks
            _addon_delete_subtree!(dock)
        end
        empty!(state.guide_docks)
        bodies=Any[state.shell.body]
        # Panels start with empty, zero-sized docks. Only a panel that owns a
        # guide needs recomposition, including when that guide is now hidden.
        # Resetting every other panel repeatedly propagates full native layout
        # and text measurements through large matrix figures.
        for guide in values(state.guides)
            guide.scope===nothing && continue
            body=state.panel_data[guide.scope].panel.layout
            body===nothing || any(existing -> existing===body,bodies) || push!(bodies,body)
        end
        for body in bodies
            # Keep the reserved 3×3 content scaffold; remove abandoned extension tracks.
            for row in size(body)[1]:-1:4
                GridLayoutBase.deleterow!(body,row)
            end
            for column in size(body)[2]:-1:4
                GridLayoutBase.deletecol!(body,column)
            end
            for i in (1,3)
                rowsize!(body,i,Fixed(0));colsize!(body,i,Fixed(0))
            end
            rowgap!(body,0);colgap!(body,0)
        end
        placements=Dict{Tuple{Any,Any},Vector{Any}}()
        for key in state.guide_order
            guide=state.guides[key]
            position=guide.position[]
            guide.kind===:colorbars && isempty(state.color_scales) && continue
            if position===nothing
                _addon_hide_guide!(p,guide)
                continue
            end
            owner=guide.kind===:legend ? _addon_legend_owner(p,guide.scope) :
                (body=state.shell.body,bounds=state.inside_bbox)
            guide.kind===:legend && !any(group -> haskey(owner.labels,group),owner.order) && continue
            if position===:inside
                _addon_place_legend!(p,guide,owner,nothing)
                continue
            end
            position===Val(:content) && begin
                _addon_place_colorbars!(p,guide,state.shell.canvas[1,1],:horizontal)
                continue
            end
            push!(get!(placements,(owner.body,position),Any[]),guide)
        end
        for ((body,position),guides) in placements
            slot,orientation=_addon_legend_slot(body,position)
            rows,columns=orientation===:vertical ? (2length(guides)+1,1) : (1,2length(guides)+1)
            dock=GridLayout(rows,columns;tellwidth=orientation===:vertical,tellheight=orientation===:horizontal)
            slot[]=dock
            push!(state.guide_docks,dock)
            row,column=_addon_dock_indices(position)
            gap=state.guide_gap
            if row!=2
                rowsize!(body,row,Auto(true));rowgap!(body,row<2 ? row : row-1,gap[row<2 ? 4 : 3])
            end
            if column!=2
                colsize!(body,column,Auto(true));colgap!(body,column<2 ? column : column-1,gap[column<2 ? 1 : 2])
            end
            for (index,guide) in enumerate(guides)
                target=orientation===:vertical ? dock[2index,1] : dock[1,2index]
                if guide.kind===:legend
                    _addon_place_legend!(p,guide,_addon_legend_owner(p,guide.scope),target;orientation)
                else
                    _addon_place_colorbars!(p,guide,target,orientation)
                end
            end
            dimension=orientation===:vertical ? 2 : 1
            bounds=body===state.shell.body ? state.inside_bbox[] :
                state.panel_data[first(guides).scope].axis.scene.viewport[]
            span=Float64(bounds.widths[dimension])
            extents=[max(0.,_addon_guide_extent(guide,dimension)) for guide in guides]
            anchors=map(guides) do guide
                object=guide.kind===:legend ? guide.object[] : guide.layout[]
                _addon_alignment(to_value(getproperty(object,dimension==1 ? :halign : :valign)),dimension)
            end
            # Fixed native tracks retain complete guide extents. Equal anchors
            # form one compact stack; different anchors pack in declaration order.
            fractions=dimension==1 ? anchors : 1 .- anchors
            total=sum(extents)
            span=max(span,total)
            getproperty(dock,dimension==1 ? :width : :height)[]=span
            assigned=dock.layoutobservables.suggestedbbox[]
            free=Float64(assigned.widths[dimension])-span
            alignment=abs(free)<1 ? 0.5 :
                clamp((Float64(bounds.origin[dimension])-Float64(assigned.origin[dimension]))/free,0.,1.)
            getproperty(dock,dimension==1 ? :halign : :valign)[]=alignment
            starts=Float64[]
            cursor=0.0
            for index in eachindex(guides)
                guide=guides[index]
                local_bounds=guide.scope===nothing ? bounds : state.panel_data[guide.scope].axis.scene.viewport[]
                local_span=Float64(local_bounds.widths[dimension])
                local_start=dimension==1 ? Float64(local_bounds.origin[1]-bounds.origin[1]) :
                    Float64(bounds.origin[2]+bounds.widths[2]-local_bounds.origin[2]-local_bounds.widths[2])
                desired=local_bounds!=bounds ? local_start+fractions[index]*max(0.,local_span-extents[index]) :
                    all(==(first(anchors)),anchors) ?
                    first(fractions)*max(0.,span-total)+sum(extents[1:index-1]) :
                    fractions[index]*max(0.,span-extents[index])
                start=clamp(desired,cursor,max(cursor,span-sum(extents[index:end])))
                push!(starts,start);cursor=start+extents[index]
            end
            previous=0.0
            for index in eachindex(guides)
                spacer=max(0.,starts[index]-previous)
                if orientation===:vertical
                    rowsize!(dock,2index-1,Fixed(spacer));rowsize!(dock,2index,Fixed(extents[index]))
                else
                    colsize!(dock,2index-1,Fixed(spacer));colsize!(dock,2index,Fixed(extents[index]))
                end
                previous=starts[index]+extents[index]
            end
            # Account for both ends of the native span. Leaving its final
            # track implicit would let GridLayout redistribute free space and
            # apply the guide's alignment a second time.
            if orientation===:vertical
                rowsize!(dock,2length(guides)+1,Fixed(max(0.,span-previous)))
            else
                colsize!(dock,2length(guides)+1,Fixed(max(0.,span-previous)))
            end
            rowgap!(dock,0);colgap!(dock,0)
        end
        for body in bodies
            occupied_rows=Set{Int}();occupied_columns=Set{Int}()
            for ((owner,position),_) in placements
                owner===body || continue
                row,column=_addon_dock_indices(position)
                push!(occupied_rows,row);push!(occupied_columns,column)
            end
            for row in 1:size(body)[1]
                row==2 || row in occupied_rows || rowsize!(body,row,Fixed(0))
            end
            for column in 1:size(body)[2]
                column==2 || column in occupied_columns || colsize!(body,column,Fixed(0))
            end
        end
        for guide in values(state.guides)
            if guide.kind===:legend && guide.position[]===:inside && guide.object[]!==nothing && !haskey(guide.attributes[],:bbox)
                guide.object[].layoutobservables.suggestedbbox[]=_addon_legend_owner(p,guide.scope).bounds[]
            end
        end
    finally
        for (axis,view) in zip(p.axes,views)
            axis.targetlimits[]==view || (axis.targetlimits[]=view)
        end
        state.composing_guides[]=false
    end
    return p
end

function _addon_watch_guide!(p,guide,object)
    attributes=object isa Colorbar ? (:ticks,:tickformat,:ticklabelsize,:ticklabelrotation,
        :ticklabelfont,:label,:labelsize,:labelfont,:vertical,:width,:height) :
        (:labelsize,:labelfont,:titlesize,:titlefont,:padding,:margin,:patchsize,:rowgap,:colgap,:orientation,:nbanks)
    append!(guide.subscriptions,_addon_watch_presentation!(p,object,attributes))
    push!(guide.subscriptions,on(object.blockscene,object.layoutobservables.autosize) do _
        _addon_compose_guides!(p)
        nothing
    end)
    return object
end

function _addon_place_legend!(p,guide,owner,target;orientation=:vertical)
    position=guide.position[]
    if guide.object[]===nothing
        if position===:inside
            guide.object[]=_addon_legend!(p.figure,owner.body,owner.groups,owner.order,owner.labels;
                dependent_plots=p.addon_state.dependent_plots,position=:inside,
                attributes=guide.attributes[],overflow=guide.overflow[],title=guide.title[],inside_bbox=owner.bounds[])
        else
            guide.layout[]=GridLayout()
            target[]=guide.layout[]
            guide.object[]=_addon_legend!(p.figure,owner.body,owner.groups,owner.order,owner.labels;
                dependent_plots=p.addon_state.dependent_plots,position,
                attributes=guide.attributes[],overflow=guide.overflow[],title=guide.title[],
                inside_bbox=owner.bounds,target=guide.layout[][1,1],target_orientation=orientation)
        end
        guide.object[]===nothing && return nothing
        _addon_watch_guide!(p,guide,guide.object[])
        if guide.overflow[]===:show_all
            guide.wrapping[]=_addon_grid_legend!(owner.bounds,guide.object[],guide.position;
                fitting_geometry=p.addon_state.fitting_geometry,automatic=!haskey(guide.attributes[],:orientation) && !haskey(guide.attributes[],:nbanks))
        end
    end
    legend=guide.object[]
    _addon_restore_guide!(guide)
    if position===:inside
        _addon_detach!(legend)
        if !haskey(guide.attributes[],:bbox)
            legend.layoutobservables.suggestedbbox[]=owner.bounds[]
        end
        if guide.placed[]!==position
            legend.tellwidth[]=false;legend.tellheight[]=false
        end
    else
        guide.layout[]===nothing && (guide.layout[]=GridLayout())
        target[]=guide.layout[]
        guide.layout[][1,1]=legend
        if guide.placed[]!==position
            guide.layout[].halign[]=position===:right ? :left : position===:left ? :right : :center
            guide.layout[].valign[]=position===:top ? :bottom : position===:bottom ? :top : :center
            legend.tellwidth[]=get(guide.attributes[],:tellwidth,true)
            legend.tellheight[]=get(guide.attributes[],:tellheight,true)
        end
    end
    guide.placed[]=position
    guide.wrapping[]===nothing || guide.wrapping[]()
    if guide.scope===nothing
        p.legend=legend
        p.addon_state.controls_enabled && (p.controls[:legend]=legend)
    else
        p.panel_legends[guide.scope]=legend
    end
    return legend
end

function _addon_place_colorbars!(p,guide,target,orientation)
    if guide.object[]===nothing
        attributes=guide.attributes[]
        if guide.position[]===Val(:content)
            # Standalone scales occupy the main canvas. Side docks instead
            # need a finite natural scale length for complete measurement.
            length_attribute=get(attributes,:vertical,false) ? (;height=Auto()) : (;width=Auto())
            attributes=merge(length_attribute,attributes)
        end
        result=_addon_colorbars!(target,p.addon_state.color_scales;
            attributes,orientation)
        guide.object[]=result.colorbars
        guide.layout[]=result.layout
        for bar in result.colorbars
            _addon_watch_guide!(p,guide,bar)
        end
    else
        target[]=guide.layout[]
    end
    guide.layout[]===nothing && return nothing
    _addon_restore_guide!(guide)
    _addon_native_guide_group!(guide.layout[],p.addon_state.colorbar_group_attributes[])
    p.colorbars=Any[guide.object[]...]
    return guide.object[]
end

function _addon_native_guide_group!(grid,attributes)
    _addon_guide_attributes(attributes)
    for (key,value) in pairs(attributes)
        if key===:margin
            grid.alignmode[]=Outside(_addon_guide_gap(to_value(value))...)
            continue
        end
        key in (:halign,:valign,:width,:height,:tellwidth,:tellheight) ||
            throw(ArgumentError("unsupported colorbar group attribute $key"))
        getproperty(grid,key)[]=value
    end
    return grid
end

function _addon_update_legend!(p,panel;position=_omitted,title=_omitted,overflow=_omitted,legend_labels=nothing,kwargs...)
    owner=_addon_legend_owner(p,panel)
    key=(:legend,panel)
    guide=get(p.addon_state.guides,key,nothing)
    resolved=position===_omitted ? (guide===nothing ? :right : guide.position[]) : position
    _addon_guide_position(resolved;legend=true)
    mode=overflow===_omitted ? (guide===nothing ? :show_all : guide.overflow[]) : overflow
    mode in (:show_all,:ellipsis) || throw(ArgumentError("legend overflow must be :show_all or :ellipsis"))
    attributes=_addon_guide_attributes((;kwargs...))
    haskey(attributes,:anchor) && throw(ArgumentError("anchor was removed; use native halign and valign"))
    _addon_relabel_legend!(owner.labels,owner.groups,owner.order,legend_labels)
    if guide===nothing
        guide=_addon_guide_state(:legend,panel,resolved,attributes;overflow=mode,title=title===_omitted ? nothing : title)
        p.addon_state.guides[key]=guide
        push!(p.addon_state.guide_order,key)
    elseif guide.object[]!==nothing && (mode!=guide.overflow[] || legend_labels!==nothing || haskey(attributes,:bbox))
        current_title,current_entries=only(guide.object[].entrygroups[])
        guide.title[]=current_title
        if legend_labels===nothing
            for (group,entry) in zip(owner.order,current_entries)
                entry.label[]=="(...)" && break
                owner.labels[group]=entry.label[]
            end
        end
        native=(; (name=>to_value(getproperty(guide.object[],name)) for name in propertynames(Legend) if name!==:entrygroups)...)
        guide.attributes[]=merge(guide.attributes[],native,attributes)
        foreach(off,guide.subscriptions);empty!(guide.subscriptions)
        _addon_remove_legend!(guide.object[])
        guide.layout[]===nothing || _addon_delete_subtree!(guide.layout[])
        guide.object[]=nothing;guide.layout[]=nothing;guide.wrapping[]=nothing
        empty!(guide.visibility);guide.hidden[]=false
    elseif guide.object[]!==nothing
        _addon_native_guide_attributes!(guide.object[],attributes)
    end
    guide.attributes[]=merge(guide.attributes[],attributes)
    guide.position[]=resolved;guide.overflow[]=mode
    if title!==_omitted
        guide.title[]=title
        if guide.object[]!==nothing
            entries=last(only(guide.object[].entrygroups[]))
            guide.object[].entrygroups[]=[(title,entries)]
        end
    end
    _addon_compose_guides!(p)
    if resolved===nothing
        panel===nothing ? (p.legend=nothing;delete!(p.controls,:legend)) : delete!(p.panel_legends,panel)
    end
    return panel===nothing ? p.legend : get(p.panel_legends,panel,nothing)
end

LineCableModels.figurelegend!(p::LineCableModels.UIPlot;kwargs...)=
    _addon_edit_presentation!(() -> _addon_update_legend!(p,nothing;kwargs...),p)
LineCableModels.panellegend!(p::LineCableModels.UIPlot,panel;kwargs...)=
    _addon_edit_presentation!(() -> _addon_update_legend!(p,panel;kwargs...),p)

function LineCableModels.figurecolorbars!(p::LineCableModels.UIPlot;position=_omitted,group_attributes=_omitted,kwargs...)
    return _addon_edit_presentation!(p) do
    guide=p.addon_state.guides[(:colorbars,nothing)]
    resolved=position===_omitted ? guide.position[] : position
    _addon_guide_position(resolved;main=position===_omitted)
    attributes=(;kwargs...)
    group_attributes===_omitted || (p.addon_state.colorbar_group_attributes[]=_addon_guide_attributes(group_attributes))
    guide.attributes[]=merge(guide.attributes[],attributes)
    guide.object[]===nothing || foreach(bar -> _addon_native_guide_attributes!(bar,attributes),guide.object[])
    guide.position[]=resolved
    _addon_compose_guides!(p)
    resolved===nothing && empty!(p.colorbars)
    return p.colorbars
    end
end
