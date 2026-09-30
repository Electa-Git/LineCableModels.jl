# One current placement for each native guide. Construction and mutations enter
# the same composition operation. Detached guides retain their native edits.
function _guide_position(position;legend=false,main=false)
    position===nothing && return nothing
    legend && position===:inside && return position
    main && position===Val(:content) && return position
    position in (:left,:right,:top,:bottom) && return position
    position isa Tuple && length(position)==2 &&
        all(v -> v isa Integer && !(v isa Bool) && v>0,position) && position!=(2,2) && return Int.(position)
    throw(ArgumentError("guide position must be a side or positive side-grid slot; (2,2) is reserved for content"))
end

function _guide_gap(value)
    values=value isa Real ? ntuple(_ -> value,4) : value
    values isa Tuple && length(values)==4 && all(x -> x isa Real && isfinite(x) && x>=0,values) ||
        throw(ArgumentError("guide_gap must be nonnegative or a nonnegative (left,right,bottom,top) tuple"))
    return Float64.(values)
end

function _guide_spacing(value,current=(rowgap=12.0,colgap=12.0))
    fields=value isa Real ? (rowgap=value,colgap=value) : value
    fields isa NamedTuple && all(key -> key in (:rowgap,:colgap),keys(fields)) ||
        throw(ArgumentError("guide_spacing must be a nonnegative number or (rowgap=..., colgap=...)"))
    result=merge(current,fields)
    all(x -> x isa Real && !(x isa Bool) && isfinite(x) && x>=0,values(result)) ||
        throw(ArgumentError("guide spacing must be finite, nonnegative, and non-Boolean"))
    return map(Float64,result)
end

function _colorbar_group_attributes(attributes,count=nothing)
    _guide_attributes(attributes)
    for (key,value) in pairs(attributes)
        if key===:layout
            value===nothing && continue
            value isa Tuple && length(value)==2 &&
                all(x -> x isa Integer && !(x isa Bool) && 0<x<=typemax(Int),value) ||
                throw(ArgumentError("colorbar group layout must be nothing or two positive non-Boolean integers"))
            count===nothing || big(value[1])*value[2]>=count ||
                throw(ArgumentError("colorbar group layout $value cannot hold $count supplied scales"))
        elseif key in (:rowgap,:colgap)
            value===nothing || _guide_spacing((;key=>value))
        elseif key in (:width,:height)
            value isa Real && (!(value isa Bool) && isfinite(value) && value>0) && continue
            value===nothing || value isa Union{Auto,Relative} ||
                throw(ArgumentError("group $key must be a positive native size"))
        elseif key in (:tellwidth,:tellheight)
            value isa Bool || throw(ArgumentError("$key must be Boolean"))
        elseif key ∉ (:halign,:valign,:margin)
            throw(ArgumentError("unsupported colorbar group attribute $key"))
        end
    end
    return attributes
end

# Ordered one-dimensional packing. Explicit spacer tracks include each adjacency
# once. Native gaps on those tracks remain zero. The caller handles overflow.
function _pack_guides(extents,span,gap,fractions;desired=nothing)
    isempty(extents) && return (;starts=Float64[],extent=0.0,span=Float64(span))
    suffix=reverse(cumsum(reverse(extents)))
    total=first(suffix)+(length(extents)-1)*gap
    span=max(span,total)
    equal=desired===nothing && all(==(first(fractions)),fractions)
    starts=Float64[]
    cursor=0.0
    for i in eachindex(extents)
        target=desired!==nothing ? desired[i] : equal ?
            first(fractions)*(span-total)+(first(suffix)-suffix[i])+(i-1)*gap :
            fractions[i]*(span-extents[i])
        upper=span-suffix[i]-(length(extents)-i)*gap
        start=clamp(target,cursor,max(cursor,upper))
        push!(starts,start)
        cursor=start+extents[i]+gap
    end
    return (;starts,extent=total,span)
end

function _guide_attributes(attributes)
    attributes isa NamedTuple || throw(ArgumentError("guide attributes must be a NamedTuple"))
    if haskey(attributes,:margin)
        _guide_gap(to_value(attributes.margin))
    end
    for (key,symbols) in ((:halign,(:left,:center,:right)),(:valign,(:bottom,:center,:top)))
        haskey(attributes,key) || continue
        value=to_value(getproperty(attributes,key))
        value in symbols || value isa Real && 0<=value<=1 ||
            throw(ArgumentError("$key must be a native alignment symbol or a fraction in [0,1]"))
    end
    return attributes
end

function _legend_owner(p,panel)
    panel===nothing && return (body=p.plot_state.shell.body,groups=p.plot_state.groups,
        order=p.plot_state.order,labels=p.plot_state.labels,bounds=p.plot_state.inside_bbox)
    panel=_panel_identity(panel)
    haskey(p.plot_state.panel_data,panel) || throw(ArgumentError("panel $(repr(panel)) is absent from this figure"))
    data=p.plot_state.panel_data[panel]
    return (body=something(data.panel.layout,p.plot_state.shell.body),groups=data.groups,order=data.order,labels=data.labels,bounds=data.axis.scene.viewport)
end

function _guide_state(kind,scope,position,attributes;max_fraction=0.5,title=nothing)
    _guide_position(position;legend=kind===:legend,main=kind===:colorbars)
    _guide_attributes(attributes)
    return (;kind,scope,position=Ref{Any}(position),attributes=Ref{Any}(attributes),
        max_fraction=Ref(_legend_fraction(max_fraction)),title=Ref{Any}(title),placed=Ref{Any}(_omitted),object=Ref{Any}(nothing),
        layout=Ref{Any}(nothing),items=Any[],fit=Ref{Any}(nothing),subscriptions=Any[],hidden=Ref(false),visibility=IdDict{Any,Bool}())
end

function _native_guide_attributes!(object,attributes)
    previous=Pair{Any,Any}[]
    try
        for (key,value) in pairs(attributes)
            key in propertynames(typeof(object)) || throw(ArgumentError("unsupported native guide attribute $key"))
            attribute=getproperty(object,key)
            push!(previous,attribute=>attribute[])
            attribute[]=value
        end
    catch
        for (attribute,value) in reverse(previous)
            attribute[]=value
        end
        rethrow()
    end
    return object
end

function _validate_native_guide(type,attributes)
    _guide_attributes(attributes)
    for (key,attribute) in pairs(attributes)
        key in propertynames(type) || throw(ArgumentError("unsupported native guide attribute $key"))
        value=to_value(attribute)
        convert(fieldtype(type,key).parameters[1],value)
        if key in (:vertical,:labelvisible,:ticklabelsvisible,:tellwidth,:tellheight)
            value isa Bool || throw(ArgumentError("$key must be Boolean"))
        elseif key in (:width,:height)
            if value isa Real
                !(value isa Bool) && isfinite(value) && value>0 ||
                    throw(ArgumentError("$key must be a finite positive native dimension"))
            else
                value===nothing || value isa Union{Auto,Relative} || throw(ArgumentError("invalid native $key"))
            end
        elseif key in (:labelsize,:ticklabelsize,:titlesize,:spinewidth)
            value isa Real && !(value isa Bool) && isfinite(value) && value>=0 ||
                throw(ArgumentError("$key must be finite and nonnegative"))
        elseif type===Legend && key===:orientation
            value in (:horizontal,:vertical) || throw(ArgumentError("legend orientation must be :horizontal or :vertical"))
        elseif type===Legend && key===:nbanks
            value isa Integer && !(value isa Bool) && value>0 || throw(ArgumentError("legend nbanks must be a positive integer"))
        elseif endswith(string(key),"color")
            Makie.to_color(value)
        elseif key===:colormap
            Makie.to_colormap(value)
        end
    end
    return attributes
end

function _hide_guide!(p,state)
    state.hidden[] && return nothing
    empty!(state.visibility)
    if state.layout[]!==nothing
        for object in p.figure.content
            if _belongs_to_slot(object,state.layout[])
                state.visibility[object]=object.blockscene.visible[]
                object.blockscene.visible[]=false
            end
        end
    end
    state.layout[]===nothing || _detach!(state.layout[])
    objects=state.kind===:legend ? (state.object[],) : something(state.object[],())
    for object in objects
        object===nothing && continue
        get!(state.visibility,object,object.blockscene.visible[])
        object.blockscene.visible[]=false
    end
    state.hidden[]=true
    return nothing
end

function _restore_guide!(state)
    state.hidden[] || return nothing
    for (object,visible) in state.visibility
        object.blockscene.visible[]=visible
    end
    empty!(state.visibility)
    state.hidden[]=false
    return nothing
end

function _guide_extent(state,dimension)
    object=state.layout[]
    object===nothing && return 0.0
    dimensions=object.layoutobservables.reporteddimensions[]
    value=dimensions.inner[dimension]
    protrusions=dimensions.outer
    extra=dimension==1 ? protrusions.left+protrusions.right : protrusions.bottom+protrusions.top
    return Float64(something(value,object.layoutobservables.computedbbox[].widths[dimension])+extra)
end

function _alignment(value,dimension)
    hasproperty(value,:x) && getproperty(value,:x) isa Real && return Float64(value.x)
    value isa Real && return Float64(value)
    value in (:left,:bottom) && return 0.0
    value in (:right,:top) && return 1.0
    return 0.5
end

function _compose_guides!(p)
    state=p.plot_state
    state.composing_guides[] && return p
    state.composing_guides[]=true
    views=[axis.targetlimits[] for axis in p.axes]
    try
        for guide in values(state.guides)
            guide.layout[]===nothing || _detach!(guide.layout[])
            if guide.kind===:legend && guide.object[]!==nothing
                _detach!(guide.object[])
            end
        end
        for dock in state.guide_docks
            _delete_subtree!(dock)
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
            # Keep the reserved 3×3 content scaffold. Remove abandoned extension tracks.
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
                _hide_guide!(p,guide)
                continue
            end
            _restore_guide!(guide)
            guide.kind===:colorbars && guide.object[]!==nothing && (p.colorbars=Any[guide.object[]...])
            guide.kind===:legend && guide.object[]!==nothing && !guide.object[].blockscene.visible[] && continue
            if guide.kind===:colorbars && !isempty(guide.items) && !any(item -> item.visible[],guide.items)
                _restore_guide!(guide)
                _reflow_colorbars!(p,guide)
                continue
            end
            owner=guide.kind===:legend ? _legend_owner(p,guide.scope) :
                (body=state.shell.body,bounds=state.inside_bbox)
            guide.kind===:legend && !any(group -> haskey(owner.labels,group),owner.order) && continue
            if position===:inside
                _place_legend!(p,guide,owner,nothing)
                continue
            end
            position===Val(:content) && begin
                _place_colorbars!(p,guide,state.shell.canvas[1,1],:horizontal)
                continue
            end
            push!(get!(placements,(owner.body,position),Any[]),guide)
        end
        for ((body,position),guides) in placements
            slot,orientation=_legend_slot(body,position)
            rows,columns=orientation===:vertical ? (2length(guides)+1,1) : (1,2length(guides)+1)
            # Both complete extents participate in content fitting. Flexible
            # construction and resize tracks still resolve the reference allocation.
            dock=GridLayout(rows,columns;tellwidth=true,tellheight=true)
            slot[]=dock
            push!(state.guide_docks,dock)
            row,column=_dock_indices(position)
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
                    _place_legend!(p,guide,_legend_owner(p,guide.scope),target;orientation)
                else
                    _place_colorbars!(p,guide,target,orientation)
                end
            end
            dimension=orientation===:vertical ? 2 : 1
            bounds=body===state.shell.body ? state.inside_bbox[] :
                state.panel_data[first(guides).scope].axis.scene.viewport[]
            span=Float64(bounds.widths[dimension])
            extents=[max(0.,_guide_extent(guide,dimension)) for guide in guides]
            anchors=map(guides) do guide
                object=guide.kind===:legend ? guide.object[] : guide.layout[]
                _alignment(to_value(getproperty(object,dimension==1 ? :halign : :valign)),dimension)
            end
            # Fixed native tracks retain complete guide extents. Equal anchors
            # form one compact stack. Different anchors pack in declaration order.
            fractions=dimension==1 ? anchors : 1 .- anchors
            gap=dimension==1 ? state.guide_spacing[].colgap : state.guide_spacing[].rowgap
            desired=if any(guide -> guide.scope!==nothing && state.panel_data[guide.scope].axis.scene.viewport[]!=bounds,guides)
                map(enumerate(guides)) do (index,guide)
                    local_bounds=guide.scope===nothing ? bounds : state.panel_data[guide.scope].axis.scene.viewport[]
                    local_span=Float64(local_bounds.widths[dimension])
                    local_start=dimension==1 ? Float64(local_bounds.origin[1]-bounds.origin[1]) :
                        Float64(bounds.origin[2]+bounds.widths[2]-local_bounds.origin[2]-local_bounds.widths[2])
                    local_start+fractions[index]*max(0.,local_span-extents[index])
                end
            else
                nothing
            end
            packed=_pack_guides(extents,span,gap,fractions;desired)
            span=packed.span
            getproperty(dock,dimension==1 ? :width : :height)[]=span
            assigned=dock.layoutobservables.suggestedbbox[]
            free=Float64(assigned.widths[dimension])-span
            # A strip wider than its data frame still follows the requested
            # alignment, including the negative free space on either side.
            anchor=all(==(first(anchors)),anchors) ? first(anchors) : 0.5
            origin=Float64(bounds.origin[dimension])-anchor*max(0.,span-Float64(bounds.widths[dimension]))
            alignment=abs(free)<1 ? 0.5 :
                clamp((origin-Float64(assigned.origin[dimension]))/free,0.,1.)
            getproperty(dock,dimension==1 ? :halign : :valign)[]=alignment
            starts=packed.starts
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
                row,column=_dock_indices(position)
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
                guide.object[].layoutobservables.suggestedbbox[]=_legend_owner(p,guide.scope).bounds[]
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

function _watch_guide!(p,guide,object)
    attributes=object isa Colorbar ? (:ticks,:tickformat,:ticklabelsize,:ticklabelrotation,
        :ticklabelfont,:label,:labelsize,:labelfont,:labelvisible,:ticklabelsvisible,
        :labelpadding,:ticklabelpad,:spinewidth,:flipaxis,:vertical,:width,:height) :
        (:labelsize,:labelfont,:titlesize,:titlefont,:padding,:margin,:patchsize,:rowgap,:colgap,:orientation,:nbanks)
    append!(guide.subscriptions,_watch_presentation!(p,object,attributes))
    push!(guide.subscriptions,on(object.blockscene,object.layoutobservables.autosize) do _
        _compose_guides!(p)
        nothing
    end)
    if object isa Legend
        push!(guide.subscriptions,on(object.blockscene,object.blockscene.visible) do _
            state=p.plot_state
            state.composing_guides[] || state.fitting_geometry[] ||
                _update_figure_layout!(() -> nothing,p)
            nothing
        end)
    end
    return object
end

function _place_legend!(p,guide,owner,target;orientation=:vertical)
    position=guide.position[]
    if guide.object[]===nothing
        if position!==:inside
            guide.layout[]=GridLayout()
            target[]=guide.layout[]
        end
        guide.object[],guide.fit[]=_legend!(p.figure,owner.groups,owner.order,owner.labels;
            dependent_plots=p.plot_state.dependent_plots,position=guide.position,
            attributes=guide.attributes[],max_fraction=guide.max_fraction,title=guide.title[],
            inside_bbox=owner.bounds,target=position===:inside ? nothing : guide.layout[][1,1],
            target_orientation=orientation,fitting_geometry=p.plot_state.fitting_geometry)
        _watch_guide!(p,guide,guide.object[])
    end
    legend=guide.object[]
    _restore_guide!(guide)
    if position===:inside
        _detach!(legend)
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
    guide.fit[]===nothing || guide.fit[]()
    if guide.scope===nothing
        p.legend=legend
        p.plot_state.controls_enabled && (p.controls[:legend]=legend)
    else
        p.panel_legends[guide.scope]=legend
    end
    return legend
end

function _watch_scale_item!(p,guide,item)
    _watch_guide!(p,guide,item.bar)
    for name in (:vertical,:width,:height,:labelvisible)
        push!(guide.subscriptions,on(item.bar.blockscene,getproperty(item.bar,name);priority=2) do value
            state=p.plot_state
            state.composing_guides[] || state.fitting_geometry[] || begin
                if name===:labelvisible
                    item.labelvisible[]=value
                else
                    getproperty(item.managed,name)[]=false
                end
            end
            nothing
        end)
    end
    push!(guide.subscriptions,on(item.bar.blockscene,item.bar.blockscene.visible) do visible
        state=p.plot_state
        state.composing_guides[] || state.fitting_geometry[] || visible==item.visible[] || begin
            _update_figure_layout!(p) do
                item.visible[]=visible
            end
        end
        nothing
    end)
    return item
end

function _reflow_colorbars!(p,guide)
    grid=guide.layout[]
    grid===nothing && return nothing
    attributes=p.plot_state.colorbar_group_attributes[]
    _colorbar_group_attributes(attributes,length(p.plot_state.color_scales))
    position=guide.position[]
    active=filter(item -> item.visible[],guide.items)
    requested=get(attributes,:layout,nothing)
    columns=requested===nothing ? (position in (:top,:bottom) ? max(1,length(active)) : 1) : Int(requested[2])
    rowgap=something(get(attributes,:rowgap,nothing),p.plot_state.guide_spacing[].rowgap)
    colgap=something(get(attributes,:colgap,nothing),p.plot_state.guide_spacing[].colgap)
    for item in guide.items
        bar=item.bar
        if item.managed.vertical[]
            vertical=position ∉ (:top,:bottom,Val(:content))
            bar.vertical[]==vertical || (bar.vertical[]=vertical)
        end
        natural=position===Val(:content) ? Auto() : _COLORBAR_DOCK_LENGTH
        # Auto keeps Makie's native `size` as the bar thickness. `nothing`
        # instead stretches that dimension across the assigned group track.
        for (key,value) in pairs(bar.vertical[] ? (width=Auto(),height=natural) : (width=natural,height=Auto()))
            item.managed[key][] && !isequal(getproperty(bar,key)[],value) && (getproperty(bar,key)[]=value)
        end
        side_label=!bar.vertical[]
        bar.labelvisible[]==(!side_label && item.labelvisible[]) || (bar.labelvisible[]=!side_label && item.labelvisible[])
        item.companion.blockscene.visible[]=item.visible[] && side_label && item.labelvisible[]
        _detach!(item.layout)
        _detach!(bar)
        _detach!(item.companion)
        if item.visible[]
            if side_label && item.labelvisible[]
                item.layout[1,1]=item.companion
                item.layout[1,2]=bar
            else
                item.layout[1,1]=bar
            end
        end
        with_updates_suspended(item.layout) do
            GridLayoutBase.trim!(item.layout)
            rowgap!(item.layout,0);colgap!(item.layout,8)
            colsize!(item.layout,1,Auto())
            # Outside grids consume raw child protrusions. Include only the extra
            # measured border and text offset. Folding all text into the bar's inner
            # size would count native perpendicular decorations twice here.
            raw=bar.layoutobservables.protrusions[]
            complete=bar.layoutobservables.reporteddimensions[].outer
            padding=ntuple(i -> max(0.,getfield(complete,i)-getfield(raw,i)),4)
            item.layout.alignmode[]=Outside(padding...)
        end
        # A detached GridContent retains its parentless placement object.
        # Settle its current native sizes before measuring local bar offsets.
        item.layout.layoutobservables.suggestedbbox[]=item.layout.layoutobservables.suggestedbbox[]
    end
    # Label and endpoint tracks are shared down each group column. Keep the
    # native bar dimensions. Only their surrounding decoration space expands.
    for column in 1:min(columns,length(active))
        members=active[column:columns:end]
        compact=filter(item -> item.companion.blockscene.visible[],members)
        isempty(compact) && continue
        labelwidth=maximum(item.companion.layoutobservables.autosize[][1] for item in compact)
        leading=[item.bar.layoutobservables.computedbbox[].origin[1]-
            item.bar.layoutobservables.suggestedbbox[].origin[1] for item in compact]
        endpoint=maximum(leading)
        for (item,before) in zip(compact,leading)
            colsize!(item.layout,1,Fixed(labelwidth))
            colgap!(item.layout,8+endpoint-before)
            item.layout.layoutobservables.suggestedbbox[]=item.layout.layoutobservables.suggestedbbox[]
        end
    end
    # Complete item boxes are for collision clearance, not alignment anchors.
    # Reserve the same measured insets about each bar in a column or row. This
    # keeps equal bar edges and baselines aligned despite unequal native decoration.
    insets=map(active) do item
        frame=item.bar.layoutobservables.computedbbox[]
        outer=item.layout.layoutobservables.computedbbox[]
        before=frame.origin-outer.origin
        after=outer.origin+outer.widths-frame.origin-frame.widths
        (before[1],after[1],before[2],after[2])
    end
    for (index,item) in enumerate(active)
        row,column=cld(index,columns),mod1(index,columns)
        shared=ntuple(4) do side
            maximum(insets[j][side] for j in eachindex(active)
                if side<=2 ? mod1(j,columns)==column : cld(j,columns)==row)
        end
        padding=item.layout.alignmode[].padding
        with_updates_suspended(item.layout) do
            item.layout.alignmode[]=Outside(ntuple(side ->
                getfield(padding,side)+max(0.,shared[side]-insets[index][side]),4)...)
        end
        grid[row,column]=item.layout
    end
    GridLayoutBase.trim!(grid)
    grid.default_rowgap=Fixed(rowgap);grid.default_colgap=Fixed(colgap)
    rowgap!(grid,rowgap);colgap!(grid,colgap)
    for (dimension,key) in enumerate((:width,:height))
        limit=to_value(getproperty(grid,key))
        required=grid.layoutobservables.autosize[][dimension]
        limit isa Real && required!==nothing && limit+1<required && throw(ArgumentError(
            "colorbar group $key=$limit cannot enclose its complete items (requires at least $required logical pixels); enlarge it or change group layout"))
    end
    return grid
end

function _place_colorbars!(p,guide,target,orientation)
    if guide.object[]===nothing
        _colorbar_group_attributes(p.plot_state.colorbar_group_attributes[],length(p.plot_state.color_scales))
        result=_colorbars!(target,p.plot_state.color_scales;
            attributes=guide.attributes[],orientation,main=guide.position[]===Val(:content),
            scene=p.figure.scene)
        guide.object[]=result.colorbars
        guide.layout[]=result.layout
        append!(guide.items,result.items)
        _native_guide_group!(guide.layout[],p.plot_state.colorbar_group_attributes[])
        for item in guide.items
            _watch_scale_item!(p,guide,item)
        end
    end
    _restore_guide!(guide)
    _reflow_colorbars!(p,guide)
    target[]=guide.layout[]
    p.colorbars=Any[guide.object[]...]
    return guide.object[]
end

function _native_guide_group!(grid,attributes)
    for (key,value) in pairs(attributes)
        key in (:layout,:rowgap,:colgap) && continue
        if key===:margin
            grid.alignmode[]=Outside(_guide_gap(to_value(value))...)
            continue
        end
        key in (:halign,:valign,:width,:height,:tellwidth,:tellheight) ||
            throw(ArgumentError("unsupported colorbar group attribute $key"))
        getproperty(grid,key)[]=value
    end
    return grid
end

function _update_legend!(p,panel;position=_omitted,title=_omitted,max_fraction=_omitted,legend_labels=nothing,kwargs...)
    owner=_legend_owner(p,panel)
    key=(:legend,panel)
    guide=get(p.plot_state.guides,key,nothing)
    resolved=position===_omitted ? (guide===nothing ? :right : guide.position[]) : position
    _guide_position(resolved;legend=true)
    fraction=_legend_fraction(max_fraction===_omitted ? (guide===nothing ? 0.5 : guide.max_fraction[]) : max_fraction)
    attributes=_validate_native_guide(Legend,(;kwargs...))
    haskey(attributes,:anchor) && throw(ArgumentError("anchor was removed; use native halign and valign"))
    previous_guard=p.plot_state.composing_guides[]
    p.plot_state.composing_guides[]=true
    try
        _relabel_legend!(owner.labels,owner.groups,owner.order,legend_labels)
        if guide===nothing
            guide=_guide_state(:legend,panel,resolved,attributes;max_fraction=fraction,title=title===_omitted ? nothing : title)
            p.plot_state.guides[key]=guide
            push!(p.plot_state.guide_order,key)
        elseif guide.object[]!==nothing && (legend_labels!==nothing || haskey(attributes,:bbox))
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
            _remove_legend!(guide.object[])
            guide.layout[]===nothing || _delete_subtree!(guide.layout[])
            guide.object[]=nothing;guide.layout[]=nothing;guide.fit[]=nothing
            empty!(guide.visibility);guide.hidden[]=false
        elseif guide.object[]!==nothing
            _native_guide_attributes!(guide.object[],attributes)
        end
        guide.attributes[]=merge(guide.attributes[],attributes)
        guide.position[]=resolved;guide.max_fraction[]=fraction
        if title!==_omitted
            guide.title[]=title
            if guide.object[]!==nothing
                entries=last(only(guide.object[].entrygroups[]))
                _legend_entries!(guide.object[],title,entries)
            end
        end
    finally
        p.plot_state.composing_guides[]=previous_guard
    end
    guide.fit[]===nothing || guide.fit[]()
    _compose_guides!(p)
    if resolved===nothing
        panel===nothing ? (p.legend=nothing;delete!(p.controls,:legend)) : delete!(p.panel_legends,panel)
    end
    return panel===nothing ? p.legend : get(p.panel_legends,panel,nothing)
end

function LineCableModels.figurelegend!(p::LineCableModels.UIPlot;guide_spacing=_omitted,kwargs...)
    state=p.plot_state
    spacing=guide_spacing===_omitted ? state.guide_spacing[] : _guide_spacing(guide_spacing,state.guide_spacing[])
    # Validate native options and placement before changing the shared setting.
    native=(; (key=>value for (key,value) in kwargs if key ∉ (:position,:title,:max_fraction,:legend_labels))...)
    _validate_native_guide(Legend,native)
    haskey(kwargs,:position) && _guide_position(kwargs[:position];legend=true)
    haskey(kwargs,:max_fraction) && _legend_fraction(kwargs[:max_fraction])
    previous=state.guide_spacing[]
    try
        return _update_figure_layout!(p) do
            state.guide_spacing[]=spacing
            _update_legend!(p,nothing;kwargs...)
        end
    catch
        state.guide_spacing[]=previous
        rethrow()
    end
end
LineCableModels.panellegend!(p::LineCableModels.UIPlot,panel;kwargs...)=
    _update_figure_layout!(() -> _update_legend!(p,panel;kwargs...),p)

function LineCableModels.figurecolorbars!(p::LineCableModels.UIPlot;position=_omitted,
        group_attributes=_omitted,guide_spacing=_omitted,kwargs...)
    state=p.plot_state
    guide=state.guides[(:colorbars,nothing)]
    resolved=position===_omitted ? guide.position[] : position
    _guide_position(resolved;main=position===_omitted)
    group_attributes===_omitted || group_attributes isa NamedTuple ||
        throw(ArgumentError("group_attributes must be a NamedTuple"))
    group=_colorbar_group_attributes(group_attributes===_omitted ? state.colorbar_group_attributes[] :
        merge(state.colorbar_group_attributes[],group_attributes),length(state.color_scales))
    spacing=guide_spacing===_omitted ? state.guide_spacing[] : _guide_spacing(guide_spacing,state.guide_spacing[])
    attributes=_validate_native_guide(Colorbar,(;kwargs...))
    previous=(position=guide.position[],attributes=guide.attributes[],group=state.colorbar_group_attributes[],
        spacing=state.guide_spacing[],size=Tuple(p.figure.scene.viewport[].widths),
        native=guide.layout[]===nothing ? nothing : (; (key=>to_value(getproperty(guide.layout[],key)) for key in
            (:width,:height,:halign,:valign,:tellwidth,:tellheight,:alignmode))...))
    saved=map(guide.items) do item
        native=(; (key=>to_value(getproperty(item.bar,key)) for key in
            unique((keys(attributes)...,:vertical,:width,:height,:labelvisible)))...)
        (;item,native,managed=map(x -> x[],item.managed),labelvisible=item.labelvisible[])
    end
    try
        _update_figure_layout!(p) do
            state.composing_guides[]=true
            try
                state.colorbar_group_attributes[]=group
                state.guide_spacing[]=spacing
                if guide.layout[]!==nothing && group_attributes!==_omitted
                    _native_guide_group!(guide.layout[],group_attributes)
                end
                guide.attributes[]=merge(guide.attributes[],attributes)
                for item in guide.items
                    for key in (:vertical,:width,:height)
                        haskey(attributes,key) && (item.managed[key][]=false)
                    end
                    haskey(attributes,:labelvisible) && (item.labelvisible[]=to_value(attributes.labelvisible))
                    _native_guide_attributes!(item.bar,attributes)
                end
                guide.position[]=resolved
            finally
                state.composing_guides[]=false
            end
        end
    catch
        state.composing_guides[]=true
        try
            guide.position[]=previous.position;guide.attributes[]=previous.attributes
            state.colorbar_group_attributes[]=previous.group;state.guide_spacing[]=previous.spacing
            if previous.native!==nothing
                for (key,value) in pairs(previous.native)
                    getproperty(guide.layout[],key)[]=value
                end
            elseif guide.layout[]!==nothing
                foreach(off,guide.subscriptions);empty!(guide.subscriptions)
                for item in guide.items
                    delete!(item.bar);delete!(item.companion)
                    _detach!(item.layout)
                end
                _detach!(guide.layout[])
                empty!(guide.items)
                guide.layout[]=nothing;guide.object[]=nothing
                empty!(p.colorbars)
            end
            for record in saved
                _native_guide_attributes!(record.item.bar,record.native)
                for key in keys(record.managed)
                    record.item.managed[key][]=record.managed[key]
                end
                record.item.labelvisible[]=record.labelvisible
            end
        finally
            state.composing_guides[]=false
        end
        _resize_preserving_views!(p,previous.size)
        rethrow()
    end
    resolved===nothing && empty!(p.colorbars)
    return p.colorbars
end
