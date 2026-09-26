# One nominal capacity, independent of quantity and actual residual extent.
function _addon_capacity(layout, matrix_positions, flow_counts)
    if layout !== nothing
        layout isa Tuple && length(layout)==2 &&
        all(v -> v isa Integer && !(v isa Bool) && v>0, layout) ||
            throw(ArgumentError("layout must be a tuple of two positive integers or nothing"))
        return Int.(layout)
    end
    if !isempty(matrix_positions)
        spans=map(positions -> _addon_panel_footprint(positions).dimensions, matrix_positions)
        return (maximum(first, spans), maximum(last, spans))
    end
    count=maximum(flow_counts; init = 1)
    columns=ceil(Int, sqrt(count))
    return (cld(count, columns), columns)
end

function _addon_figure_size(fig_size, dimensions)
    if fig_size !== nothing
        fig_size isa Tuple{Int, Int} && all(>(0), fig_size) || throw(ArgumentError(
            "fig_size must be a tuple of two positive integers or nothing",
        ))
        return fig_size
    end
    rows, columns = dimensions
    return (max(680, 390columns + 180), max(440, 290rows + 100))
end

function _addon_panel_footprint(positions)
    origin=ntuple(d -> minimum(p -> p[d], positions), 2)
    dimensions=ntuple(d -> maximum(p -> p[d], positions)-origin[d]+1, 2)
    return (; origin, dimensions)
end

# The original block identity and the occupied rectangle are independent.
# Explicit capacities use origin (1,1); automatic overviews use the selected
# minimum coordinate. Empty perimeter tracks disappear, internal holes remain.
function _addon_matrix_pages(positions, extent, capacity; origin = (1, 1))
    all(p -> all(d -> 1<=p[d]<=extent[d], 1:2), positions) ||
        throw(ArgumentError("selected panel coordinates exceed the retained matrix extent"))
    blocks=map(p -> ntuple(d -> fld(p[d]-origin[d], capacity[d])+1, 2), positions)
    pages=NamedTuple[]
    for index in sort!(unique(blocks))
        members=findall(==(index), blocks)
        selected=positions[members]
        footprint=_addon_panel_footprint(selected)
        local_positions=Tuple(ntuple(d -> p[d]-footprint.origin[d]+1, 2) for p in selected)
        push!(pages, (; members, positions = local_positions, footprint..., index))
    end
    return pages
end

function _addon_flow_pages(items, capacity)
    rows, columns=capacity
    pages=NamedTuple[]
    for start in 1:prod(capacity):length(items)
        selected=items[start:min(start + prod(capacity) - 1, length(items))]
        positions=Tuple((cld(i, columns), mod1(i, columns)) for i in eachindex(selected))
        dimensions=(maximum(first, positions), maximum(last, positions))
        push!(pages, (; facets = selected, positions, dimensions,
            index = (cld(start, prod(capacity)), 1)))
    end
    return pages
end

# Assign retained coordinates and gridpoints to panels and curves, then paginate.
function _addon_observation_facets(published, ydata, orientations, retained, slots, reference_position, labels,
        coordinate_labels, coordinate_identities; explicit_labels=false)
    facets = NamedTuple[]
    for request_index in eachindex(ydata)
        product = first(published).observations[request_index]
        c = product.coordinates
        rows, columns, _ = first(published).coordinates[request_index]
        extent=c.kind===:matrix ?
               ntuple(
            d -> maximum(source.observations[request_index].coordinates.extent[d]
            for source in published),
            length(c.extent)) : c.extent
        orientation=orientations[request_index]
        if orientation===:rows
            for index in retained[request_index]
                source_product=published[index].observations[request_index]
                source_coordinates=source_product.coordinates
                row_labels=index==1 || explicit_labels ? coordinate_labels[request_index] :
                    _coordinate_labels(source_product,published[index].coordinates[request_index],orientation)
                point=slots[index]==0 ? reference_position : slots[index]
                for (local_column, column) in enumerate(columns)
                    curves=[(source_index=index, local_row, local_column, local_sample=nothing,
                        slot=row, role=:candidate,
                        group=Symbol("row_", row), label=row_labels[local_row],
                        position=local_row) for (local_row, row) in enumerate(rows)]
                    push!(facets, (;request_index, local_row=1, local_column,
                        row=point, column, row_label="",
                        column_label=get(source_coordinates,:column_labels,source_coordinates.labels)[column],
                        panel_identity=(point, column), curves, orientation, kind=c.kind, extent,
                        domain=get(source_coordinates,:domain,:unspecified),
                        column_domain=get(source_coordinates,:column_domain,nothing),
                        quantity=product.quantity, identity=request_identity(ydata[request_index]),
                        source_index=index, source_label=labels[request_index][index]))
                end
            end
            continue
        end
        coordinates=[(local_row=local_row,local_column=local_column,local_sample=nothing,row,
            column=c.kind===:diagonal ? row : column,
            row_label=hasproperty(c,:labels) ? c.labels[row] : string(row),
            column_label=haskey(c,:column_labels) ? c.column_labels[column] : string(column),
            style_slot=c.kind===:matrix ? (row-1)*c.extent[2]+column :
                c.kind in (:vector,:diagonal) ? row : ordinal,
            position=ordinal)
            for (ordinal,(local_row,local_column,row,column)) in enumerate(
                (local_row,local_column,row,column) for (local_row,row) in enumerate(rows)
                for (local_column,column) in enumerate(columns))]
        if orientation===:gridpoints
            for coordinate in coordinates
                curves=[(source_index=index,local_row=coordinate.local_row,
                    local_column=coordinate.local_column,local_sample=coordinate.local_sample,slot=slots[index],
                    role=slots[index]==0 ? :reference : :candidate,
                    group=Symbol("result_",index),label=labels[request_index][index])
                    for index in retained[request_index]]
                panel_identity=c.kind in (:matrix,:diagonal) ?
                    (coordinate.row,coordinate.column) :
                    c.kind===:vector ? coordinate.row : coordinate.position
                push!(facets,(;request_index,coordinate.local_row,coordinate.local_column,
                    coordinate.row,coordinate.column,coordinate.row_label,coordinate.column_label,
                    panel_identity,curves,orientation,kind=c.kind,extent,
                    domain=get(c,:domain,:unspecified),column_domain=get(c,:column_domain,nothing),
                    quantity=product.quantity,identity=request_identity(ydata[request_index])))
            end
        else
            for index in retained[request_index]
                source_coordinates=published[index].observations[request_index].coordinates
                panel_coordinates=c.kind===:assemblies ?
                    [begin
                        slot=findfirst(isequal(source_coordinates.labels[assembly]),
                            coordinate_identities[request_index])
                        (local_row=1,local_column=1,local_sample=ordinal,
                            style_slot=slot,position=slot)
                    end
                        for (ordinal,assembly) in enumerate(source_coordinates.assemblies)] : coordinates
                curves=[(source_index=index,local_row=coordinate.local_row,
                    local_column=coordinate.local_column,local_sample=coordinate.local_sample,
                    slot=coordinate.style_slot,
                    role=:candidate,group=Symbol("coordinate_",coordinate.style_slot),
                    label=coordinate_labels[request_index][coordinate.position],
                    position=coordinate.position)
                    for coordinate in panel_coordinates]
                push!(facets,(;request_index,local_row=1,local_column=1,
                    row=slots[index]==0 ? reference_position : slots[index],column=1,
                    row_label="",column_label="",panel_identity=slots[index]==0 ?
                        reference_position : slots[index],curves,orientation,
                    kind=c.kind,extent,domain=get(c,:domain,:unspecified),
                    column_domain=get(c,:column_domain,nothing),quantity=product.quantity,
                    identity=request_identity(ydata[request_index]),
                    source_label=labels[request_index][index]))
            end
        end
    end
    return facets
end

function _addon_observation_pages(facets, capacities; automatic = false)
    pages=NamedTuple[]
    for request_index in unique(facet.request_index for facet in facets)
        selected=filter(facet -> facet.request_index==request_index, facets)
        first_facet=first(selected)
        capacity=capacities[request_index]
        if first_facet.orientation===:rows
            for index in unique(f.source_index for f in selected)
                point_panels=filter(f -> f.source_index==index, selected)
                append!(pages, (merge(page,(;capacity))
                    for page in _addon_flow_pages(point_panels, capacity)))
            end
            continue
        end
        if first_facet.kind!==:matrix || first_facet.orientation===:coordinates
            append!(pages, (merge(page,(;capacity)) for page in _addon_flow_pages(selected, capacity)))
            continue
        end
        positions=[(f.row, f.column) for f in selected]
        origin=automatic ? _addon_panel_footprint(positions).origin : (1, 1)
        for page in _addon_matrix_pages(positions, first_facet.extent, capacity; origin)
            push!(pages,
                (; facets = selected[page.members], positions = page.positions,capacity,
                    dimensions = page.dimensions, origin = page.origin, index = page.index))
        end
    end
    return pages
end

# Equalize actual native panel frames at construction. Native protrusions keep
# their scientific labels/ticks; only the space allocated around them changes.
# The original alignmode remains available for fresh local measurements after a
# title or guide edit. A caller's subsequent alignmode edit releases this padding.
function _addon_panel_padding!(p; minimum_decoration = (0.0, 0.0, 0.0, 0.0))
    managed=p.addon_state.frame_padding
    for (axis, record) in managed
        axis.alignmode[]==record.applied && (axis.alignmode[]=record.original)
    end
    measurements=map(collect(values(p.addon_state.panel_data))) do data
        axis=data.axis
        data.panel.layout===nothing && return nothing
        panel=data.panel.layout.layoutobservables.computedbbox[]
        frame=axis.layoutobservables.computedbbox[]
        lower=frame.origin-panel.origin
        upper=panel.origin+panel.widths-frame.origin-frame.widths
        (; axis, decoration = (lower[1], upper[1], lower[2], upper[2]))
    end
    filter!(!isnothing, measurements)
    isempty(measurements) && return (0.0, 0.0, 0.0, 0.0)
    maximum_decoration=ntuple(i -> max(minimum_decoration[i], maximum(m -> m.decoration[i], measurements)), 4)
    for m in measurements
        axis=m.axis
        original=haskey(managed, axis) ? managed[axis].original : axis.alignmode[]
        # Respect constructor/native alignmode choices beyond ordinary Inside.
        axis.alignmode[] isa GridLayoutBase.Inside || continue
        native=axis.layoutobservables.protrusions[]
        values=ntuple(i -> Float32(getfield(native, i)+max(0.0, maximum_decoration[i]-m.decoration[i])), 4)
        applied=Mixed(;
            (key=>GridLayoutBase.Protrusion(value)
        for (key, value) in zip((:left, :right, :bottom, :top), values))...)
        axis.alignmode[]=applied
        managed[axis]=(; original, applied)
    end
    return maximum_decoration
end

# Measure the physical frame before native integer-pixel rounding. The allocated
# Axis rectangle can include unused space imposed by DataAspect/AxisAspect.
function _addon_frame_size(axis, allocation = Tuple(axis.layoutobservables.computedbbox[].widths))
    w, h=Float64.(allocation)
    aspect=axis.aspect[]
    ratio=aspect isa DataAspect ?
          Float64(axis.finallimits[].widths[1]/axis.finallimits[].widths[2]) :
          aspect isa Makie.AxisAspect ? Float64(aspect.aspect) : nothing
    return ratio===nothing ? (w, h) : (min(w, h*ratio), min(h, w/ratio))
end

# Native width/height refer to the Axis inner size. Outside/Mixed margins
# consume part of its assigned track; Relative sizes consume a fraction of it.
function _addon_axis_insets(axis)
    mode=axis.alignmode[]
    native=axis.layoutobservables.protrusions[]
    sides=ntuple(4) do i
        mode isa Outside && return Float64(getfield(native, i)+getfield(mode.padding, i))
        mode isa Mixed || return 0.0
        side=getfield(mode.sides, i)
        return side isa Real ? Float64(getfield(native, i)+side) : 0.0
    end
    return (sides[1]+sides[2], sides[3]+sides[4])
end

function _addon_axis_span(axis, frame)
    insets=_addon_axis_insets(axis)
    return ntuple(2) do d
        size=getproperty(axis, d==1 ? :width : :height)[]
        inner=size isa Real ? Float64(size) : frame[d]
        (inner+insets[d])/(size isa Relative ? size.x : 1.0)
    end
end

# Native autosize includes protrusions/padding without tight_bbox's one-pixel
# minimum for zero-sized reserved guide tracks. Opaque native content can leave
# a dimension undetermined; its current declared allocation is then retained.
function _addon_required_size(grid)
    automatic=grid.layoutobservables.autosize[]
    allocated=grid.layoutobservables.computedbbox[].widths
    return ntuple(d -> Float64(something(automatic[d], allocated[d])), 2)
end

function _addon_panel_size(data)
    required=_addon_required_size(data.panel.layout)
    # GridLayout's automatic track measurement reads native protrusions; its
    # placement also honors Mixed protrusion overrides. Include the latter's
    # actual extra reservation in a fitted cell, without counting either twice.
    axis=data.axis
    reported=axis.layoutobservables.reporteddimensions[].outer
    native=axis.layoutobservables.protrusions[]
    extra=(reported.left-native.left+reported.right-native.right,
        reported.bottom-native.bottom+reported.top-native.top)
    return required .+ extra
end

function _addon_frame_measurements(p)
    state=p.addon_state
    canvas=state.shell.canvas
    decorations=map(collect(values(state.panel_data))) do data
        panel=_addon_panel_size(data)
        frame=_addon_frame_size(data.axis)
        panel .- _addon_axis_span(data.axis, frame)
    end
    decoration=ntuple(i -> maximum(x -> x[i], decorations; init = 0.0), 2)
    gaps=(Float64(canvas.default_colgap.x), Float64(canvas.default_rowgap.x))
    # Both grids now report constrained content. Flexible allocation outside
    # that content is not a decoration and cannot become a sizing budget.
    outer=_addon_required_size(p.figure.layout) .- _addon_required_size(canvas)
    return (; decoration, gaps, outer)
end

function _addon_resize_preserving_views!(p, size)
    views=[axis.targetlimits[] for axis in p.axes]
    try
        resize!(p.figure, round.(Int, size)...)
        _addon_compose_guides!(p)
    finally
        for (axis, view) in zip(p.axes, views)
            axis.targetlimits[]==view || (axis.targetlimits[]=view)
        end
    end
    return p
end

# Constrain only shell-owned tracks around current per-panel frames. Native
# grids measure axis protrusions, titles, panel guides, figure guides and chrome.
# Empty internal matrix tracks inherit the common cell span; no axis is added.
function _addon_fit_frames!(p, frames; canvas = nothing)
    state=p.addon_state
    shell=state.shell
    previous=state.fitting_geometry[]
    state.fitting_geometry[]=true
    try
        if state.panel_page!==nothing
            rows, columns=state.panel_page.dimensions
            widths=zeros(columns)
            heights=zeros(rows)
            for data in values(state.panel_data)
                panel=data.panel.layout
                frame=_addon_axis_span(data.axis, frames[data.axis])
                panel.width[]=Auto()
                panel.height[]=Auto()
                panel.tellwidth[]=true
                panel.tellheight[]=true
                colsize!(panel, 2, Fixed(frame[1]))
                rowsize!(panel, 2, Fixed(frame[2]))
            end
            _addon_compose_guides!(p)
            for data in values(state.panel_data)
                panel=data.panel.layout
                gc=GridLayoutBase.gridcontent(panel)
                bounds=_addon_panel_size(data)
                panel.width[]=bounds[1]
                panel.height[]=bounds[2]
                row=first(gc.span.rows)
                column=first(gc.span.cols)
                widths[column]=max(widths[column], bounds[1])
                heights[row]=max(heights[row], bounds[2])
            end
            # Calibrated matrix frames share decorated cell spans. Physical
            # aspect galleries retain their individual row/column requirements.
            if !isempty(state.frame_padding)
                widths.=maximum(widths)
                heights.=maximum(heights)
            end
            for column in 1:columns
                colsize!(shell.canvas, column, Fixed(iszero(widths[column]) ?
                                                     maximum(widths) : widths[column]))
            end
            for row in 1:rows
                rowsize!(shell.canvas, row, Fixed(iszero(heights[row]) ? maximum(heights) :
                                                  heights[row]))
            end
            shell.canvas.width[]=Auto()
            shell.canvas.height[]=Auto()
        else
            # plotwindow's declared native canvas is opaque. Preserve its own
            # topology and allocation; only fit the surrounding managed shell.
            canvas===nothing &&
                (canvas=Tuple(shell.canvas.layoutobservables.computedbbox[].widths))
            shell.canvas.width[]=canvas[1]
            shell.canvas.height[]=canvas[2]
        end
        shell.canvas.tellwidth[]=true
        shell.canvas.tellheight[]=true
        shell.body.width[]=Auto()
        shell.body.height[]=Auto()
        shell.body.tellwidth[]=true
        shell.body.tellheight[]=true
        colsize!(shell.body, 2, Auto(true))
        rowsize!(shell.body, 2, Auto(true))
        body_row=first(GridLayoutBase.gridcontent(shell.body).span.rows)
        rowsize!(shell.root, body_row, Auto(true))
        colsize!(shell.root, 1, Auto(true))
        _addon_compose_guides!(p)
        required=_addon_required_size(shell.root)
        # Width settles measured wrapping before the native required height.
        _addon_resize_preserving_views!(p, (
            required[1], p.figure.scene.viewport[].widths[2]))
        required=_addon_required_size(shell.root)
        _addon_resize_preserving_views!(p, required)
    finally
        state.fitting_geometry[]=previous
    end
    return p
end

# A native resize is a new local allocation, not another content-fit request.
# Release the shell's fitted tracks before the existing aspect/reflow behavior.
function _addon_release_frames!(p)
    state=p.addon_state
    shell=state.shell
    if state.panel_page!==nothing
        for data in values(state.panel_data)
            panel=data.panel.layout
            panel.width[]=Relative(1)
            panel.height[]=Relative(1)
            panel.tellwidth[]=false
            panel.tellheight[]=false
            colsize!(panel, 2, Auto(false, 1))
            rowsize!(panel, 2, Auto(false, 1))
        end
        for row in 1:size(shell.canvas)[1]
            rowsize!(shell.canvas, row, Auto(false, 1))
        end
        for column in 1:size(shell.canvas)[2]
            colsize!(shell.canvas, column, Auto(false, 1))
        end
    end
    shell.canvas.width[]=Auto()
    shell.canvas.height[]=Auto()
    shell.canvas.tellwidth[]=false
    shell.canvas.tellheight[]=false
    shell.body.width[]=Auto()
    shell.body.height[]=Auto()
    shell.body.tellwidth[]=false
    shell.body.tellheight[]=false
    colsize!(shell.body, 2, Auto(false, 1))
    rowsize!(shell.body, 2, Auto(false, 1))
    rowsize!(shell.root, first(GridLayoutBase.gridcontent(shell.body).span.rows), Auto(false, 1))
    colsize!(shell.root, 1, Auto(false, 1))
    return p
end

function _addon_frame_budget(p, capacity)
    br, bc=capacity
    shell=p.addon_state.shell
    measured=_addon_frame_measurements(p)
    frames=((shell.reference_size[1]-measured.outer[1]-(bc-1)*measured.gaps[1])/bc-measured.decoration[1],
        (shell.reference_size[2]-measured.outer[2]-(br-1)*measured.gaps[2])/br-measured.decoration[2])
    all(x -> isfinite(x) && x>1, frames) || throw(ArgumentError(
        "the nominal figure size cannot fit layout=$capacity and its measured decorations; increase the figure size or reduce layout"))
    return frames
end

function _addon_calibrate_frames!(pages, capacity)
    isempty(pages) && return pages
    if all(p -> all(axis -> axis.aspect[]===nothing, p.axes), pages)
        measured=map(_addon_panel_padding!, pages)
        padding=ntuple(i -> maximum(x -> x[i], measured), 4)
        for p in pages
            _addon_panel_padding!(p; minimum_decoration = padding)
            _addon_fit_frames!(p, _addon_frame_snapshot(p).frames)
        end
    end
    candidates=map(p -> _addon_frame_budget(p, capacity), pages)
    frames=ntuple(i -> minimum(x -> x[i], candidates), 2)
    for p in pages
        physical=IdDict{Any, Tuple{Float64, Float64}}()
        for axis in p.axes
            insets=_addon_axis_insets(axis)
            inner=ntuple(2) do d
                size=getproperty(axis, d==1 ? :width : :height)[]
                size isa Real ? Float64(size) :
                frames[d]*(size isa Relative ? size.x : 1.0)-insets[d]
            end
            physical[axis]=_addon_frame_size(axis, inner)
        end
        _addon_fit_frames!(p, physical)
    end
    return pages
end

# Physical aspect is a panel requirement. Fit native panel cells into the
# available canvas using each panel's data ratio, including wide systems and
# one-row collections. No figure or window aspect ratio is imposed.
function _addon_fit_panel_aspects!(p)
    isempty(p.axes) && return p
    all(axis -> axis.aspect[] isa DataAspect, p.axes) || return p
    state=p.addon_state
    shell=state.shell
    panels=collect(values(state.panel_data))
    all(data -> data.panel.layout!==nothing, panels) || return p
    rows, columns=size(shell.canvas)
    column_ratios=zeros(columns)
    column_decoration=zeros(columns)
    row_decoration=zeros(rows)
    for data in panels
        gc=GridLayoutBase.gridcontent(data.panel.layout)
        row=first(gc.span.rows)
        column=first(gc.span.cols)
        view=data.axis.targetlimits[]
        ratio=Float64(view.widths[1]/view.widths[2])
        isfinite(ratio) && ratio>0 || continue
        column_ratios[column]=max(column_ratios[column], ratio)
        panel=data.panel.layout.layoutobservables.computedbbox[]
        frame=data.axis.layoutobservables.computedbbox[]
        decoration=panel.widths-frame.widths
        column_decoration[column]=max(column_decoration[column], decoration[1])
        row_decoration[row]=max(row_decoration[row], decoration[2])
    end
    any(iszero, column_ratios) && return p
    body=shell.body.layoutobservables
    canvas=shell.canvas.layoutobservables.computedbbox[]
    outside=max.(0.0, Float64.(body.computedbbox[].widths-canvas.widths))
    available=Float64.(body.suggestedbbox[].widths) .- outside
    gaps=(Float64(shell.canvas.default_colgap.x), Float64(shell.canvas.default_rowgap.x))
    frame_height=min(
        (available[1]-sum(column_decoration)-(columns-1)*gaps[1])/sum(column_ratios),
        (available[2]-sum(row_decoration)-(rows-1)*gaps[2])/rows)
    isfinite(frame_height) && frame_height>1 || return p # native zero-sized initialization
    column_widths=frame_height .* column_ratios .+ column_decoration
    row_heights=frame_height .+ row_decoration
    shell.canvas.width[]=sum(column_widths)+(columns-1)*gaps[1]
    shell.canvas.height[]=sum(row_heights)+(rows-1)*gaps[2]
    shell.canvas.tellwidth[]=true
    shell.canvas.tellheight[]=true
    shell.body.width[]=Auto()
    shell.body.height[]=Auto()
    shell.body.tellwidth[]=true
    shell.body.tellheight[]=true
    shell.body.halign[]=:center
    shell.body.valign[]=:center
    colsize!(shell.body, 2, Auto(true))
    rowsize!(shell.body, 2, Auto(true))
    for column in 1:columns
        colsize!(shell.canvas, column, Fixed(column_widths[column]))
    end
    for row in 1:rows
        rowsize!(shell.canvas, row, Fixed(row_heights[row]))
    end
    return p
end

function _addon_frame_snapshot(p)
    frames=IdDict(axis=>_addon_frame_size(axis) for axis in p.axes)
    canvas=Tuple(p.addon_state.shell.canvas.layoutobservables.computedbbox[].widths)
    views=[axis.targetlimits[] for axis in p.axes]
    protrusions=[axis.layoutobservables.protrusions[] for axis in p.axes]
    return (; frames, canvas, views, protrusions)
end

Base.@noinline function _addon_edit_presentation!(action, p; before = nothing)
    state=p.addon_state
    state.presentation_ready[] || return action()
    state.fitting_geometry[] && return action()
    (; frames, canvas, views, protrusions)=before===nothing ? _addon_frame_snapshot(p) :
                                           before
    state.fitting_geometry[]=true
    try
        result=action()
        _addon_compose_guides!(p)
        if state.panel_page!==nothing && all(axis -> axis.aspect[]===nothing, p.axes) &&
           protrusions!=[axis.layoutobservables.protrusions[] for axis in p.axes]
            _addon_panel_padding!(p)
        end
        _addon_fit_frames!(p, frames; canvas)
        return result
    finally
        for (axis, view) in zip(p.axes, views)
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
    on(p.figure.scene, p.figure.scene.viewport) do _
        p.addon_state.fitting_geometry[] && return nothing
        grid=state.shell.canvas
        box=state.shell.body.layoutobservables.suggestedbbox[]
        all(>(0), box.widths) || return nothing
        columns=clamp(round(Int, sqrt(length(panels)*box.widths[1]/box.widths[2])), 1, length(panels))
        columns==current_columns[] && return nothing
        views=[data.axis.targetlimits[] for data in panels]
        state.fitting_geometry[]=true
        try
            # Empty matrix-selection cells are not part of flow membership.
            if haskey(state, :page_cells)
                used=Set(data.panel.layout for data in panels)
                for cell in values(state.page_cells)
                    cell.layout in used || _addon_delete_subtree!(cell.layout)
                end
                empty!(state.page_cells)
            end
            rows=cld(length(panels), columns)
            for (index, data) in enumerate(panels)
                row, column=cld(index, columns), mod1(index, columns)
                grid[row, column]=data.panel.layout
                for key in (:xlabelvisible, :xticklabelsvisible, :xticksvisible)
                    haskey(state.shell.axis_attributes, key) ||
                        (getproperty(data.axis, key)[]=row==rows)
                end
            end
            GridLayoutBase.trim!(grid)
            for row in 1:rows
                rowsize!(grid, row, Auto(false, 1))
            end
            for column in 1:columns
                colsize!(grid, column, Auto(false, 1))
            end
            current_columns[]=columns
            p.addon_state=merge(p.addon_state,
                (panel_page = merge(p.addon_state.panel_page, (dimensions = (
                    rows, columns),)),))
            _addon_fit_panel_aspects!(p)
            _addon_compose_guides!(p)
        finally
            for (data, view) in zip(panels, views)
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
function _addon_watch_presentation!(p, block, names)
    subscriptions=Any[]
    for name in names
        name in propertynames(typeof(block)) || continue
        attribute=getproperty(block, name)
        attribute isa Observable || continue
        pending=Ref{Any}(nothing)
        push!(subscriptions,
            on(block.blockscene, attribute; priority = 1) do _
                state=p.addon_state
                pending[]=state.presentation_ready[] && !state.fitting_geometry[] &&
                          !state.composing_guides[] ? _addon_frame_snapshot(p) : nothing
                nothing
            end)
        push!(subscriptions, on(block.blockscene, attribute; priority = -100) do _
            before=pending[]
            pending[]=nothing
            before===nothing || _addon_edit_presentation!(() -> nothing, p; before)
            nothing
        end)
    end
    return subscriptions
end
