# Mesh rendering and inspection own only detached mesh data and native Makie
# resources. Display indices never replace the original node and element tags.
function _fem_mesh_drawing(mesh)
    _fem_planar(mesh.coordinates)
    elements = Tuple{Int, Int}[]
    edges = Tuple{Int, Int}[]
    edge_elements = Vector{Int}[]
    lookup = Dict{Tuple{Int, Int}, Int}()
    for (b, block) in enumerate(mesh.blocks)
        block.dimension <= 2 || throw(ArgumentError("volume mesh plotting is not supported"))
        block.dimension == 0 || block.order == 1 || throw(ArgumentError(
            "mesh plotting requires first-order elements; higher-order geometry needs tessellation"))
        block.primary_nodes in (1, 2, 3, 4) || throw(ArgumentError("unsupported planar mesh element"))
        for (column, nodes) in enumerate(eachcol(block.connectivity))
            push!(elements, (b, column))
            block.dimension == 0 && continue
            for i in 1:(length(nodes) == 2 ? 1 : length(nodes))
                edge = minmax(nodes[i], nodes[mod1(i+1, length(nodes))])
                index = get!(lookup, edge) do
                    push!(edges, edge)
                    push!(edge_elements, Int[])
                    length(edges)
                end
                push!(edge_elements[index], length(elements))
            end
        end
    end
    # An actual line element takes precedence over adjacent display edges.
    for owners in edge_elements
        sort!(owners; by = i -> begin
            b, c = elements[i]
            (mesh.blocks[b].dimension, mesh.blocks[b].element_tags[c])
        end)
    end
    return (; elements, edges, edge_elements)
end

_fem_mesh_point(mesh, i) = GeometryBasics.Point2d(mesh.coordinates[1, i], mesh.coordinates[2, i])

function _fem_mesh_memberships(mesh, dimension, tags)
    isempty(tags) && return "Ungrouped"
    return join([begin
        name = get(mesh.physical_names, (dimension, tag), "")
        "$dimension-D #$tag: " * (isempty(name) ? "Unnamed" : name)
    end for tag in sort(tags)], "; ")
end

function _fem_mesh_overview(mesh)
    counts = Dict{Tuple{Int, Int}, Int}()
    for block in mesh.blocks
        key = (block.dimension, block.element_type)
        counts[key] = get(counts, key, 0) + length(block.element_tags)
    end
    groups = Set((b.dimension, tag) for b in mesh.blocks for tag in b.physical_tags)
    rows = ["$d-D, type $t: $(counts[(d,t)])" for (d,t) in sort!(collect(keys(counts)))]
    return join([basename(mesh.source), "Nodes: $(length(mesh.node_tags))", rows...,
        "Physical groups: $(length(groups))", "Click an element or node to inspect it.",
        "Mesh x/y = displayed y/z [m]."], '\n')
end

function _fem_mesh_incidence!(view)
    view.incidence[] === nothing || return view.incidence[]
    incident = [Int[] for _ in view.mesh.node_tags]
    for (i, (b, c)) in enumerate(view.drawing.elements)
        for node in view.mesh.blocks[b].connectivity[:, c]
            push!(incident[node], i)
        end
    end
    view.incidence[] = incident
    return incident
end

function _fem_mesh_details(view, selection)
    mesh = view.mesh
    selection === nothing && return _fem_mesh_overview(mesh)
    kind, index = selection
    if kind === :node
        incident = _fem_mesh_incidence!(view)[index]
        refs = view.drawing.elements[incident]
        tags = sort!([mesh.blocks[b].element_tags[c] for (b,c) in refs])
        entities = sort!(unique((mesh.blocks[b].dimension, mesh.blocks[b].entity) for (b,_) in refs))
        groups = sort!(unique((mesh.blocks[b].dimension, tag) for (b,_) in refs
            for tag in mesh.blocks[b].physical_tags))
        return join(["Node #$(mesh.node_tags[index])",
            "Native (x, y, z) [m]: " * join((@sprintf("%.12g", x) for x in mesh.coordinates[:, index]), ", "),
            "Incident elements ($(length(tags))): " * join(tags, ", "),
            "Incident entities (dimension, tag): " * join(entities, ", "),
            "Physical groups:", isempty(groups) ? "Ungrouped" : join(
                [_fem_mesh_memberships(mesh, d, [tag]) for (d,tag) in groups], '\n'),
            "Mesh x/y = displayed y/z [m]."], '\n')
    end
    b, c = view.drawing.elements[index]
    block = mesh.blocks[b]
    nodes = block.connectivity[:, c]
    points = mesh.coordinates[:, nodes]
    name = block.dimension == 2 ? (length(nodes) == 3 ? "Triangle" : "Quadrangle") :
           block.dimension == 1 ? "Line" : "Point"
    rows = ["$name #$(block.element_tags[c])",
        "Type $(block.element_type); order $(block.order); $(block.dimension)-D entity #$(block.entity)",
        "Nodes: " * join(mesh.node_tags[nodes], ", "),
        "Centroid (x, y, z) [m]: " * join((@sprintf("%.12g", x) for x in vec(mean(points; dims=2))), ", "),
        "Physical groups:", isempty(block.physical_tags) ? "Ungrouped" : join(
            [_fem_mesh_memberships(mesh, block.dimension, [tag]) for tag in sort(block.physical_tags)], '\n')]
    if block.dimension > 0
        count = block.dimension == 1 ? 1 : length(nodes)
        lengths = [sqrt(sum(abs2, points[:, i]-points[:, mod1(i+1, length(nodes))])) for i in 1:count]
        push!(rows, @sprintf("%s [m]: %.12g", block.dimension == 1 ? "Length" : "Perimeter", sum(lengths)))
        if block.dimension == 2
            # Translate before the shoelace sum to avoid cancellation far from zero.
            xy = points[1:2, :] .- points[1:2, 1]
            area = abs(sum(xy[1,i]*xy[2,mod1(i+1,count)] - xy[2,i]*xy[1,mod1(i+1,count)] for i in 1:count))/2
            push!(rows, @sprintf("Area [m²]: %.12g", area))
            push!(rows, @sprintf("Min / max edge [m]: %.12g / %.12g", extrema(lengths)...))
        end
    end
    return join(rows, '\n')
end

function _fem_mesh_categories(mesh, mode)
    mode in (:uniform, :physical, :entity) || throw(ArgumentError(
        "color_by must be :uniform, :physical, or :entity"))
    dimension = maximum(b.dimension for b in mesh.blocks)
    keys = [(b.dimension, mode === :entity ? (b.entity,) : Tuple(sort(b.physical_tags))) for b in mesh.blocks]
    categories = sort!(unique(keys[i] for i in eachindex(keys) if mesh.blocks[i].dimension == dimension))
    indices = [something(findfirst(==(key), categories), 0) for key in keys]
    colors = Makie.to_colormap(:tab20)
    palette = [i <= length(colors) ? colors[i] : RGBA(RGB(HSV(mod(137.508i, 360),
        .55+.15isodd(i), .75+.1mod(i,3))), 1) for i in eachindex(categories)]
    labels = [mode === :entity ? "$d-D entity #$(only(tags))" : _fem_mesh_memberships(mesh, d, tags)
              for (d,tags) in categories]
    return (; indices, palette, labels)
end

function _fem_mesh_faces!(view)
    view.faces[] === nothing || return view.faces[]
    points = GeometryBasics.Point2d[]
    faces = GeometryBasics.TriangleFace{Int}[]
    for (i, (b,c)) in enumerate(view.drawing.elements)
        block = view.mesh.blocks[b]
        block.dimension == 2 || continue
        nodes = block.connectivity[:, c]
        offset = length(points)
        append!(points, [_fem_mesh_point(view.mesh, n) for n in nodes])
        append!(view.face_owners, fill(i, length(nodes)))
        push!(faces, GeometryBasics.TriangleFace{Int}(offset+1, offset+2, offset+3))
        length(nodes) == 4 && push!(faces, GeometryBasics.TriangleFace{Int}(offset+1, offset+3, offset+4))
    end
    isempty(faces) && return nothing
    view.faces[] = Makie.mesh!(view.axis, points, faces; color=fill(RGBA(.6,.6,.6,.12), length(points)),
        shading=Makie.NoShading, inspectable=false)
    translate!(view.faces[], 0, 0, view.depth-1)
    return view.faces[]
end

function _fem_mesh_nodes!(view)
    view.nodes[] === nothing || return view.nodes[]
    # Draw the lowest native tag last at coincident points, preserving identities.
    append!(view.node_owners, sortperm(view.mesh.node_tags; rev=true))
    points = [_fem_mesh_point(view.mesh, i) for i in view.node_owners]
    view.nodes[] = scatter!(view.axis, points; color=:black, markersize=5, inspectable=false)
    translate!(view.nodes[], 0, 0, view.depth+1)
    return view.nodes[]
end

function _fem_mesh_update!(view)
    colored = view.color_by[] !== :uniform
    active = view.mode[] !== :none
    categories = _fem_mesh_categories(view.mesh, view.color_by[])
    if colored || view.mode[] === :element
        _fem_mesh_faces!(view)
    end
    if view.faces[] !== nothing
        view.faces[].color = colored ? [categories.palette[categories.indices[view.drawing.elements[i][1]]]
            for i in view.face_owners] : fill(RGBA(.6,.6,.6,.12), length(view.face_owners))
        view.faces[].visible = view.wireframe.visible[] && (colored || view.mode[] === :element)
    end
    if view.mode[] === :node
        _fem_mesh_nodes!(view)
    end
    view.nodes[] === nothing || (view.nodes[].visible = view.wireframe.visible[] && view.mode[] === :node)
    dimension = maximum(b.dimension for b in view.mesh.blocks)
    if dimension == 1
        view.wireframe.color = colored ? [categories.palette[categories.indices[
            view.drawing.elements[first(owners)][1]]] for owners in view.drawing.edge_elements for _ in 1:2] :
            fill(Makie.to_color(view.color), 2length(view.drawing.edges))
    elseif dimension == 0
        colors = fill(Makie.to_color(view.color), length(view.mesh.node_tags))
        if colored
            for (b, block) in enumerate(view.mesh.blocks), node in block.connectivity
                colors[node] = categories.palette[categories.indices[b]]
            end
        end
        view.wireframe.color = colors
    end
    view.categories[] = categories
    view.outline.visible = view.wireframe.visible[] && active && view.selection[] !== nothing && first(view.selection[]) === :element
    point_selected = view.selection[] !== nothing && (first(view.selection[]) === :node ||
        view.mesh.blocks[view.drawing.elements[last(view.selection[])][1]].dimension == 0)
    view.marker.visible = view.wireframe.visible[] && active && point_selected
    return nothing
end

function _fem_mesh_select!(view, selection)
    view.selection[] = selection
    points = GeometryBasics.Point2d[]
    marker = GeometryBasics.Point2d[]
    if selection !== nothing
        kind, i = selection
        if kind === :node
            push!(marker, _fem_mesh_point(view.mesh, i))
        else
            b, c = view.drawing.elements[i]
            nodes = view.mesh.blocks[b].connectivity[:, c]
            append!(length(nodes)==1 ? marker : points, [_fem_mesh_point(view.mesh, n) for n in nodes])
            length(nodes) > 2 && push!(points, first(points))
        end
    end
    view.outline[1] = points
    view.marker[1] = marker
    view.outline.visible = !isempty(points) && view.wireframe.visible[]
    view.marker.visible = !isempty(marker) && view.wireframe.visible[]
    view.details[] = _fem_mesh_details(view, selection)
    return selection
end

function _fem_mesh_pick(view, pixel)
    view.wireframe.visible[] || return nothing
    for (plot, index) in Makie.pick_sorted(Makie.root(view.axis.scene), pixel, 6)
        # A highlight can cover its own primitive in the native picking buffer.
        if (view.mode[] !== :none && plot === view.marker) ||
                (view.mode[] === :element && plot === view.outline)
            return view.selection[]
        end
        if view.mode[] === :node && plot === view.nodes[]
            return (:node, view.node_owners[index])
        elseif view.mode[] === :element
            plot === view.faces[] && return (:element, view.face_owners[index])
            if plot === view.wireframe && !isempty(view.drawing.edges)
                return (:element, first(view.drawing.edge_elements[cld(index, 2)]))
            elseif plot === view.wireframe
                owners = _fem_mesh_incidence!(view)[index]
                isempty(owners) || return (:element, first(owners))
            end
        end
    end
    return nothing
end

function _spatial_layer!(axis, mesh::FEMMesh;
        color=(:black, .35), linewidth=.5, visible=true,
        color_by=:uniform, inspect=:none, depth=0)
    inspect in (:none, :element, :node) || throw(ArgumentError("inspect must be :none, :element, or :node"))
    _fem_mesh_categories(mesh, color_by)
    drawing = _fem_mesh_drawing(mesh)
    segments = [_fem_mesh_point(mesh, n) for edge in drawing.edges for n in edge]
    wire_color = maximum(b.dimension for b in mesh.blocks) < 2 ?
        fill(Makie.to_color(color), isempty(segments) ? length(mesh.node_tags) : length(segments)) : color
    wireframe = isempty(segments) ? scatter!(axis, [_fem_mesh_point(mesh, i) for i in eachindex(mesh.node_tags)];
        color=wire_color, visible, markersize=3, inspectable=false) :
        Makie.linesegments!(axis, segments; color=wire_color, linewidth, visible, inspectable=false)
    translate!(wireframe, 0, 0, depth)
    outline = lines!(axis, GeometryBasics.Point2d[]; color=:magenta, linewidth=3,
        inspectable=false, xautolimits=false, yautolimits=false, visible=false)
    marker = scatter!(axis, GeometryBasics.Point2d[]; color=:magenta, markersize=12,
        inspectable=false, xautolimits=false, yautolimits=false, visible=false)
    translate!(outline, 0, 0, depth+2)
    translate!(marker, 0, 0, depth+2)
    view = (; axis, mesh, drawing, wireframe, color, depth,
        faces=Ref{Any}(nothing), face_owners=Int[], nodes=Ref{Any}(nothing), node_owners=Int[],
        color_by=Observable{Symbol}(color_by), mode=Observable{Symbol}(inspect),
        selection=Observable{Union{Nothing, Tuple{Symbol, Int}}}(nothing),
        incidence=Ref{Any}(nothing), categories=Observable{Any}(nothing),
        details=Observable(_fem_mesh_overview(mesh)), outline, marker,
        interaction=gensym(:fem_mesh_inspection))
    on(axis.scene, view.mode) do mode
        _fem_mesh_select!(view, nothing)
        _fem_mesh_update!(view)
    end
    on(axis.scene, view.color_by) do mode
        _fem_mesh_update!(view)
    end
    on(axis.scene, wireframe.visible) do visible
        visible || _fem_mesh_select!(view, nothing)
        _fem_mesh_update!(view)
    end
    Makie.register_interaction!(axis, view.interaction) do event, _
        if event isa Makie.MouseEvent && event.type === Makie.MouseEventTypes.leftclick && view.mode[] !== :none
            # Axis mouse events use local pixels. picking addresses the root screen.
            _fem_mesh_select!(view, _fem_mesh_pick(view, event.px + axis.scene.viewport[].origin))
            return Makie.Consume(true)
        end
        Makie.Consume(false)
    end
    _fem_mesh_update!(view)
    return (; layers=Any[wireframe], scales=(), excluded=0, colorscale=identity,
        group=:mesh, label="Mesh", mesh_view=view)
end

# The sidebar is part of this spatial recipe, not a second shell or a global
# selection registry. Its fixed width bounds long physical names and node lists.
function _mesh_sidebar!(shell, panel, view; controls)
    _addon_detach!(view.axis)
    layout = GridLayout(panel.content; tellwidth=false, tellheight=false)
    layout[1,1] = view.axis
    aside = GridLayout(layout[1,2]; tellwidth=false, tellheight=false, valign=:top)
    colgap!(layout, 12)
    colsize!(layout, 1, Auto(false, 1))
    page = Observable(1)
    legend = Ref{Any}(nothing)
    pager = GridLayout(aside[2,1]; tellwidth=false, height=26)
    previous = Button(pager[1,1]; label="‹", width=28)
    page_label = Label(pager[1,2], ""; fontsize=10, tellwidth=false)
    following = Button(pager[1,3]; label="›", width=28)
    readout = Label(aside[3,1], ""; fontsize=11, halign=:left, valign=:top,
        justification=:left, tellwidth=false)
    GridLayout(aside[4,1]; tellheight=false)
    rowsize!(aside, 3, Auto(true))
    rowsize!(aside, 4, Auto(false, 1))
    # Full text remains available through Copy details. Summaries never pretend
    # an abbreviated connectivity or membership list is complete.
    function summarize(text)
        lines = split(text, '\n')
        shortened = length(lines) > 18 || any(length(line)>52 for line in lines)
        result = join([length(line)>52 ? first(line,49)*"…" : line for line in first(lines,min(18,length(lines)))], '\n')
        return shortened ? result*"\n… Copy details for the complete record." : result
    end
    on(view.axis.scene, view.details; update=true) do text
        readout.text = summarize(text)
    end
    function refresh()
        colored = view.color_by[] !== :uniform
        active = colored || view.mode[] !== :none
        colsize!(layout, 2, Fixed(active ? 310 : 0))
        colgap!(layout, active ? 12 : 0)
        readout.blockscene.visible[] = active
        legend[] === nothing || (delete!(legend[]); legend[] = nothing)
        categories = view.categories[]
        count = colored ? length(categories.labels) : 0
        pages = max(1, cld(count, 8))
        selected_page = clamp(page[], 1, pages)
        if count > 0
            indices = ((selected_page-1)*8+1):min(selected_page*8, count)
            labels = [length(categories.labels[i])>43 ? first(categories.labels[i],40)*"…" : categories.labels[i] for i in indices]
            entries = [Makie.PolyElement(color=categories.palette[i]) for i in indices]
            title = "Mesh " * (view.color_by[] === :physical ? "physical groups" : "entities")
            pages > 1 && (title *= " ($(first(indices))–$(last(indices)) of $count)")
            legend[] = Legend(aside[1,1], entries, labels, title;
                labelsize=10, titlesize=12, tellwidth=false, tellheight=true,
                framevisible=false, halign=:left,
                patchsize=(14,12), padding=(0,0,2,2))
        end
        rowsize!(aside, 1, count > 0 ? Auto(true) : Fixed(0))
        rowsize!(aside, 2, controls && pages > 1 ? Fixed(26) : Fixed(0))
        for control in (previous, following, page_label)
            control.blockscene.visible[] = controls && active && pages > 1
        end
        page_label.text = "$selected_page / $pages"
        return nothing
    end
    on(view.axis.scene, view.categories) do _
        refresh()
    end
    on(view.axis.scene, page) do _
        refresh()
    end
    on(view.axis.scene, previous.clicks) do _
        page[] = max(1, page[]-1)
    end
    on(view.axis.scene, following.clicks) do _
        page[] = min(max(1,cld(length(view.categories[].labels),8)),page[]+1)
    end
    refresh()
    return (; layout, aside, legend, readout, page, previous, following)
end

function _mesh_controls!(p, view, sidebar)
    # This recipe contains an axis and a sidebar in one native canvas. Use the
    # shell's existing custom-canvas path so export preserves both allocations
    # instead of constraining the entire cell to the axis frame alone.
    _addon_release_frames!(p)
    p.addon_state = merge(p.addon_state, (; mesh_inspection=(; view, sidebar), panel_page=nothing))
    # The spatial shell fits the data aspect before these controls are attached.
    # Restore the requested canvas so the sidebar does not consume that fitted axis.
    (view.color_by[] !== :uniform || view.mode[] !== :none) &&
        resize!(p.figure, p.addon_state.shell.reference_size...)
    p.addon_state.controls_enabled || return p
    _addon_widget!((_,slot) -> Makie.Menu(slot; width=175, options=[
        ("Color: Uniform",:uniform), ("Color: Physical groups",:physical), ("Color: Entities",:entity)],
        default=findfirst(==(view.color_by[]), (:uniform,:physical,:entity))), p, :mesh_color_by;
        event=menu -> menu.selection, callback=(_,value) -> (view.color_by[]=value))
    _addon_widget!((_,slot) -> Makie.Menu(slot; width=145, options=[
        ("Inspect: Off",:none), ("Inspect: Elements",:element), ("Inspect: Nodes",:node)],
        default=findfirst(==(view.mode[]), (:none,:element,:node))), p, :mesh_inspect;
        event=menu -> menu.selection, callback=(_,value) -> (view.mode[]=value))
    _addon_widget!((_,slot) -> Button(slot; label="Copy details"), p, :mesh_copy;
        event=button -> button.clicks, callback=(_,_) -> InteractiveUtils.clipboard(view.details[]),
        success="Mesh details copied")
    return p
end
