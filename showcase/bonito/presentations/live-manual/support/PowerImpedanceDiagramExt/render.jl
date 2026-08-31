const _DEFAULT_STYLE = (
    figure_size = (900, 560),
    background = :white,
    ac_bus_color = :steelblue4,
    dc_bus_color = :darkorange2,
    line_color = :gray35,
    converter_color = :mediumpurple3,
    source_color = :seagreen3,
    machine_color = :dodgerblue3,
    load_color = :darkorange2,
    shunt_color = :goldenrod2,
    generic_color = :gray55,
    inactive_color = :gray75,
    highlight_color = :orangered2,
    label_color = :gray15
)

struct NetworkDiagram{F, A, P, M, Q, S, D, H}
    figure::F
    axis::A
    plots::P
    model::M
    positions::Q
    selected::S
    details::D
    highlight::H
end

function _diagram_style(overrides::NamedTuple)
    unknown = setdiff(propertynames(overrides), propertynames(_DEFAULT_STYLE))
    isempty(unknown) || throw(
        ArgumentError("unknown diagram style field(s): $(join(string.(unknown), ", "))"),
    )
    return merge(_DEFAULT_STYLE, overrides)
end

function _external_position_key(key)
    key isa DiagramEntityKey && return key
    key isa Symbol && return ComponentKey(key)
    if key isa Tuple && length(key) == 2 && key[1] isa Integer && key[2] isa Integer
        return BusKey(Int(key[1]), Int(key[2]))
    end
    throw(
        ArgumentError(
        "explicit diagram position keys must be BusKey, ComponentKey, " *
        "(domain, bus), or an element Symbol; received $(repr(key))",
    ),
    )
end

function _point2(position, key)
    applicable(length, position) && length(position) == 2 || throw(
        ArgumentError("explicit position for $(repr(key)) must contain two coordinates"),
    )
    x, y = position[1], position[2]
    x isa Real && y isa Real || throw(
        ArgumentError("explicit position for $(repr(key)) must be numeric"),
    )
    all(isfinite, (x, y)) || throw(
        ArgumentError("explicit position for $(repr(key)) must be finite"),
    )
    return Makie.Point2f(x, y)
end

function _explicit_positions(model::DiagramModel, supplied)
    supplied isa AbstractDict || throw(
        ArgumentError("positions must be nothing or a complete explicit position mapping"),
    )
    copied = Dict{DiagramEntityKey, Makie.Point2f}()
    for (raw_key, value) in pairs(supplied)
        key = _external_position_key(raw_key)
        haskey(model.key_to_vertex, key) || throw(
            ArgumentError("explicit positions contain unknown diagram key $(repr(raw_key))"),
        )
        haskey(copied, key) && throw(
            ArgumentError("explicit positions specify $(repr(raw_key)) more than once"),
        )
        copied[key] = _point2(value, raw_key)
    end
    missing = [key for key in model.vertex_keys if !haskey(copied, key)]
    isempty(missing) || throw(
        ArgumentError("explicit positions are incomplete; missing $(join(repr.(missing), ", "))"),
    )
    return copied
end

function _layout_positions(model::DiagramModel, layout)
    isempty(model.vertex_keys) && return Dict{DiagramEntityKey, Makie.Point2f}()
    raw = if layout isa AbstractVector
        copy(layout)
    elseif applicable(layout, model.graph)
        layout(model.graph)
    else
        throw(ArgumentError("layout must accept the diagram's Graphs.jl graph"))
    end
    length(raw) == length(model.vertex_keys) || throw(
        ArgumentError(
        "layout returned $(length(raw)) positions for " *
        "$(length(model.vertex_keys)) diagram vertices",
    ),
    )
    return Dict{DiagramEntityKey, Makie.Point2f}(
        key => _point2(raw[index], key)
    for (index, key) in enumerate(model.vertex_keys)
    )
end

function _resolve_positions(model::DiagramModel, layout, positions)
    return positions === nothing ? _layout_positions(model, layout) :
           _explicit_positions(model, positions)
end

function _component_lookup(model)
    return Dict(component.key => component for component in model.components)
end

function _edge_style(model, style)
    components = _component_lookup(model)
    colors = Any[]
    widths = Float32[]
    for edge in Graphs.edges(model.graph)
        left = model.vertex_keys[Graphs.src(edge)]
        right = model.vertex_keys[Graphs.dst(edge)]
        component_key = left isa ComponentKey ? left : right
        component = components[component_key]
        color = if component.active === false
            style.inactive_color
        elseif component.role == :dc_branch
            style.dc_bus_color
        elseif component.role == :converter
            style.converter_color
        else
            style.line_color
        end
        push!(colors, color)
        push!(widths, component.role in (:ac_branch, :dc_branch) ? 2.4f0 : 1.8f0)
    end
    return colors, widths
end

_bus_label(key::BusKey) = "$(key.domain == 1 ? "AC" : "DC") $(key.bus)"

function _bus_plot!(axis, descriptor, position, style)
    points = _bus_points(descriptor, position)
    if descriptor.key.domain == 1
        return Makie.linesegments!(axis, points; color = style.ac_bus_color, linewidth = 8)
    end
    return Makie.linesegments!(axis, points; color = style.dc_bus_color, linewidth = 4)
end

function _bus_points(descriptor, position)
    halfwidth = 0.13f0
    if descriptor.key.domain == 1
        return [
            position - Makie.Vec2f(halfwidth, 0),
            position + Makie.Vec2f(halfwidth, 0)
        ]
    end
    offset = 0.025f0
    return [
        position - Makie.Vec2f(halfwidth, offset),
        position + Makie.Vec2f(halfwidth, -offset),
        position - Makie.Vec2f(halfwidth, -offset),
        position + Makie.Vec2f(halfwidth, offset)
    ]
end

function _component_appearance(component, style)
    component.active === false &&
        return (:circle, 9, style.inactive_color, style.inactive_color, 0)
    component.role == :converter &&
        return (:diamond, 24, :white, style.converter_color, 3)
    component.role == :transformer &&
        return (:circle, 22, :white, style.ac_bus_color, 3)
    component.role == :source &&
        return (:circle, 20, style.source_color, style.line_color, 1.5)
    component.role == :machine &&
        return (:circle, 20, style.machine_color, style.line_color, 1.5)
    component.role == :load &&
        return (:dtriangle, 20, style.load_color, style.line_color, 1)
    component.role == :shunt &&
        return (:utriangle, 18, style.shunt_color, style.line_color, 1)
    component.role == :generic &&
        return (:rect, 16, style.generic_color, style.line_color, 1)
    line_color = component.role == :dc_branch ? style.dc_bus_color : style.line_color
    return (:circle, 3, line_color, line_color, 0)
end

function _component_plot!(axis, component, position, style)
    marker, markersize, color, strokecolor, strokewidth = _component_appearance(component, style)
    return Makie.scatter!(
        axis,
        [position];
        marker,
        markersize,
        color,
        strokecolor,
        strokewidth
    )
end

const _ROLE_LABELS = Dict(
    :ac_branch => "AC line / cable",
    :converter => "Converter",
    :dc_branch => "DC branch",
    :generic => "Component",
    :load => "Load",
    :machine => "Machine",
    :shunt => "Shunt",
    :source => "Source",
    :transformer => "Transformer"
)

function _legend_element(role::Symbol, style)
    if role in (:ac_branch, :dc_branch)
        color = role === :dc_branch ? style.dc_bus_color : style.line_color
        return Makie.LineElement(color = color, linewidth = 3)
    end
    descriptor = (; role, active = true)
    marker, markersize, color, strokecolor, strokewidth =
        _component_appearance(descriptor, style)
    return Makie.MarkerElement(
        marker = marker,
        markersize = max(markersize, 12),
        color = color,
        strokecolor = strokecolor,
        strokewidth = strokewidth
    )
end

function _component_legend!(figure, model, style)
    roles = unique(
        component.role for component in model.components
        if component.active !== false
    )
    isempty(roles) && return nothing
    elements = [_legend_element(role, style) for role in roles]
    labels = [get(_ROLE_LABELS, role, string(role)) for role in roles]
    return Makie.Legend(
        figure[2, 1],
        elements,
        labels;
        orientation = :horizontal,
        nbanks = 1,
        tellwidth = false,
        framevisible = false,
        labelcolor = style.label_color,
        labelsize = 11,
        padding = (0, 0, 0, 0),
        patchsize = (22, 12),
        colgap = 12
    )
end

function _entity_details(model::DiagramModel, key::BusKey)
    descriptor = only(bus for bus in model.buses if bus.key == key)
    return (
        kind = :bus,
        key = descriptor.key,
        domain = descriptor.key.domain,
        bus = descriptor.key.bus,
        nodes = descriptor.nodes,
        input = descriptor.input,
        solution = descriptor.solution
    )
end

function _entity_details(model::DiagramModel, key::ComponentKey)
    descriptor = only(component for component in model.components if component.key == key)
    setpoint = model.operating_point === nothing ? nothing :
               get(model.operating_point, key.element, nothing)
    return (
        kind = :component,
        key = descriptor.key,
        element = descriptor.key.element,
        role = descriptor.role,
        connected_buses = descriptor.connected_buses,
        grounded_sides = descriptor.grounded_sides,
        terminal_rows = descriptor.terminal_rows,
        model_type = descriptor.model_type,
        pmtype = descriptor.pmtype,
        pmindex = descriptor.pmindex,
        active = descriptor.active,
        input = descriptor.input,
        solution = descriptor.solution,
        operating_point = setpoint
    )
end

function select!(view::NetworkDiagram, key::Union{Nothing, DiagramEntityKey})
    key === nothing || haskey(view.model.key_to_vertex, key) ||
        throw(
            ArgumentError("cannot select unknown diagram entity $(repr(key))"),
        )
    view.selected[] = key
    view.details[] = key === nothing ? (; kind = :none) : _entity_details(view.model, key)
    view.highlight[] = key === nothing ? Makie.Point2f[] : [view.positions[key]]
    return view
end

function _render_diagram_on_axis!(
        figure,
        axis,
        model::DiagramModel,
        resolved_positions;
        interactive::Bool,
        style::NamedTuple
)
    position_vector = [resolved_positions[key] for key in model.vertex_keys]
    Makie.hidedecorations!(axis)
    Makie.hidespines!(axis)

    edge_colors, edge_widths = _edge_style(model, style)
    graph_plot = GraphMakie.graphplot!(
        axis,
        model.graph;
        layout = position_vector,
        node_color = Makie.RGBAf(1, 1, 1, 0.001),
        node_size = fill(30, length(model.vertex_keys)),
        node_strokewidth = 0,
        edge_color = edge_colors,
        edge_width = edge_widths,
        arrow_show = false
    )

    bus_plots = Any[]
    bus_labels = Any[]
    for bus in model.buses
        position = resolved_positions[bus.key]
        push!(bus_plots, _bus_plot!(axis, bus, position, style))
        push!(
            bus_labels,
            Makie.text!(
                axis,
                [position + Makie.Vec2f(0, 0.09)];
                text = [_bus_label(bus.key)],
                align = (:center, :bottom),
                fontsize = 12,
                color = style.label_color
            )
        )
    end

    component_plots = Any[]
    component_labels = Any[]
    for component in model.components
        position = resolved_positions[component.key]
        push!(
            component_plots,
            _component_plot!(axis, component, position, style)
        )
        push!(
            component_labels,
            Makie.text!(
                axis,
                [position + Makie.Vec2f(0, 0.055)];
                text = [String(component.key.element)],
                align = (:center, :bottom),
                fontsize = 11,
                color = component.active === false ?
                        style.inactive_color : style.label_color
            )
        )
    end

    highlight = Makie.Observable(Makie.Point2f[])
    highlight_plot = Makie.scatter!(
        axis,
        highlight;
        marker = :circle,
        markersize = 38,
        color = :transparent,
        strokecolor = style.highlight_color,
        strokewidth = 2.5
    )
    legend = _component_legend!(figure, model, style)
    Makie.translate!(GraphMakie.get_node_plot(graph_plot), 0, 0, 100)

    selected = Makie.Observable{Union{Nothing, BusKey, ComponentKey}}(nothing)
    details = Makie.Observable{Any}((; kind = :none))
    plots = (
        graph = graph_plot,
        buses = bus_plots,
        components = component_plots,
        bus_labels = bus_labels,
        component_labels = component_labels,
        highlight = highlight_plot,
        legend
    )
    view = NetworkDiagram(
        figure,
        axis,
        plots,
        model,
        resolved_positions,
        selected,
        details,
        highlight
    )

    if interactive
        action = (index, _...) -> begin
            select!(view, model.vertex_keys[index])
            return true
        end
        Makie.register_interaction!(
            axis,
            :powerimpedance_diagram_select,
            GraphMakie.NodeClickHandler(action)
        )
    end
    Makie.autolimits!(axis)
    return view
end

"""Update selected diagram entity positions without rebuilding its Makie scene."""
function update_positions!(view::NetworkDiagram, updates; autolimits::Bool = false)
    for (raw_key, raw_position) in pairs(updates)
        key = _external_position_key(raw_key)
        haskey(view.positions, key) || throw(KeyError(key))
        view.positions[key] = _point2(raw_position, raw_key)
    end

    position_vector = [view.positions[key] for key in view.model.vertex_keys]
    view.plots.graph.layout[] = position_vector

    for (descriptor, plot, label) in zip(
            view.model.buses,
            view.plots.buses,
            view.plots.bus_labels
    )
        position = view.positions[descriptor.key]
        Makie.update!(plot; arg1 = _bus_points(descriptor, position))
        Makie.update!(label; arg1 = [position + Makie.Vec2f(0, 0.09)])
    end
    for (descriptor, plot, label) in zip(
            view.model.components,
            view.plots.components,
            view.plots.component_labels
    )
        position = view.positions[descriptor.key]
        Makie.update!(plot; arg1 = [position])
        Makie.update!(label; arg1 = [position + Makie.Vec2f(0, 0.055)])
    end
    selected = view.selected[]
    selected === nothing || (view.highlight[] = [view.positions[selected]])
    autolimits && Makie.autolimits!(view.axis)
    return view
end

function _diagram_plots(view::NetworkDiagram)
    return Any[
        view.plots.graph,
        view.plots.buses...,
        view.plots.components...,
        view.plots.bus_labels...,
        view.plots.component_labels...,
        view.plots.highlight
    ]
end

function _render_diagram(
        model::DiagramModel;
        layout,
        positions,
        interactive::Bool,
        style::NamedTuple
)
    resolved_style = _diagram_style(style)
    resolved_positions = _resolve_positions(model, layout, positions)
    figure = Makie.Figure(
        size = resolved_style.figure_size,
        backgroundcolor = resolved_style.background
    )
    axis = Makie.Axis(
        figure[1, 1];
        aspect = Makie.DataAspect(),
        xautolimitmargin = (0.12, 0.12),
        yautolimitmargin = (0.12, 0.12)
    )
    return _render_diagram_on_axis!(
        figure,
        axis,
        model,
        resolved_positions;
        interactive,
        style = resolved_style
    )
end
