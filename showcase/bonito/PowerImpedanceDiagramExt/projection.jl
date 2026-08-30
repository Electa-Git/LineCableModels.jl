const DiagramTerminal = @NamedTuple{
    node::Symbol,
    bus::Int,
    element::Symbol,
    side::Int,
    terminal::Int,
    domain::Int
}

"Identify one physical AC or DC bus without allowing domain collisions."
struct BusKey
    domain::Int
    bus::Int

    function BusKey(domain::Int, bus::Int)
        domain in (1, 2) || throw(ArgumentError("diagram bus domain must be 1 or 2"))
        bus > 0 || throw(ArgumentError("diagram buses must have positive bus numbers"))
        return new(domain, bus)
    end
end

"Identify one named PowerImpedance component vertex."
struct ComponentKey
    element::Symbol
end

const DiagramEntityKey = Union{BusKey, ComponentKey}
const PowerModelsPayload = Union{Nothing, AbstractDict}

struct BusDescriptor
    key::BusKey
    nodes::Tuple{Vararg{Symbol}}
    input::PowerModelsPayload
    solution::PowerModelsPayload
end

struct SideBus
    side::Int
    bus::BusKey
end

struct ComponentDescriptor
    key::ComponentKey
    role::Symbol
    connected_buses::Tuple{Vararg{SideBus}}
    grounded_sides::Tuple{Vararg{Int}}
    terminal_rows::Tuple{Vararg{DiagramTerminal}}
    model_type::DataType
    pmtype::Union{Nothing, String}
    pmindex::Union{Nothing, Int}
    input::PowerModelsPayload
    solution::PowerModelsPayload
    active::Union{Nothing, Bool}
end

struct IncidenceDescriptor
    component::ComponentKey
    side::Int
    bus::BusKey
end

struct DiagramDiagnostic
    kind::Symbol
    entity::Union{Nothing, DiagramEntityKey}
    message::String
end

struct DiagramModel{O, D}
    buses::Vector{BusDescriptor}
    components::Vector{ComponentDescriptor}
    incidences::Vector{IncidenceDescriptor}
    graph::Graphs.SimpleGraph{Int}
    vertex_keys::Vector{DiagramEntityKey}
    key_to_vertex::Dict{DiagramEntityKey, Int}
    diagnostics::Vector{DiagramDiagnostic}
    operating_point::O
    powerflow_diagnostics::D
end

# PowerModels/PowerModelsACDC adapter. All schema-specific dictionary access is
# intentionally confined to this section.

function _record_get(record::AbstractDict, field::String, default = nothing)
    haskey(record, field) && return record[field]
    symbol = Symbol(field)
    haskey(record, symbol) && return record[symbol]
    return default
end

_record_get(::Nothing, ::String, default = nothing) = default

function _pm_table(container, table::String)
    container isa AbstractDict || return nothing
    return _record_get(container, table, nothing)
end

function _pm_indexed_record(container, table::String, index::Int)
    records = _pm_table(container, table)
    records isa AbstractDict || return nothing
    string_index = string(index)
    haskey(records, string_index) && return records[string_index]
    haskey(records, index) && return records[index]
    return nothing
end

function _pm_solution_root(result)
    result isa AbstractDict || return nothing
    return _record_get(result, "solution", nothing)
end

function _pm_mapping(mapping, element::Symbol)
    haskey(mapping, element) || return nothing
    entry = mapping[element]
    pmtype = if hasproperty(entry, :pmtype)
        getproperty(entry, :pmtype)
    elseif entry isa AbstractDict
        _record_get(entry, "pmtype", nothing)
    end
    pmindex = if hasproperty(entry, :compkey)
        getproperty(entry, :compkey)
    elseif hasproperty(entry, :pmindex)
        getproperty(entry, :pmindex)
    elseif entry isa AbstractDict
        compkey = _record_get(entry, "compkey", nothing)
        compkey === nothing ? _record_get(entry, "pmindex", nothing) : compkey
    end
    pmtype isa Union{String, Symbol} || throw(
        ArgumentError("power-flow mapping for :$element has no valid pmtype"),
    )
    pmindex isa Integer || throw(
        ArgumentError("power-flow mapping for :$element has no valid component index"),
    )
    return (pmtype = String(pmtype), pmindex = Int(pmindex))
end

_pm_bus_table(key::BusKey) = key.domain == 1 ? "bus" : "busdc"

function _pm_bus_record(container, key::BusKey)
    return _pm_indexed_record(container, _pm_bus_table(key), key.bus)
end

const _PM_ENDPOINT_FIELDS = Dict(
    "branch" => ((1, "f_bus"), (1, "t_bus")),
    "branchdc" => ((2, "fbusdc"), (2, "tbusdc")),
    "convdc" => ((1, "busac_i"), (2, "busdc_i")),
    "gen" => ((1, "gen_bus"),),
    "gendc" => ((2, "gen_bus"),),
    "load" => ((1, "load_bus"),),
    "shunt" => ((1, "shunt_bus"),),
    "im" => ((1, "im_bus"),)
)

const _PM_STATUS_FIELDS = Dict(
    "branch" => "br_status",
    "branchdc" => "status",
    "convdc" => "status",
    "gen" => "gen_status",
    "gendc" => "gen_status",
    "load" => "status",
    "shunt" => "status",
    "im" => "status"
)

function _pm_component_endpoints(pmtype::String, record::AbstractDict)
    fields = get(_PM_ENDPOINT_FIELDS, pmtype, ())
    endpoints = BusKey[]
    for (domain, field) in fields
        bus = _record_get(record, field, nothing)
        bus isa Integer || return BusKey[]
        bus > 0 && push!(endpoints, BusKey(domain, Int(bus)))
    end
    return endpoints
end

function _pm_component_active(pmtype::String, record::AbstractDict)
    field = get(_PM_STATUS_FIELDS, pmtype, nothing)
    field === nothing && return nothing
    status = _record_get(record, field, nothing)
    status isa Bool && return status
    status isa Number && return !iszero(status)
    return nothing
end

function _pm_transformer(record::AbstractDict)
    value = _record_get(record, "transformer", false)
    value isa Bool && return value
    value isa Number && return !iszero(value)
    return false
end

# End PowerModels adapter.

function _expected_node_bus(row::DiagramTerminal)
    row.bus == 0 && return (:ground, 0)
    return (row.domain == 1 ? :ac : :dc, row.bus)
end

function _validate_powerflow_network(network, powerflow)
    topology_nodes = Set(network.topology.connections.node)
    mapped_nodes = Set(keys(powerflow.nodes2bus))
    topology_nodes == mapped_nodes || throw(
        ArgumentError(
        "power-flow result does not belong to the supplied network: " *
        "node mappings differ",
    ),
    )
    for row in network.topology.connections
        mapped = powerflow.nodes2bus[row.node]
        mapped == _expected_node_bus(row) || throw(
            ArgumentError(
            "power-flow result does not belong to the supplied network: " *
            "node :$(row.node) maps to $mapped, expected $(_expected_node_bus(row))",
        ),
        )
    end

    mapped_elements = Set(keys(powerflow.elem2comp))
    network_elements = Set(keys(network.elements))
    if powerflow.data === nothing
        isempty(mapped_elements) || throw(
            ArgumentError(
            "power-flow result has component mappings but no PowerModels input data",
        ),
        )
    elseif mapped_elements != network_elements
        throw(
            ArgumentError(
            "power-flow result does not belong to the supplied network: " *
            "element mappings differ",
        ),
        )
    end
    return nothing
end

function _powerflow_records(powerflow, element::Symbol)
    mapping = _pm_mapping(powerflow.elem2comp, element)
    mapping === nothing && return (nothing, nothing, nothing, nothing)
    powerflow.data === nothing && throw(
        ArgumentError(
        "power-flow mapping for :$element cannot be resolved because data is nothing",
    ),
    )
    input = _pm_indexed_record(powerflow.data, mapping.pmtype, mapping.pmindex)
    input isa AbstractDict || throw(
        ArgumentError(
        "power-flow mapping for :$element expects PowerModels component " *
        "$(mapping.pmtype)[$(mapping.pmindex)], but its input record is missing",
    ),
    )
    solution = _pm_indexed_record(
        _pm_solution_root(powerflow.result),
        mapping.pmtype,
        mapping.pmindex
    )
    return (mapping.pmtype, mapping.pmindex, input, solution)
end

function _topology_role(element, connected_buses, grounded_sides)
    model = element.element_model
    PowerImpedance.is_converter(element) && return :converter
    PowerImpedance.is_source(element) && return :source
    (PowerImpedance.is_generator(element) ||
     PowerImpedance.is_inductionmachine(element)) && return :machine
    model isa PowerImpedance.Transformer && return :transformer
    if model isa PowerImpedance.Impedance
        isempty(grounded_sides) || return :shunt
        return all(connection -> connection.bus.domain == 2, connected_buses) ?
               :dc_branch : :ac_branch
    end
    if model isa PowerImpedance.Transmission_line
        return all(connection -> connection.bus.domain == 2, connected_buses) ?
               :dc_branch : :ac_branch
    end
    return :generic
end

function _refined_role(role::Symbol, element, pmtype, input)
    pmtype === nothing && return role
    pmtype == "convdc" && return :converter
    pmtype == "branchdc" && return :dc_branch
    pmtype == "branch" && return _pm_transformer(input) ? :transformer : :ac_branch
    pmtype == "load" && return :load
    pmtype == "shunt" && return :shunt
    pmtype == "im" && return :machine
    if pmtype in ("gen", "gendc")
        return PowerImpedance.is_source(element) ? :source : :machine
    end
    return role
end

function _validate_component_endpoints(component::ComponentDescriptor)
    component.input === nothing && return nothing
    component.pmtype in ("branch", "branchdc", "convdc", "load", "shunt") ||
        return nothing
    expected = Set(connection.bus for connection in component.connected_buses)
    mapped = Set(_pm_component_endpoints(component.pmtype, component.input))
    expected == mapped || throw(
        ArgumentError(
        "power-flow mapping for :$(component.key.element) connects $(collect(mapped)); " *
        "the supplied network connects $(collect(expected))",
    ),
    )
    return nothing
end

function _bus_descriptors(bus_order, bus_nodes, powerflow)
    descriptors = BusDescriptor[]
    solution_root = powerflow === nothing ? nothing : _pm_solution_root(powerflow.result)
    for key in bus_order
        input = powerflow === nothing || powerflow.data === nothing ?
                nothing : _pm_bus_record(powerflow.data, key)
        if powerflow !== nothing && powerflow.data !== nothing && input === nothing
            throw(
                ArgumentError(
                "power-flow input is missing $(_pm_bus_table(key))[$(key.bus)] " *
                "for diagram bus $(key)",
            ),
            )
        end
        solution = solution_root === nothing ? nothing : _pm_bus_record(solution_root, key)
        push!(descriptors, BusDescriptor(key, Tuple(bus_nodes[key]), input, solution))
    end
    return descriptors
end

function _projection_parts(network, powerflow)
    powerflow === nothing || _validate_powerflow_network(network, powerflow)
    rows = network.topology.connections
    bus_order = BusKey[]
    bus_nodes = Dict{BusKey, Vector{Symbol}}()
    for row in rows
        row.bus == 0 && continue
        key = BusKey(row.domain, row.bus)
        if !haskey(bus_nodes, key)
            bus_nodes[key] = Symbol[]
            push!(bus_order, key)
        end
        row.node in bus_nodes[key] || push!(bus_nodes[key], row.node)
    end

    buses = _bus_descriptors(bus_order, bus_nodes, powerflow)
    components = ComponentDescriptor[]
    incidences = IncidenceDescriptor[]
    diagnostics = DiagramDiagnostic[]

    for (element_name, element) in pairs(network.elements)
        terminal_rows = DiagramTerminal[row for row in rows if row.element == element_name]
        if isempty(terminal_rows)
            push!(
                diagnostics,
                DiagramDiagnostic(
                    :disconnected,
                    ComponentKey(element_name),
                    "element :$element_name has no topology rows and was omitted"
                )
            )
            continue
        end

        connected = SideBus[]
        seen = Set{Tuple{Int, Int, Int}}()
        grounded_sides = Int[]
        for row in terminal_rows
            if row.bus == 0
                row.side in grounded_sides || push!(grounded_sides, row.side)
                continue
            end
            incidence_key = (row.side, row.domain, row.bus)
            incidence_key in seen && continue
            push!(seen, incidence_key)
            side_bus = SideBus(row.side, BusKey(row.domain, row.bus))
            push!(connected, side_bus)
            push!(incidences, IncidenceDescriptor(ComponentKey(element_name), row.side, side_bus.bus))
        end

        pmtype, pmindex, input, solution = powerflow === nothing ?
                                           (nothing, nothing, nothing, nothing) :
                                           _powerflow_records(powerflow, element_name)
        role = _refined_role(
            _topology_role(element, connected, grounded_sides),
            element,
            pmtype,
            input
        )
        active = input === nothing ? nothing : _pm_component_active(pmtype, input)
        component = ComponentDescriptor(
            ComponentKey(element_name),
            role,
            Tuple(connected),
            Tuple(grounded_sides),
            Tuple(terminal_rows),
            typeof(element.element_model),
            pmtype,
            pmindex,
            input,
            solution,
            active
        )
        _validate_component_endpoints(component)
        push!(components, component)
    end
    return buses, components, incidences, diagnostics
end

function _layout_graph(buses, components, incidences)
    vertex_keys = DiagramEntityKey[
        (bus.key for bus in buses)...,
        (component.key for component in components)...
    ]
    key_to_vertex = Dict{DiagramEntityKey, Int}(
        key => index for (index, key) in enumerate(vertex_keys)
    )
    graph = Graphs.SimpleGraph(length(vertex_keys))
    for incidence in incidences
        Graphs.add_edge!(
            graph,
            key_to_vertex[incidence.bus],
            key_to_vertex[incidence.component]
        )
    end
    return graph, vertex_keys, key_to_vertex
end

function project_diagram(network::NetworkBuilder.NetworkState, powerflow = nothing)
    powerflow === nothing || powerflow isa PowerFlowResult ||
        throw(
            ArgumentError("diagram enrichment requires a PowerFlowResult"),
        )
    buses, components, incidences, diagnostics = _projection_parts(network, powerflow)
    graph, vertex_keys, key_to_vertex = _layout_graph(buses, components, incidences)
    operating_point = powerflow === nothing ? nothing : powerflow.operating_point
    powerflow_diagnostics = powerflow === nothing ? nothing : powerflow.diagnostics
    return DiagramModel(
        buses,
        components,
        incidences,
        graph,
        vertex_keys,
        key_to_vertex,
        diagnostics,
        operating_point,
        powerflow_diagnostics
    )
end
