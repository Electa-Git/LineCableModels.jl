const _PSCAD_FREQUENCY_BINDING = "master:Line_FrePhase_Options"
const _PSCAD_CABLE_BINDING = "master:Cable_Coax"
const _PSCAD_GROUND_BINDING = "master:Line_Ground"

function _pscad_attributes!(node::EzXML.Node, attributes)
    for (name, value) in pairs(attributes)
        node[string(name)] = string(value)
    end
    return node
end

function _pscad_paramlist!(
        parent::EzXML.Node,
        parameters;
        name = nothing,
        link = nothing,
        crc = nothing
)
    list = addelement!(parent, "paramlist")
    name === nothing || (list["name"] = string(name))
    link === nothing || (list["link"] = string(link))
    crc === nothing || (list["crc"] = string(crc))
    for (name, value) in parameters
        parameter = addelement!(list, "param")
        parameter["name"] = string(name)
        parameter["value"] = string(value)
    end
    return list
end

function _pscad_params!(parent::EzXML.Node, parameters)
    return _pscad_paramlist!(parent, parameters; name = "", link = -1, crc = -1)
end

function _pscad_user!(
        parent::EzXML.Node,
        id::Integer,
        binding::AbstractString,
        parameters;
        x::Integer = 0,
        y::Integer = 0
)
    user = addelement!(parent, "User")
    _pscad_attributes!(user,
        (
            id = id,
            name = binding,
            classid = "UserCmp",
            x = x,
            y = y,
            w = 0,
            h = 0,
            z = -1,
            orient = 0,
            defn = binding,
            link = -1,
            q = 4,
            disable = "false"
        ))
    _pscad_params!(user, parameters)
    return user
end

function _pscad_wire!(
        parent::EzXML.Node,
        wire_id::Integer,
        user_id::Integer,
        kind::AbstractString,
        name::AbstractString,
        definition::AbstractString,
        user_binding::AbstractString,
        parameters;
        x::Integer,
        y::Integer,
        width::Integer,
        height::Integer,
        crc = nothing
)
    wire = addelement!(parent, "Wire")
    attributes = (
        id = wire_id,
        name,
        classid = kind,
        x,
        y,
        w = width,
        h = height,
        orient = 0,
        defn = definition,
        recv = -1,
        send = -1,
        back = -1,
        disable = "false"
    )
    _pscad_attributes!(wire, attributes)
    crc === nothing || (wire["crc"] = string(crc))
    for (vertex_x, vertex_y) in ((0, 0), (0, 18), (54, 54), (54, 72))
        vertex = addelement!(wire, "vertex")
        _pscad_attributes!(vertex, (x = vertex_x, y = vertex_y))
    end
    _pscad_user!(wire, user_id, user_binding, parameters)
    return wire
end

function _pscad_canvas_params!(canvas::EzXML.Node; user::Bool = false)
    parameters = Pair{String, String}[
        "show_grid" => "0",
        "size" => "0",
        "orient" => "1",
        "show_border" => "0",
        "monitor_bus_voltage" => "0",
        "show_signal" => "0",
        "show_virtual" => "0",
        "show_sequence" => "0",
        "auto_sequence" => "1",
        "bus_expand_x" => "8",
        "bus_expand_y" => "8",
        "bus_length" => "4"
    ]
    user && append!(parameters, [
        "show_terminals" => "0",
        "virtual_filter" => "",
        "animation_freq" => "500"
    ])
    return _pscad_paramlist!(canvas, parameters)
end

function _pscad_value(value; sigdigits::Integer = 6, minimum = -Inf, maximum = Inf)
    scalar = clamp(round(nominal(value); sigdigits), minimum, maximum)
    iszero(scalar) && return "0.0"
    return string(scalar)
end
