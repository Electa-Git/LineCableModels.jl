const _PSCAD_FREQUENCY_BINDING = "master:Line_FrePhase_Options"
const _PSCAD_CABLE_BINDING = "master:Cable_Coax"
const _PSCAD_GROUND_BINDING = "master:Line_Ground"

function _pscad_attributes!(node::EzXML.Node, attributes)
    for (name, value) in pairs(attributes)
        node[string(name)] = string(value)
    end
    return node
end

function _pscad_params!(parent::EzXML.Node, parameters)
    list = addelement!(parent, "paramlist")
    list["name"] = ""
    list["link"] = "-1"
    list["crc"] = "-1"
    for (name, value) in parameters
        parameter = addelement!(list, "param")
        parameter["name"] = string(name)
        parameter["value"] = string(value)
    end
    return list
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

function _pscad_value(value; sigdigits::Integer = 6, minimum = -Inf, maximum = Inf)
    scalar = clamp(round(nominal(value); sigdigits), minimum, maximum)
    iszero(scalar) && return "0.0"
    return string(scalar)
end
