const JSON_SCHEMA_VERSION = 2
const MATERIALS_SCHEMA = "LineCableModels.MaterialsLibrary"
const CABLES_SCHEMA = "LineCableModels.CablesLibrary"

_type_tag(::Material) = "Material"
_type_tag(::NominalData) = "NominalData"
_type_tag(::CircStrands) = "CircStrands"
_type_tag(::RectStrands) = "RectStrands"
_type_tag(::Strip) = "Strip"
_type_tag(::Tubular) = "Tubular"
_type_tag(::Insulator) = "Insulator"
_type_tag(::Semicon) = "Semicon"
_type_tag(::ConductorGroup) = "ConductorGroup"
_type_tag(::InsulatorGroup) = "InsulatorGroup"
_type_tag(::CableComponent) = "CableComponent"
_type_tag(::CableDesign) = "CableDesign"

_scalar_tag(::Type{Float16}) = "Float16"
_scalar_tag(::Type{Float32}) = "Float32"
_scalar_tag(::Type{Float64}) = "Float64"
_scalar_tag(::Type{BigFloat}) = "BigFloat"

"""
$(TYPEDSIGNATURES)

Encode a supported scalar, collection, or model value for the versioned JSON
schema. Package extensions add methods for their owned scalar types.
"""
function serialize_value(value::AbstractFloat)
    tag = _scalar_tag(typeof(value))
    payload = value isa BigFloat ? string(value) : value
    if isfinite(value)
        return Dict("__type__" => tag, "value" => payload)
    end
    text = isnan(value) ? "NaN" : signbit(value) ? "-Inf" : "Inf"
    return Dict("__type__" => tag, "special" => text)
end

serialize_value(value::Integer) = value
serialize_value(value::Union{Nothing, String, Bool}) = value
serialize_value(value::Symbol) = Dict("__type__" => "Symbol", "value" => string(value))
function serialize_value(value::Complex)
    Dict(
        "__type__" => "Complex",
        "re" => serialize_value(real(value)),
        "im" => serialize_value(imag(value))
    )
end
function serialize_value(value::AbstractDict)
    Dict(string(key) => serialize_value(item) for (key, item) in value)
end
function serialize_value(value::Union{AbstractVector, Tuple})
    [serialize_value(item) for item in value]
end
serialize_value(value) = _serialize_object(value)

function _object(tag::AbstractString; fields...)
    result = Dict{String, Any}("type" => String(tag))
    for (name, value) in pairs(fields)
        result[string(name)] = serialize_value(value)
    end
    return result
end

function _serialize_object(value::Material)
    _object(
        _type_tag(value);
        rho = value.rho,
        eps_r = value.eps_r,
        mu_r = value.mu_r,
        T0 = value.T0,
        alpha = value.alpha
    )
end

function _serialize_object(value::NominalData)
    _object(
        _type_tag(value);
        designation_code = value.designation_code,
        U0 = value.U0,
        U = value.U,
        conductor_cross_section = value.conductor_cross_section,
        screen_cross_section = value.screen_cross_section,
        armor_cross_section = value.armor_cross_section,
        resistance = value.resistance,
        capacitance = value.capacitance,
        inductance = value.inductance
    )
end

function _serialize_object(value::CircStrands)
    _object(
        _type_tag(value);
        r_in = value.r_in,
        r_ex = value.r_ex,
        radius_wire = value.radius_wire,
        num_wires = value.num_wires,
        lay_ratio = value.lay_ratio,
        lay_direction = value.lay_direction,
        material = value.material_props
    )
end

function _serialize_object(value::RectStrands)
    _object(
        _type_tag(value);
        r_in = value.r_in,
        r_ex = value.r_ex,
        thickness = value.thickness,
        width = value.width,
        num_wires = value.num_wires,
        lay_ratio = value.lay_ratio,
        lay_direction = value.lay_direction,
        material = value.material_props
    )
end

function _serialize_object(value::Strip)
    _object(
        _type_tag(value);
        r_in = value.r_in,
        r_ex = value.r_ex,
        width = value.width,
        lay_ratio = value.lay_ratio,
        lay_direction = value.lay_direction,
        material = value.material_props
    )
end

function _serialize_object(value::Tubular)
    return _object(
        _type_tag(value);
        r_in = value.r_in,
        r_ex = value.r_ex,
        material = value.material_props
    )
end

function _serialize_object(value::Insulator)
    return _object(
        _type_tag(value);
        r_in = value.r_in,
        r_ex = value.r_ex,
        material = value.material_props
    )
end

function _serialize_object(value::Semicon)
    return _object(
        _type_tag(value);
        r_in = value.r_in,
        r_ex = value.r_ex,
        material = value.material_props
    )
end

_serialize_object(value::ConductorGroup) = _object(_type_tag(value); layers = value.layers)
function _serialize_object(value::InsulatorGroup)
    _object(
        _type_tag(value);
        layers = value.layers,
        reference_frequency = value.reference_frequency
    )
end
function _serialize_object(value::CableComponent)
    _object(
        _type_tag(value);
        id = value.id,
        conductor_group = value.conductor_group,
        insulator_group = value.insulator_group
    )
end
function _serialize_object(value::CableDesign)
    _object(
        _type_tag(value);
        cable_id = value.cable_id,
        nominal_data = value.nominal_data,
        components = value.components
    )
end

function _json_document(library::MaterialsLibrary)
    return Dict(
        "schema" => MATERIALS_SCHEMA,
        "version" => JSON_SCHEMA_VERSION,
        "materials" => serialize_value(library.data)
    )
end

function _json_document(library::CablesLibrary)
    return Dict(
        "schema" => CABLES_SCHEMA,
        "version" => JSON_SCHEMA_VERSION,
        "cables" => serialize_value(library.data)
    )
end
