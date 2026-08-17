function _deserialize_extension end

_float_type(::Val{:Float16}) = Float16
_float_type(::Val{:Float32}) = Float32
_float_type(::Val{:Float64}) = Float64
_float_type(::Val{:BigFloat}) = BigFloat

function _decode_float(type_name::AbstractString, value::AbstractDict)
    T = _float_type(Val(Symbol(type_name)))
    if haskey(value, "special")
        special = value["special"]
        special == "Inf" && return T(Inf)
        special == "-Inf" && return T(-Inf)
        special == "NaN" && return T(NaN)
        throw(ArgumentError("unknown special floating-point value '$special'"))
    end
    haskey(value, "value") || throw(ArgumentError(
        "serialized $type_name requires a value field",
    ))
    payload = value["value"]
    return T === BigFloat ? parse(BigFloat, String(payload)) : convert(T, payload)
end

function _deserialize_value(value)
    value isa AbstractVector && return [_deserialize_value(item) for item in value]
    value isa AbstractDict || return value
    if haskey(value, "__type__")
        marker = String(value["__type__"])
        marker in ("Float16", "Float32", "Float64", "BigFloat") &&
            return _decode_float(marker, value)
        marker == "Symbol" && return Symbol(_required(value, "value", marker))
        marker == "Complex" && return complex(
            _deserialize_value(_required(value, "re", marker)),
            _deserialize_value(_required(value, "im", marker))
        )
        if marker == "Measurement"
            applicable(_deserialize_extension, Val(:Measurement), value) || throw(
                ArgumentError("deserializing Measurement values requires loading Measurements.jl"),
            )
            return _deserialize_extension(Val(:Measurement), value)
        end
        throw(ArgumentError("unsupported serialized scalar tag '$marker'"))
    end
    haskey(value, "type") && return _decode_object(Val(Symbol(value["type"])), value)
    return Dict(String(key) => _deserialize_value(item) for (key, item) in value)
end

function _required(value, name::AbstractString, owner)
    haskey(value, name) || throw(ArgumentError(
        "serialized $owner requires a '$name' field",
    ))
    return value[name]
end

function _field(value, name::AbstractString)
    _deserialize_value(
        _required(value, name, get(value, "type", "object")),
    )
end
function _optional(value, name::AbstractString)
    haskey(value, name) ? _deserialize_value(value[name]) : nothing
end

function _decode_object(::Val{:Material}, value)
    return Material(
        _field(value, "rho"), _field(value, "eps_r"), _field(value, "mu_r"),
        _field(value, "T0"), _field(value, "alpha")
    )
end

function _decode_object(::Val{:NominalData}, value)
    return NominalData(;
        designation_code = _optional(value, "designation_code"),
        U0 = _optional(value, "U0"),
        U = _optional(value, "U"),
        conductor_cross_section = _optional(value, "conductor_cross_section"),
        screen_cross_section = _optional(value, "screen_cross_section"),
        armor_cross_section = _optional(value, "armor_cross_section"),
        resistance = _optional(value, "resistance"),
        capacitance = _optional(value, "capacitance"),
        inductance = _optional(value, "inductance")
    )
end

function _decode_object(::Val{:CircStrands}, value)
    return CircStrands(
        _field(value, "r_in"), _field(value, "r_ex"),
        _field(value, "radius_wire"), _field(value, "num_wires"),
        _field(value, "lay_ratio"), _field(value, "material");
        lay_direction = _field(value, "lay_direction")
    )
end

function _decode_object(::Val{:RectStrands}, value)
    return RectStrands(
        _field(value, "r_in"), _field(value, "r_ex"),
        _field(value, "thickness"), _field(value, "width"),
        _field(value, "num_wires"), _field(value, "lay_ratio"),
        _field(value, "material"); lay_direction = _field(value, "lay_direction")
    )
end

function _decode_object(::Val{:Strip}, value)
    return Strip(
        _field(value, "r_in"), _field(value, "r_ex"), _field(value, "width"),
        _field(value, "lay_ratio"), _field(value, "material");
        lay_direction = _field(value, "lay_direction")
    )
end

function _decode_object(::Val{:Tubular}, value)
    Tubular(
        _field(value, "r_in"), _field(value, "r_ex"), _field(value, "material")
    )
end
function _decode_object(::Val{:Insulator}, value)
    Insulator(
        _field(value, "r_in"), _field(value, "r_ex"), _field(value, "material")
    )
end
function _decode_object(::Val{:Semicon}, value)
    Semicon(
        _field(value, "r_in"), _field(value, "r_ex"), _field(value, "material")
    )
end

function _decode_group(::Type{G}, layers; kwargs...) where {G}
    isempty(layers) && throw(ArgumentError("serialized $G requires at least one layer"))
    group = G(first(layers); kwargs...)
    for layer in Iterators.drop(layers, 1)
        add!(group, layer)
    end
    return validate(group)
end

function _decode_object(::Val{:ConductorGroup}, value)
    return _decode_group(ConductorGroup, _field(value, "layers"))
end
function _decode_object(::Val{:InsulatorGroup}, value)
    return _decode_group(
        InsulatorGroup, _field(value, "layers");
        reference_frequency = _field(value, "reference_frequency")
    )
end
function _decode_object(::Val{:CableComponent}, value)
    return CableComponent(
        _field(value, "id"), _field(value, "conductor_group"),
        _field(value, "insulator_group")
    )
end
function _decode_object(::Val{:CableDesign}, value)
    return CableDesign(
        _field(value, "cable_id"), _field(value, "components");
        nominal_data = _optional(value, "nominal_data")
    )
end

function _decode_object(::Val{tag}, value) where {tag}
    throw(ArgumentError("unsupported serialized object tag '$tag'"))
end

function _read_document(file_name::AbstractString, expected_schema::AbstractString)
    document = open(file_name, "r") do io
        JSON3.read(io, Dict{String, Any})
    end
    if !haskey(document, "schema") || !haskey(document, "version")
        retired_legacy_json()
    end
    document["schema"] == expected_schema || throw(ArgumentError(
        "expected schema '$expected_schema', found '$(document["schema"])'",
    ))
    document["version"] == JSON_SCHEMA_VERSION || throw(ArgumentError(
        "unsupported $expected_schema schema version $(document["version"]); " *
        "supported version is $JSON_SCHEMA_VERSION",
    ))
    return document
end
