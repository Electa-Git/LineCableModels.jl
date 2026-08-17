"""
$(TYPEDSIGNATURES)

Save a cable library as versioned JSON or trusted Julia serialization (`.jls`).
JLS input must come from a trusted source and use matching package types.
"""
function save(
        library::CablesLibrary;
        file_name::String = "cables_library.json"
)
    extension = lowercase(splitext(file_name)[2])
    if extension == ".jls"
        Serialization.serialize(file_name, library.data)
        return abspath(file_name)
    end
    path = _json_path(file_name)
    open(path, "w") do io
        JSON3.pretty(io, _json_document(library); allow_inf = true)
    end
    return abspath(path)
end

"""
$(TYPEDSIGNATURES)

Atomically replace a cable library from supported JSON or trusted JLS data.
"""
function load!(
        library::CablesLibrary;
        file_name::String = "cables_library.json"
)
    isfile(file_name) || throw(ArgumentError(
        "cables library file not found: '$(_display_path(file_name))'",
    ))
    extension = lowercase(splitext(file_name)[2])
    candidate = if extension == ".jls"
        _trusted_cable_data(Serialization.deserialize(file_name))
    elseif extension == ".json"
        document = _read_document(file_name, CABLES_SCHEMA)
        _decoded_cable_data(_deserialize_value(document["cables"]))
    else
        throw(ArgumentError("CablesLibrary loading requires a .json or .jls file"))
    end
    library.data = candidate
    return library
end

function _decoded_cable_data(decoded)
    decoded isa AbstractDict || throw(ArgumentError(
        "the cables field must be a JSON object",
    ))
    candidate = Dict{String, CableDesign}()
    for (name, design) in decoded
        design isa CableDesign || throw(ArgumentError(
            "cable '$name' decoded as $(typeof(design)), not CableDesign",
        ))
        String(name) == design.cable_id || throw(ArgumentError(
            "cable key '$name' differs from cable_id '$(design.cable_id)'",
        ))
        candidate[String(name)] = validate(design)
    end
    return candidate
end

function _trusted_cable_data(decoded)
    decoded isa AbstractDict || throw(ArgumentError(
        "trusted JLS cable data must be a dictionary",
    ))
    return _decoded_cable_data(decoded)
end
