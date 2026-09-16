"""
$(TYPEDSIGNATURES)

Save a material library as versioned JSON or trusted Julia serialization (`.jls`).
JLS input must come from a trusted source and use matching package types.
"""
function save(
        library::MaterialsLibrary;
        file_name::String = "materials_library.json"
)
    extension = lowercase(splitext(file_name)[2])
    if extension == ".jls"
        Serialization.serialize(file_name, library.data)
        return abspath(file_name)
    end
    path = _json_path(file_name)
    open(path, "w") do io
        #! explicit-imports: off
        # JSON3 exposes this established writer without a public marker.
        JSON3.pretty(io, _json_document(library); allow_inf = true)
        #! explicit-imports: on
    end
    return abspath(path)
end

"""
$(TYPEDSIGNATURES)

Atomically replace a material library from supported JSON or trusted JLS data.
The original library remains unchanged if parsing or validation fails.
"""
function load!(
        library::MaterialsLibrary;
        file_name::String = "materials_library.json"
)
    isfile(file_name) || throw(ArgumentError(
        "materials library file not found: '$(_display_path(file_name))'",
    ))
    extension = lowercase(splitext(file_name)[2])
    candidate = if extension == ".jls"
        _trusted_material_data(Serialization.deserialize(file_name))
    elseif extension == ".json"
        document = _read_document(file_name, MATERIALS_SCHEMA)
        Dict(
            name => validate(material) for (name, material) in _document_materials(document)
        )
    else
        throw(ArgumentError(
            "MaterialsLibrary loading requires a .json or .jls file",
        ))
    end
    library.data = candidate
    return library
end

function _trusted_material_data(decoded)
    decoded isa AbstractDict || throw(ArgumentError(
        "trusted JLS material data must be a dictionary",
    ))
    candidate = Dict{String, Material}()
    for (name, material) in decoded
        material isa Material || throw(ArgumentError(
            "material '$name' must be Material, not $(typeof(material))",
        ))
        candidate[String(name)] = validate(material)
    end
    return candidate
end
