"""
$(TYPEDSIGNATURES)

Save a material library using the current versioned JSON schema.
"""
function save(
        library::MaterialsLibrary;
        file_name::String = "materials_library.json"
)
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

Atomically replace a material library from a supported versioned JSON file.
The original library remains unchanged if parsing or validation fails.
"""
function load!(
        library::MaterialsLibrary;
        file_name::String = "materials_library.json"
)
    isfile(file_name) || throw(ArgumentError(
        "materials library file not found: '$(_display_path(file_name))'",
    ))
    lowercase(splitext(file_name)[2]) == ".json" || throw(ArgumentError(
        "MaterialsLibrary loading requires a .json file",
    ))
    document = _read_document(file_name, MATERIALS_SCHEMA)
    candidate = Dict(
        name => validate(material) for (name, material) in _document_materials(document)
    )
    library.data = candidate
    return library
end
