"""
$(TYPEDSIGNATURES)

Save a cable library as versioned JSON or trusted Julia serialisation (`.jls`).
JLS input must come from a trusted source and use matching package types.
"""
function save(
        library::CablesLibrary;
        file_name::String = "cables_library.json"
)
    extension = lowercase(splitext(file_name)[2])
    if extension == ".jls"
        Serialization.serialize(file_name, (
            designs = library.data,
            catalogues = library.catalogues
        ))
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
        materials = _document_materials(document)
        root = _required(document, "root", CABLES_SCHEMA)
        get(root, "kind", nothing) == "cable_library" || throw(ArgumentError(
            "cable document root must have kind 'cable_library'"
        ))
        raw_cables = _required(root, "cables", "cable_library")
        raw_cables isa AbstractDict || throw(ArgumentError(
            "cable_library cables must be an object"
        ))
        _decoded_cable_library(raw_cables, materials)
    else
        throw(ArgumentError("CablesLibrary loading requires a .json or .jls file"))
    end
    library.data = candidate.designs
    library.catalogues = candidate.catalogues
    return library
end

function _decode_catalogue(value)
    value isa AbstractDict || throw(ArgumentError(
        "a cable catalogue record must be an object"
    ))
    names = sort!(Symbol.(collect(keys(value))); by = String)
    entries = Tuple(deserialize_value(value[String(name)]) for name in names)
    return NamedTuple{Tuple(names)}(entries)
end

function _decoded_cable_library(raw_cables, materials)
    designs = Dict{String, CableDesign}()
    catalogues = Dict{String, NamedTuple}()
    for (name, value) in raw_cables
        cable_id = String(name)
        design = _decode_design(value, materials)
        design isa CableDesign || throw(ArgumentError(
            "cable '$cable_id' decoded as $(typeof(design)), not CableDesign"
        ))
        cable_id == design.cable_id || throw(ArgumentError(
            "cable key '$cable_id' differs from cable_id '$(design.cable_id)'"
        ))
        designs[cable_id] = validate(design)
        catalogues[cable_id] = _decode_catalogue(
            _required(value, "catalogue", "cable_design")
        )
    end
    return (; designs, catalogues)
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
    decoded isa NamedTuple && keys(decoded) == (:designs, :catalogues) || throw(
        ArgumentError(
            "trusted JLS cable data must contain designs and catalogues"
        )
    )
    designs = _decoded_cable_data(decoded.designs)
    decoded.catalogues isa AbstractDict || throw(ArgumentError(
        "trusted JLS catalogue data must be a dictionary",
    ))
    catalogues = Dict{String, NamedTuple}()
    for cable_id in keys(designs)
        record = get(decoded.catalogues, cable_id) do
            throw(KeyError(cable_id))
        end
        record isa NamedTuple || throw(ArgumentError(
            "catalogue '$cable_id' must be a NamedTuple"
        ))
        catalogues[cable_id] = record
    end
    return (; designs, catalogues)
end
