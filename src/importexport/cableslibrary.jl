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
        Serialization.serialize(file_name, (
            designs = library.data,
            datasheets = library.datasheets
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
    decoded_library = if extension == ".jls"
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
    library.data = decoded_library.designs
    library.datasheets = decoded_library.datasheets
    return library
end

function _decode_datasheet(value)
    value isa AbstractDict || throw(ArgumentError(
        "a cable datasheet record must be an object"
    ))
    names = sort!(Symbol.(collect(keys(value))); by = String)
    entries = Tuple(deserialize_value(value[String(name)]) for name in names)
    return DatasheetInfo(NamedTuple{Tuple(names)}(entries))
end

function _decoded_cable_library(raw_cables, materials)
    designs = Dict{String, CableDesign}()
    datasheets = Dict{String, DatasheetInfo}()
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
        datasheets[cable_id] = _decode_datasheet(
            _required(value, "datasheet", "cable_design")
        )
    end
    return (; designs, datasheets)
end

function _decoded_cable_data(decoded)
    decoded isa AbstractDict || throw(ArgumentError(
        "the cables field must be a JSON object",
    ))
    designs = Dict{String, CableDesign}()
    for (name, design) in decoded
        design isa CableDesign || throw(ArgumentError(
            "cable '$name' decoded as $(typeof(design)), not CableDesign",
        ))
        String(name) == design.cable_id || throw(ArgumentError(
            "cable key '$name' differs from cable_id '$(design.cable_id)'",
        ))
        designs[String(name)] = validate(design)
    end
    return designs
end

function _trusted_cable_data(decoded)
    decoded isa NamedTuple && keys(decoded) == (:designs, :datasheets) || throw(
        ArgumentError(
            "trusted JLS cable data must contain designs and datasheets"
        )
    )
    designs = _decoded_cable_data(decoded.designs)
    decoded.datasheets isa AbstractDict || throw(ArgumentError(
        "trusted JLS datasheet data must be a dictionary",
    ))
    datasheets = Dict{String, DatasheetInfo}()
    for cable_id in keys(designs)
        record = get(decoded.datasheets, cable_id) do
            throw(KeyError(cable_id))
        end
        record isa Union{DatasheetInfo, NamedTuple} || throw(ArgumentError(
            "datasheet '$cable_id' must be DatasheetInfo or a named tuple"
        ))
        datasheets[cable_id] =
            record isa DatasheetInfo ? record : DatasheetInfo(record)
    end
    return (; designs, datasheets)
end
