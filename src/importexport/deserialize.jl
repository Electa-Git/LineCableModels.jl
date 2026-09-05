"Decode an extension-owned tagged value from the v1 JSON format."
function deserialize_extension end

_float_type(::Val{:Float16}) = Float16
_float_type(::Val{:Float32}) = Float32
_float_type(::Val{:Float64}) = Float64
_float_type(::Val{:BigFloat}) = BigFloat

function _required(value, name::AbstractString, owner)
    haskey(value, name) || throw(ArgumentError(
        "$owner requires a '$name' field"
    ))
    return value[name]
end

function _field(value, name::AbstractString)
    deserialize_value(
        _required(value, name, get(value, "kind", get(value, "type", "object")))
    )
end
function _optional(value, name::AbstractString)
    haskey(value, name) ? deserialize_value(value[name]) : nothing
end

function _decode_named_tuple(value)
    value isa AbstractVector || throw(ArgumentError(
        "named-tuple data must be an ordered array"
    ))
    names = Tuple(Symbol(_required(entry, "name", "named-tuple entry"))
    for entry in value)
    allunique(names) || throw(ArgumentError("named-tuple names must be unique"))
    values = Tuple(deserialize_value(
                       _required(entry, "value", "named-tuple entry")
                   ) for entry in value)
    return NamedTuple{names}(values)
end

function _decode_float(type_name::AbstractString, value::AbstractDict)
    T = _float_type(Val(Symbol(type_name)))
    if haskey(value, "special")
        special = value["special"]
        special == "Inf" && return T(Inf)
        special == "-Inf" && return T(-Inf)
        special == "NaN" && return T(NaN)
        throw(ArgumentError("unknown special floating-point value '$special'"))
    end
    payload = _required(value, "value", type_name)
    return T === BigFloat ? parse(BigFloat, String(payload)) : convert(T, payload)
end

function _decode_material_record(value)
    raw_kind = deserialize_value(_required(value, "kind", "material"))
    kind = raw_kind isa AbstractString ? Symbol(raw_kind) : raw_kind
    fields = (
        kind,
        _field(value, "rho"),
        _field(value, "eps_r"),
        _field(value, "mu_r"),
        _field(value, "T0"),
        _field(value, "alpha"),
        _field(value, "rho_thermal"),
        _field(value, "theta_max"),
        _field(value, "tan_delta"),
        _field(value, "sigma_solar")
    )
    build = (selected_kind,
        properties...) -> Material(
        selected_kind isa Symbol ? selected_kind : Symbol(selected_kind),
        properties...
    )
    return _decoded_target(Material, build, fields)
end

_decoded_source(value) = value isa Union{AbstractGrid, Gridspace}
function _decoded_target(::Type{Target}, build, values::Tuple) where {Target}
    any(_decoded_source, values) || return build(values...)
    sources = map(values) do value
        _decoded_source(value) ? value : Grid((value,))
    end
    return Gridspace{Target}(build, sources)
end

"Decode a supported scalar, collection, Grid, or v1 declaration."
function deserialize_value(value)
    value isa AbstractVector && return [deserialize_value(item) for item in value]
    value isa AbstractDict || return value
    if haskey(value, "__type__")
        marker = String(value["__type__"])
        marker in ("Float16", "Float32", "Float64", "BigFloat") &&
            return _decode_float(marker, value)
        marker == "Symbol" && return Symbol(_required(value, "value", marker))
        marker == "Complex" && return complex(
            deserialize_value(_required(value, "re", marker)),
            deserialize_value(_required(value, "im", marker))
        )
        if marker == "Measurement"
            applicable(deserialize_extension, Val(:Measurement), value) || throw(
                ArgumentError("deserialising Measurement values requires Measurements.jl")
            )
            return deserialize_extension(Val(:Measurement), value)
        end
        throw(ArgumentError("unsupported serialised scalar tag '$marker'"))
    end
    if haskey(value, "grid")
        values = deserialize_value(value["grid"])
        haskey(value, "rel") && return Grid(values, deserialize_value(value["rel"]))
        haskey(value, "abs") && return Grid(
            values,
            AbsoluteError(deserialize_value(value["abs"]))
        )
        return Grid(values)
    end
    if get(value, "type", nothing) == "material"
        return _decode_material_record(_required(value, "value", "material"))
    end
    haskey(value, "type") && throw(ArgumentError(
        "unsupported serialised object type '$(value["type"])'"
    ))
    haskey(value, "kind") && return _decode_node(Val(Symbol(value["kind"])), value)
    return Dict(String(key) => deserialize_value(item) for (key, item) in value)
end

_decode_node(::Val{:disk}, value) = _decoded_target(Disk, Disk, (_field(value, "r"),))
function _decode_node(::Val{:rectangle}, value)
    _decoded_target(
        Rectangle,
        Rectangle,
        (_field(value, "w"), _field(value, "h"))
    )
end
function _decode_node(::Val{:ellipse}, value)
    _decoded_target(
        Ellipse,
        Ellipse,
        (_field(value, "a"), _field(value, "b"))
    )
end
function _decode_node(::Val{:sector}, value)
    _decoded_target(
        Sector,
        Sector,
        (
            _field(value, "span"),
            _field(value, "r_base"),
            _field(value, "r_back"),
            _field(value, "fillet")
        )
    )
end
function _decode_node(::Val{:annulus}, value)
    _decoded_target(
        Annulus,
        Annulus,
        (_field(value, "ri"), _field(value, "ro"))
    )
end
_decode_node(::Val{:shell}, value) = _decoded_target(Shell, Shell, (_field(value, "t"),))
function _decode_node(::Val{:polygon}, value)
    _decoded_target(Polygon, Polygon, (_field(value, "points"),))
end
function _decode_node(::Val{:pose2}, value)
    _decoded_target(
        Pose2,
        Pose2,
        (_field(value, "x"), _field(value, "y"), _field(value, "φ"))
    )
end
function _decode_node(::Val{:earth_layer}, value)
    EarthLayer(
        _field(value, "rho"),
        _field(value, "eps_r"),
        _field(value, "mu_r"),
        _field(value, "thickness")
    )
end
function _decode_node(::Val{:earth_model}, value)
    raw_layers = _required(value, "layers", "earth_model")
    raw_layers isa AbstractVector || throw(ArgumentError(
        "earth_model layers must be an array"
    ))
    layers = EarthLayer[deserialize_value(layer) for layer in raw_layers]
    isempty(layers) && throw(ArgumentError("earth_model layers cannot be empty"))
    T = promote_type(map(eltype, layers)...)
    converted = Tuple(convert(EarthLayer{T}, layer) for layer in layers)
    length(converted) >= 2 || throw(ArgumentError(
        "earth_model requires air and at least one earth layer"
    ))
    return build(
        EarthModel,
        Base.tail(converted);
        vertical_layers = Bool(_required(value, "vertical_layers", "earth_model")),
        air_layer = first(converted)
    )
end

function _decode_node(::Val{:ring}, value)
    values = (
        _field(value, "n"), _field(value, "r"),
        _field(value, "φ0"), _field(value, "span"),
        haskey(value, "gap_frac") ? _field(value, "gap_frac") : 0
    )
    build = (n, r, φ0, span, gap_frac) -> Ring(n; r, φ0, span, gap_frac)
    return _decoded_target(Ring, build, values)
end
function _decode_node(::Val{:polar}, value)
    values = (
        _field(value, "nr"), _field(value, "nφ"),
        _field(value, "r0"), _field(value, "dr"),
        _field(value, "φ0"), _field(value, "span")
    )
    build = (nr, nφ, r0, dr, φ0, span) -> Polar(; nr, nφ, r0, dr, φ0, span)
    return _decoded_target(Polar, build, values)
end
function _decode_node(::Val{:fill}, value)
    values = (
        _field(value, "r"), _field(value, "φ"),
        _field(value, "φ0"), _field(value, "span")
    )
    build = (r, φ, φ0, span) -> Fill(; r, φ, φ0, span)
    return _decoded_target(Fill, build, values)
end
function _decode_node(::Val{:lattice}, value)
    values = (
        _field(value, "nx"), _field(value, "ny"),
        _field(value, "dx"), _field(value, "dy")
    )
    build = (nx, ny, dx, dy) -> Lattice(; nx, ny, dx, dy)
    return _decoded_target(Lattice, build, values)
end
function _decode_node(::Val{:fill_factor}, value)
    _decoded_target(
        FillFactor, FillFactor, (_field(value, "η"),)
    )
end
_decode_node(::Val{:capacity}, value) = capacity()
function _decode_node(::Val{:lay_ratio}, value)
    _decoded_target(LayRatio, LayRatio, (_field(value, "q"),))
end
_decode_node(::Val{:pitch}, value) = _decoded_target(Pitch, Pitch, (_field(value, "p"),))
function _decode_node(::Val{:lay_angle}, value)
    _decoded_target(LayAngle, LayAngle, (_field(value, "α"),))
end
function _decode_node(::Val{:helix}, value)
    values = (
        _field(value, "lay"), _field(value, "dir"), _field(value, "φ0")
    )
    build = (lay, dir, φ0) -> Helix(lay; dir, φ0)
    return _decoded_target(Helix, build, values)
end

function _material_reference(value, materials)
    value isa AbstractString && return get(materials, String(value)) do
        throw(KeyError(String(value)))
    end
    decoded = deserialize_value(value)
    decoded isa Union{Material, Gridspace{Material}} || throw(ArgumentError(
        "region material must decode as Material"
    ))
    return decoded
end

function _decode_part(value, materials)
    value isa AbstractDict || throw(ArgumentError("physical nodes must be objects"))
    kind = Symbol(_required(value, "kind", "physical node"))
    return _decode_part(Val(kind), value, materials)
end

function _decode_part(::Val{:region}, value, materials)
    values = (
        Symbol(_required(value, "tag", "region")),
        deserialize_value(_required(value, "primitive", "region")),
        _material_reference(_required(value, "material", "region"), materials)
    )
    return _decoded_target(Region, Region, values)
end
function _decode_part(::Val{:stack}, value, materials)
    items = _required(value, "items", "stack")
    items isa AbstractVector && !isempty(items) || throw(ArgumentError(
        "stack items must be a nonempty array"
    ))
    decoded = Tuple(_decode_part(item, materials) for item in items)
    return _decoded_target(Stack, Stack, decoded)
end
function _decode_part(::Val{:group}, value, materials)
    path = deserialize_value(_required(value, "path", "group"))
    path isa AbstractVector && (path = Tuple(path))
    values = (
        Symbol(_required(value, "name", "group")),
        deserialize_value(_required(value, "at", "group")),
        _decode_part(_required(value, "item", "group"), materials),
        deserialize_value(_required(value, "pattern", "group")),
        path,
        deserialize_value(_required(value, "compact", "group")),
        deserialize_value(_required(value, "boundary", "group"))
    )
    return _decoded_target(Group, Group, values)
end
function _decode_part(::Val{:assembly}, value, materials)
    if haskey(value, "members")
        raw_members = _required(value, "members", "assembly")
        raw_members isa AbstractVector && !isempty(raw_members) || throw(
            ArgumentError("explicit assembly members must be a nonempty array")
        )
        members = Tuple(map(raw_members) do member
            member isa AbstractDict || throw(ArgumentError(
                "explicit assembly members must be objects"
            ))
            DataModel.AssemblyMember(
                _decode_part(_required(member, "item", "assembly member"), materials),
                deserialize_value(_required(member, "at", "assembly member"))
            )
        end)
        return Assembly(
            deserialize_value(_required(value, "at", "assembly")),
            members,
            nothing,
            nothing,
            nothing,
            nothing
        )
    end
    raw_names = _required(value, "names", "assembly")
    names = raw_names === nothing ? nothing : Symbol.(raw_names)
    values = (
        deserialize_value(_required(value, "at", "assembly")),
        _decode_part(_required(value, "item", "assembly"), materials),
        deserialize_value(_required(value, "pattern", "assembly")),
        deserialize_value(_required(value, "path", "assembly")),
        deserialize_value(_required(value, "compact", "assembly")),
        names
    )
    return _decoded_target(Assembly, Assembly, values)
end
function _decode_part(::Val{:enclosure}, value, materials)
    raw_fill = _required(value, "fill", "enclosure")
    fill = raw_fill isa AbstractString ? _material_reference(raw_fill, materials) :
           get(raw_fill, "kind", nothing) == "region" ?
           _decode_part(raw_fill, materials) : deserialize_value(raw_fill)
    raw_wall = _required(value, "wall", "enclosure")
    wall = raw_wall === nothing ? nothing : _decode_part(raw_wall, materials)
    values = (
        Symbol(_required(value, "tag", "enclosure")),
        deserialize_value(_required(value, "at", "enclosure")),
        deserialize_value(_required(value, "primitive", "enclosure")),
        _decode_part(_required(value, "item", "enclosure"), materials),
        fill,
        wall
    )
    return _decoded_target(Enclosure, Enclosure, values)
end
function _decode_part(::Val{kind}, value, materials) where {kind}
    throw(ArgumentError("unsupported physical node kind '$kind'"))
end

_decode_node(::Val{:region}, value) = _decode_part(value, Dict{String, Material}())
_decode_node(::Val{:stack}, value) = _decode_part(value, Dict{String, Material}())
_decode_node(::Val{:group}, value) = _decode_part(value, Dict{String, Material}())
_decode_node(::Val{:assembly}, value) = _decode_part(value, Dict{String, Material}())
_decode_node(::Val{:enclosure}, value) = _decode_part(value, Dict{String, Material}())

function _decode_design_resolved(value, materials)
    get(value, "kind", nothing) == "cable_design" || throw(ArgumentError(
        "cable declaration must have kind 'cable_design'"
    ))
    values = (
        String(_required(value, "cable_id", "cable_design")),
        _decode_part(_required(value, "root", "cable_design"), materials),
        haskey(value, "nominal_data") ?
        _decode_named_tuple(_required(value, "nominal_data", "cable_design")) : (;)
    )
    caller = (
        cable_id, root, nominal_data) -> build(
        CableDesign, cable_id, root; nominal_data
    )
    return _decoded_target(CableDesign, caller, values)
end

function _decode_design(value, materials = Dict{String, Material}())
    varying = Tuple(
        (name, material)
    for (name, material) in sort!(collect(materials); by = first)
    if material isa Gridspace{Material}
    )
    isempty(varying) && return _decode_design_resolved(value, materials)
    names = first.(varying)
    sources = last.(varying)
    build = function (selected...)
        resolved = Dict{String, Any}(materials)
        for (name, material) in zip(names, selected)
            resolved[name] = material
        end
        design = _decode_design_resolved(value, resolved)
        design isa CableDesign || throw(ArgumentError(
            "a named-material parameter space cannot be combined with another " *
            "unresolved design field in the same JSON declaration"
        ))
        return design
    end
    return Gridspace{CableDesign}(build, sources)
end
_decode_node(::Val{:cable_design}, value) = _decode_design(value)

function _decode_node(::Val{:line_cable_system}, value)
    designs = CableDesign[deserialize_value(item)
                          for item in _required(value, "designs", "line_cable_system")]
    return build(
        LineCableSystem,
        designs,
        deserialize_value(_required(value, "positions", "line_cable_system"));
        system_id = String(_required(value, "system_id", "line_cable_system")),
        line_length = _field(value, "line_length"),
        connections = deserialize_value(_required(value, "connections", "line_cable_system")),
        environment = _optional(value, "environment")
    )
end

function _decode_node(::Val{:line_parameters_problem}, value)
    system = _field(value, "system")
    system isa LineCableSystem || throw(ArgumentError(
        "line_parameters_problem system must decode as LineCableSystem"
    ))
    earth = _field(value, "earth_props")
    earth isa EarthModel || throw(ArgumentError(
        "line_parameters_problem earth_props must decode as EarthModel"
    ))
    return Engine.LineParametersProblem(
        system;
        temperature = _field(value, "temperature"),
        earth_props = earth,
        frequencies = _field(value, "frequencies"),
        Γ = _optional(value, "Gamma")
    )
end

function _decode_node(::Val{kind}, value) where {kind}
    throw(ArgumentError("unsupported declaration kind '$kind'"))
end

function _read_document(file_name::AbstractString, expected_format::AbstractString)
    document = open(file_name, "r") do io
        #! explicit-imports: off
        JSON3.read(io, Dict{String, Any})
        #! explicit-imports: on
    end
    _required(document, "\$schema", "LineCableModels document") ==
    JSON_SCHEMA_DIALECT || throw(ArgumentError(
        "unsupported JSON Schema dialect"
    ))
    format = _required(document, "format", "LineCableModels document")
    format == expected_format || throw(ArgumentError(
        "expected format '$expected_format', found '$format'"
    ))
    version = _required(document, "version", expected_format)
    version == JSON_SCHEMA_VERSION || throw(ArgumentError(
        "unsupported $expected_format version $version; supported version is " *
        JSON_SCHEMA_VERSION
    ))
    return document
end

function _document_materials(document)
    raw = _required(document, "materials", "LineCableModels document")
    raw isa AbstractDict || throw(ArgumentError("materials must be an object"))
    return Dict{String, Any}(
        String(name) => _decode_material_record(value) for (name, value) in raw
    )
end
