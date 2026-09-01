const JSON_SCHEMA_VERSION = "1.0.0"
const JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"
const MATERIALS_SCHEMA = "linecablemodels.materials"
const CABLES_SCHEMA = "linecablemodels.cable"

_scalar_tag(::Type{Float16}) = "Float16"
_scalar_tag(::Type{Float32}) = "Float32"
_scalar_tag(::Type{Float64}) = "Float64"
_scalar_tag(::Type{BigFloat}) = "BigFloat"

"""Encode a supported scalar, collection, Grid, or v1 declaration."""
function serialize_value(value::AbstractFloat)
    tag = _scalar_tag(typeof(value))
    payload = value isa BigFloat ? string(value) : value
    if isfinite(value)
        return Dict("__type__" => tag, "value" => payload)
    end
    special = isnan(value) ? "NaN" : signbit(value) ? "-Inf" : "Inf"
    return Dict("__type__" => tag, "special" => special)
end

serialize_value(value::Integer) = value
serialize_value(value::Union{Nothing, String, Bool}) = value
serialize_value(value::Symbol) = Dict("__type__" => "Symbol", "value" => String(value))
function serialize_value(value::Complex)
    return Dict(
        "__type__" => "Complex",
        "re" => serialize_value(real(value)),
        "im" => serialize_value(imag(value))
    )
end
function serialize_value(value::AbstractDict)
    Dict(string(key) => serialize_value(item) for (key, item) in value)
end
function serialize_value(value::NamedTuple)
    Dict(string(key) => serialize_value(item) for (key, item) in pairs(value))
end
function serialize_value(value::Union{AbstractVector, Tuple})
    [serialize_value(item) for item in value]
end
serialize_value(value) = _serialize_object(value)

function _node(kind::AbstractString; fields...)
    Dict{String, Any}(
        "kind" => String(kind),
        (String(name) => serialize_value(value) for (name, value) in pairs(fields))...
    )
end

function _material_record(value::Material)
    return Dict{String, Any}(
        "kind" => String(value.kind),
        "rho" => serialize_value(value.rho),
        "eps_r" => serialize_value(value.eps_r),
        "mu_r" => serialize_value(value.mu_r),
        "T0" => serialize_value(value.T0),
        "alpha" => serialize_value(value.alpha),
        "rho_thermal" => serialize_value(value.rho_thermal),
        "theta_max" => serialize_value(value.theta_max),
        "tan_delta" => serialize_value(value.tan_delta),
        "sigma_solar" => serialize_value(value.sigma_solar)
    )
end

function _serialize_object(value::Material)
    Dict(
        "type" => "material",
        "value" => _material_record(value)
    )
end

_serialize_object(value::Disk) = _node("disk"; r = value.r)
_serialize_object(value::Rectangle) = _node("rectangle"; w = value.w, h = value.h)
_serialize_object(value::Ellipse) = _node("ellipse"; a = value.a, b = value.b)
function _serialize_object(value::Sector)
    return _node("sector"; ri = value.ri, ro = value.ro, φ0 = value.φ0,
        span = value.span)
end
_serialize_object(value::Annulus) = _node("annulus"; ri = value.ri, ro = value.ro)
_serialize_object(value::Shell) = _node("shell"; t = value.t)
_serialize_object(value::Polygon) = _node("polygon"; points = value.points)
function _serialize_object(value::RoundedSector)
    return _node(
        "rounded_sector";
        span = value.span,
        r_base = value.r_base,
        r_back = value.r_back,
        fillet = value.fillet
    )
end
_serialize_object(value::Pose2) = _node("pose2"; x = value.x, y = value.y, φ = value.φ)

function _serialize_object(value::EarthLayer)
    return _node(
        "earth_layer";
        rho = value.rho,
        eps_r = value.eps_r,
        mu_r = value.mu_r,
        thickness = value.thickness
    )
end
function _serialize_object(value::EarthModel)
    return _node(
        "earth_model";
        vertical_layers = value.vertical_layers,
        layers = value.layers
    )
end

function _serialize_object(value::Ring)
    return _node("ring"; n = value.n, r = value.r, φ0 = value.φ0,
        span = value.span, gap_frac = value.gap_frac)
end
function _serialize_object(value::Hexa)
    return _node(
        "hexagonal_course";
        course = value.course,
        φ0 = value.φ0,
        gap_frac = value.gap_frac
    )
end
function _serialize_object(value::Polar)
    return _node("polar"; nr = value.nr, nφ = value.nφ, r0 = value.r0,
        dr = value.dr, φ0 = value.φ0, span = value.span)
end
function _serialize_object(value::Fill)
    _node("fill"; r = value.r, φ = value.φ, φ0 = value.φ0, span = value.span)
end
function _serialize_object(value::Lattice)
    _node("lattice"; nx = value.nx, ny = value.ny, dx = value.dx, dy = value.dy)
end
_serialize_object(value::DiameterFactor) = _node("diameter_factor"; k = value.k)
_serialize_object(value::FillFactor) = _node("fill_factor"; η = value.η)
function _serialize_object(value::TabulatedCompaction)
    _node("tabulated_compaction"; data = value.data)
end
_serialize_object(value::AffineCompaction) = _node("affine_compaction"; map = value.map)
_serialize_object(::typeof(capacity())) = _node("capacity")
_serialize_object(value::LayRatio) = _node("lay_ratio"; q = value.q)
_serialize_object(value::Pitch) = _node("pitch"; p = value.p)
_serialize_object(value::LayAngle) = _node("lay_angle"; α = value.α)
function _serialize_object(value::Helix)
    _node("helix"; lay = value.lay, dir = value.dir, φ0 = value.φ0)
end

function _serialize_part(value::Region, material_name)
    material = material_name === nothing ? serialize_value(value.material) :
               String(material_name(value.material))
    return Dict(
        "kind" => "region",
        "tag" => String(value.tag),
        "primitive" => serialize_value(value.primitive),
        "material" => material
    )
end
function _serialize_part(value::Stack, material_name)
    return Dict(
        "kind" => "stack",
        "items" => [_serialize_part(item, material_name) for item in value.items]
    )
end
function _serialize_part(value::Group, material_name)
    return Dict(
        "kind" => "group",
        "name" => String(value.name),
        "at" => serialize_value(value.at),
        "item" => _serialize_part(value.item, material_name),
        "pattern" => serialize_value(value.pattern),
        "path" => serialize_value(value.path),
        "compact" => serialize_value(value.compact)
    )
end
function _serialize_part(
        value::Assembly{<:Any, <:AbstractCablePart},
        material_name
)
    names = value.names === nothing ? nothing : String.(value.names)
    return Dict(
        "kind" => "assembly",
        "at" => serialize_value(value.at),
        "item" => _serialize_part(value.item, material_name),
        "pattern" => serialize_value(value.pattern),
        "path" => serialize_value(value.path),
        "compact" => serialize_value(value.compact),
        "names" => names
    )
end
function _serialize_part(value::Assembly{<:Any, <:Tuple}, material_name)
    return Dict(
        "kind" => "assembly",
        "at" => serialize_value(value.at),
        "members" => [Dict(
             "at" => serialize_value(member.at),
             "item" => _serialize_part(member.item, material_name)
         ) for member in value.item]
    )
end
function _serialize_part(value::Enclosure, material_name)
    fill = value.fill isa Material ?
           (material_name === nothing ? serialize_value(value.fill) :
            String(material_name(value.fill))) :
           _serialize_part(value.fill, material_name)
    wall = value.wall === nothing ? nothing : _serialize_part(value.wall, material_name)
    return Dict(
        "kind" => "enclosure",
        "tag" => String(value.tag),
        "at" => serialize_value(value.at),
        "primitive" => serialize_value(value.primitive),
        "item" => _serialize_part(value.item, material_name),
        "fill" => fill,
        "wall" => wall
    )
end

_serialize_object(value::AbstractCablePart) = _serialize_part(value, nothing)

function _serialize_design(value::CableDesign, material_name = nothing)
    return Dict(
        "kind" => "cable_design",
        "cable_id" => value.cable_id,
        "nominal_data" => [Dict(
             "name" => String(name),
             "value" => serialize_value(item)
         ) for (name, item) in pairs(value.nominal_data)],
        "root" => _serialize_part(value.root, material_name)
    )
end
_serialize_object(value::CableDesign) = _serialize_design(value)

function _serialize_object(value::LineCableSystem)
    return Dict(
        "kind" => "line_cable_system",
        "system_id" => value.system_id,
        "line_length" => serialize_value(value.line_length),
        "designs" => [_serialize_design(design) for design in value.designs],
        "positions" => serialize_value(value.positions),
        "connections" => serialize_value(value.connections),
        "environment" => serialize_value(value.environment)
    )
end

function _serialize_object(value::Engine.LineParametersProblem)
    return Dict(
        "kind" => "line_parameters_problem",
        "system" => serialize_value(value.system),
        "temperature" => serialize_value(value.temperature),
        "earth_props" => serialize_value(value.earth_props),
        "frequencies" => serialize_value(value.frequencies),
        "Gamma" => serialize_value(value.Γ)
    )
end

function _serialize_object(grid::DeterministicGrid)
    return Dict("grid" => serialize_value(collect(grid.vals)))
end
function _serialize_object(grid::RelativeGrid)
    return Dict(
        "grid" => serialize_value(collect(grid.vals)),
        "rel" => serialize_value(collect(grid.rel_err))
    )
end
function _serialize_object(grid::AbsoluteGrid)
    return Dict(
        "grid" => serialize_value(collect(grid.vals)),
        "abs" => serialize_value(collect(grid.abs_err))
    )
end

function _collect_materials!(materials::Vector{Material}, part::Region)
    any(item -> isequal(item, part.material), materials) || push!(materials, part.material)
    return materials
end
function _collect_materials!(materials::Vector{Material}, part::Stack)
    foreach(item -> _collect_materials!(materials, item), part.items)
    return materials
end
function _collect_materials!(materials::Vector{Material}, part::Group)
    return _collect_materials!(materials, part.item)
end
function _collect_materials!(
        materials::Vector{Material},
        part::Assembly{<:Any, <:AbstractCablePart}
)
    return _collect_materials!(materials, part.item)
end
function _collect_materials!(
        materials::Vector{Material},
        part::Assembly{<:Any, <:Tuple}
)
    foreach(member -> _collect_materials!(materials, member.item), part.item)
    return materials
end
function _collect_materials!(materials::Vector{Material}, part::Enclosure)
    _collect_materials!(materials, part.item)
    part.fill isa Material ?
    (any(item -> isequal(item, part.fill), materials) || push!(materials, part.fill)) :
    _collect_materials!(materials, part.fill)
    part.wall === nothing || _collect_materials!(materials, part.wall)
    return materials
end

function _json_document(library::MaterialsLibrary)
    return Dict(
        "\$schema" => JSON_SCHEMA_DIALECT,
        "format" => MATERIALS_SCHEMA,
        "version" => JSON_SCHEMA_VERSION,
        "materials" => Dict(
            name => _material_record(material)
        for (name, material) in sort!(collect(library.data); by = first)
        )
    )
end

function _json_document(library::CablesLibrary)
    materials = Material[]
    cable_ids = sort!(collect(keys(library.data)))
    for cable_id in cable_ids
        _collect_materials!(materials, library.data[cable_id].root)
    end
    names = Dict{Int, String}(
        index => "$(material.kind)_$index" for (index, material) in enumerate(materials)
    )
    material_name(material) = names[only(findall(item -> isequal(item, material), materials))]
    return Dict(
        "\$schema" => JSON_SCHEMA_DIALECT,
        "format" => CABLES_SCHEMA,
        "version" => JSON_SCHEMA_VERSION,
        "materials" => Dict(
            names[index] => _material_record(material)
        for (index, material) in enumerate(materials)
        ),
        "root" => Dict(
            "kind" => "cable_library",
            "cables" => Dict(
                cable_id => merge(
                    _serialize_design(library.data[cable_id], material_name),
                    Dict("catalogue" => Dict(
                        String(name) => serialize_value(item)
                    for (name, item) in pairs(library.catalogues[cable_id])
                    ))
                )
            for cable_id in cable_ids
            )
        )
    )
end

function _json_document(system::LineCableSystem)
    return Dict(
        "\$schema" => JSON_SCHEMA_DIALECT,
        "format" => CABLES_SCHEMA,
        "version" => JSON_SCHEMA_VERSION,
        "materials" => Dict{String, Any}(),
        "root" => serialize_value(system)
    )
end
