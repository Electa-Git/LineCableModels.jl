const JSON_SCHEMA_VERSION = "1.0.0"
const JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"
const MATERIALS_SCHEMA = "linecablemodels.materials"
const CABLES_SCHEMA = "linecablemodels.cable"
const PROBLEM_SCHEMA = "linecablemodels.line_parameters_problem"

_scalar_tag(::Type{Float16}) = "Float16"
_scalar_tag(::Type{Float32}) = "Float32"
_scalar_tag(::Type{Float64}) = "Float64"
_scalar_tag(::Type{BigFloat}) = "BigFloat"

"""Encode a supported scalar, collection, Grid, or v1 declaration."""
function serialize_value(value::AbstractFloat)
    tag = _scalar_tag(typeof(value))
    isfinite(value) && return Dict("__type__" => tag, "value" => value)
    special = isnan(value) ? "NaN" : signbit(value) ? "-Inf" : "Inf"
    return Dict("__type__" => tag, "special" => special)
end
# A finite BigFloat is encoded as text with its precision.
function serialize_value(value::BigFloat)
    isfinite(value) || return Dict("__type__" => "BigFloat",
        "special" => isnan(value) ? "NaN" : signbit(value) ? "-Inf" : "Inf")
    return Dict("__type__" => "BigFloat", "value" => string(value), "precision" => precision(value))
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

# With `material_name`, a material is encoded by its library name.
function serialize_value(value::AbstractMaterial; material_name = nothing)
    material_name === nothing || return String(material_name(value))
    Dict(
        "type" => "material",
        "value" => _material_record(value)
    )
end

serialize_value(value::Disk) = _node("disk"; r = value.r)
serialize_value(value::Rectangle) = _node("rectangle"; w = value.w, h = value.h)
serialize_value(value::Ellipse) = _node("ellipse"; a = value.a, b = value.b)
function serialize_value(value::Sector)
    return _node(
        "sector";
        span = value.span,
        r_base = value.r_base,
        r_back = value.r_back,
        fillet = value.fillet
    )
end
serialize_value(value::Annulus) = _node("annulus"; ri = value.ri, ro = value.ro)
serialize_value(value::Shell) = _node("shell"; t = value.t)
serialize_value(value::Polygon) = _node("polygon"; points = value.points)
serialize_value(value::Pose2) = _node("pose2"; x = value.x, y = value.y, φ = value.φ)

function serialize_value(value::EarthLayer)
    return _node(
        "earth_layer";
        rho = value.rho,
        eps_r = value.eps_r,
        mu_r = value.mu_r,
        thickness = value.thickness
    )
end
function serialize_value(value::EarthModel)
    return _node(
        "earth_model";
        vertical_layers = value.vertical_layers,
        layers = value.layers
    )
end

function serialize_value(value::Ring)
    return _node("ring"; n = value.n, r = value.r, φ0 = value.φ0,
        span = value.span, gap_frac = value.gap_frac)
end
serialize_value(value::FillFactor) = _node("fill_factor"; η = value.η)
serialize_value(::typeof(capacity())) = _node("capacity")
serialize_value(value::LayRatio) = _node("lay_ratio"; q = value.q)
serialize_value(value::Pitch) = _node("pitch"; p = value.p)
serialize_value(value::LayAngle) = _node("lay_angle"; α = value.α)
function serialize_value(value::Helix)
    _node("helix"; lay = value.lay, dir = value.dir, φ0 = value.φ0)
end

function serialize_value(value::Region; material_name = nothing)
    material = serialize_value(value.material; material_name)
    return Dict(
        "kind" => "region",
        "tag" => String(value.tag),
        "primitive" => serialize_value(value.primitive),
        "material" => material
    )
end
function serialize_value(value::Stack; material_name = nothing)
    return Dict(
        "kind" => "stack",
        "items" => [serialize_value(item; material_name) for item in value.items]
    )
end
function serialize_value(value::Group; material_name = nothing)
    return Dict(
        "kind" => "group",
        "name" => String(value.name),
        "at" => serialize_value(value.at),
        "item" => serialize_value(value.item; material_name),
        "pattern" => serialize_value(value.pattern),
        "path" => serialize_value(value.path),
        "compact" => serialize_value(value.compact),
        "boundary" => serialize_value(value.boundary)
    )
end
function serialize_value(
        value::Assembly{<:Pose2, <:AbstractCablePart};
        material_name = nothing
)
    names = value.names === nothing ? nothing : String.(value.names)
    return Dict(
        "kind" => "assembly",
        "at" => serialize_value(value.at),
        "item" => serialize_value(value.item; material_name),
        "pattern" => serialize_value(value.pattern),
        "path" => serialize_value(value.path),
        "compact" => serialize_value(value.compact),
        "names" => names
    )
end
function serialize_value(
        value::Assembly{<:Pose2, <:Tuple{Vararg{DataModel.AssemblyMember}}};
        material_name = nothing
)
    return Dict(
        "kind" => "assembly",
        "at" => serialize_value(value.at),
        "members" => [Dict(
             "at" => serialize_value(member.at),
             "item" => serialize_value(member.item; material_name)
         ) for member in value.item]
    )
end
function serialize_value(value::Enclosure; material_name = nothing)
    fill = serialize_value(value.fill; material_name)
    wall = value.wall === nothing ? nothing : serialize_value(value.wall; material_name)
    return Dict(
        "kind" => "enclosure",
        "tag" => String(value.tag),
        "at" => serialize_value(value.at),
        "primitive" => serialize_value(value.primitive),
        "item" => serialize_value(value.item; material_name),
        "fill" => fill,
        "wall" => wall
    )
end

function serialize_value(value::CableDesign; material_name = nothing)
    return Dict(
        "kind" => "cable_design",
        "cable_id" => value.cable_id,
        "nominal_data" => [Dict(
             "name" => String(name),
             "value" => serialize_value(item)
         ) for (name, item) in pairs(value.nominal_data)],
        "origin" => serialize_value(value.origin; material_name)
    )
end

function serialize_value(value::LineCableSystem)
    return Dict(
        "kind" => "line_cable_system",
        "system_id" => value.system_id,
        "line_length" => serialize_value(value.line_length),
        "designs" => [serialize_value(design) for design in value.designs],
        "positions" => serialize_value(value.positions),
        "connections" => serialize_value(value.connections),
        "environment" => serialize_value(value.environment)
    )
end

function serialize_value(value::Engine.LineParametersProblem)
    return Dict(
        "kind" => "line_parameters_problem",
        "system" => serialize_value(value.system),
        "temperature" => serialize_value(value.temperature),
        "earth_props" => serialize_value(value.earth_props),
        "frequencies" => serialize_value(value.frequencies)
    )
end

function serialize_value(grid::DeterministicGrid)
    return Dict("grid" => serialize_value(collect(grid.vals)))
end
function serialize_value(grid::RelativeGrid)
    return Dict(
        "grid" => serialize_value(collect(grid.vals)),
        "rel" => serialize_value(collect(grid.rel_err))
    )
end
function serialize_value(grid::AbsoluteGrid)
    return Dict(
        "grid" => serialize_value(collect(grid.vals)),
        "abs" => serialize_value(collect(grid.abs_err))
    )
end

function _collect_materials!(materials::Vector{AbstractMaterial}, part::Region)
    return _collect_materials!(materials, part.material)
end
function _collect_materials!(materials::Vector{AbstractMaterial}, material::AbstractMaterial)
    any(item -> isequal(item, material), materials) || push!(materials, material)
    return materials
end
function _collect_materials!(materials::Vector{AbstractMaterial}, part::Stack)
    foreach(item -> _collect_materials!(materials, item), part.items)
    return materials
end
function _collect_materials!(materials::Vector{AbstractMaterial}, part::Group)
    return _collect_materials!(materials, part.item)
end
function _collect_materials!(
        materials::Vector{AbstractMaterial},
        part::Assembly{<:Pose2, <:AbstractCablePart}
)
    return _collect_materials!(materials, part.item)
end
function _collect_materials!(
        materials::Vector{AbstractMaterial},
        part::Assembly{<:Pose2, <:Tuple{Vararg{DataModel.AssemblyMember}}}
)
    foreach(member -> _collect_materials!(materials, member.item), part.item)
    return materials
end
function _collect_materials!(materials::Vector{AbstractMaterial}, part::Enclosure)
    _collect_materials!(materials, part.item)
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
    materials = AbstractMaterial[]
    cable_ids = sort!(collect(keys(library.data)))
    for cable_id in cable_ids
        _collect_materials!(materials, library.data[cable_id].origin)
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
                    serialize_value(library.data[cable_id]; material_name),
                    Dict("datasheet" => Dict(
                        String(name) => serialize_value(item)
                    for (name, item) in pairs(library.datasheets[cable_id])
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

function _json_document(problem::Engine.LineParametersProblem)
    return Dict(
        "\$schema" => JSON_SCHEMA_DIALECT,
        "format" => PROBLEM_SCHEMA,
        "version" => JSON_SCHEMA_VERSION,
        "root" => serialize_value(problem)
    )
end
