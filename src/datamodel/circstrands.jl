"""
$(TYPEDEF)

Represent one circular-wire layer at the material reference temperature.

$(TYPEDFIELDS)
"""
struct CircStrands{T <: Real, U <: Integer} <: AbstractStrandsLayer{T}
    "Inner radial boundary \\[m\\]."
    r_in::T
    "Outer radial boundary \\[m\\]."
    r_ex::T
    "Wire radius \\[m\\]."
    radius_wire::T
    "Number of wires \\[dimensionless\\]."
    num_wires::U
    "Lay-length ratio \\[dimensionless\\]."
    lay_ratio::T
    "Mean helix diameter \\[m\\]."
    mean_diameter::T
    "Helix pitch length \\[m\\]."
    pitch_length::T
    "Lay direction: `1` or `-1` \\[dimensionless\\]."
    lay_direction::U
    "Wire material at its reference state."
    material_props::Material{T}
    "Total metallic cross-section \\[m²\\]."
    cross_section::T
    "Layer resistance \\[Ω/m\\]."
    resistance::T
    "Geometric mean radius \\[m\\]."
    gmr::T
end

function BaseParams.gmd_elements(layer::CircStrands)
    coordinates = wire_coordinates(
        layer.num_wires, layer.radius_wire, layer.r_in
    )
    return (coordinates, layer.radius_wire, π * layer.radius_wire^2)
end

function Validation.rules(::Type{CircStrands})
    return (
        Finite(:r_in), Finite(:r_ex), Finite(:radius_wire), Finite(:lay_ratio),
        Nonnegative(:r_in), Positive(:r_ex), Less(:r_in, :r_ex),
        Positive(:radius_wire), IntegerField(:num_wires), Positive(:num_wires),
        Nonnegative(:lay_ratio), OneOf(:lay_direction, (-1, 1)),
        PhysicalFillLimit(:num_wires, (:r_in, :radius_wire))
    )
end

function validate(layer::CircStrands)
    Validation.check(CircStrands, layer)
    isapprox(layer.r_ex,
        layer.num_wires == 1 ? layer.radius_wire :
        layer.r_in + 2 * layer.radius_wire;
        rtol = sqrt(eps(typeof(float(nominal(layer.r_ex))))),
        atol = sqrt(eps(Float64))) || throw(DomainError(
        layer.r_ex,
        "CircStrands outer radius is inconsistent with its wire geometry"
    ))
    return layer
end

"""
$(TYPEDSIGNATURES)

Construct a circular-wire layer from resolved radial geometry. Electrical
properties are evaluated at `material.T0`.
"""
function CircStrands(
        r_in::Real,
        r_ex::Real,
        radius_wire::Real,
        num_wires::Integer,
        lay_ratio::Real,
        material::Material;
        lay_direction::Integer = 1
)
    T = promote_type(
        Float32,
        typeof(float(r_in)), typeof(float(r_ex)), typeof(float(radius_wire)),
        typeof(float(lay_ratio)), eltype(material)
    )
    rin, rex, rw, ratio = convert.(T, (r_in, r_ex, radius_wire, lay_ratio))
    props = convert(Material{T}, material)
    candidate = (; r_in = rin, r_ex = rex, radius_wire = rw, num_wires,
        lay_ratio = ratio, lay_direction)
    Validation.check(CircStrands, candidate)
    expected = num_wires == 1 ? rw : rin + 2 * rw
    isapprox(rex, expected; rtol = sqrt(eps(T)), atol = sqrt(eps(Float64))) ||
        throw(DomainError(
            rex,
            "CircStrands outer radius must equal $expected for the supplied wires"
        ))
    mean_diameter, pitch_length, overlength = helix(rin, rex, ratio)
    cross_section = num_wires * (one(rin) * π) * rw^2
    resistance = tubular_resistance(
        zero(T), rw, props.rho, props.alpha, props.T0, props.T0
    ) * overlength / num_wires
    layer = CircStrands(
        rin, rex, rw, num_wires, ratio, mean_diameter, pitch_length,
        lay_direction, props, cross_section, resistance,
        strand_gmr(rin + rw, num_wires, rw, props.mu_r)
    )
    return validate(layer)
end

function maxfill(::Type{CircStrands}, r_in::Real, wire_radius::Real)
    r_in >= zero(r_in) || throw(DomainError(r_in, "inner radius must be nonnegative"))
    wire_radius > zero(wire_radius) ||
        throw(DomainError(wire_radius, "wire radius must be positive"))
    iszero(r_in) && return 1
    count = π / asin(wire_radius / (r_in + wire_radius))
    value = float(nominal(count))
    return floor(Int, value + 8eps(value))
end

function Base.convert(
        ::Type{AbstractConductorPart{T}},
        layer::CircStrands
) where {T <: Real}
    rin = convert(T, layer.r_in)
    radius = convert(T, layer.radius_wire)
    return CircStrands(
        rin, convert(T, layer.r_ex), radius, layer.num_wires,
        convert(T, layer.lay_ratio), convert(Material{T}, layer.material_props);
        lay_direction = layer.lay_direction
    )
end
