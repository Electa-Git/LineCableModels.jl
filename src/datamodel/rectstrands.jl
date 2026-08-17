"Pure geometric state of one rectangular-strand layer."
struct RectStrandsShape{T <: Real, U <: Integer} <: AbstractShapeGeometry
    thickness::T
    width::T
    num_wires::U
    lay_ratio::T
    lay_direction::U
    mean_diameter::T
    pitch_length::T
    cross_section::T
end

"""
$(TYPEDEF)

Represent one area-preserving rectangular-strand layer at material reference
temperature.

$(TYPEDFIELDS)
"""
struct RectStrands{T <: Real, S <: RectStrandsShape} <: AbstractStrandsLayer{T}
    r_in::T
    r_ex::T
    material_props::Material{T}
    resistance::T
    gmr::T
    shape::S
end

function Validation.rules(::Type{RectStrands})
    return (
        Finite(:r_in), Finite(:r_ex), Finite(:thickness), Finite(:width),
        Finite(:lay_ratio),
        Nonnegative(:r_in), Positive(:r_ex), Less(:r_in, :r_ex),
        Positive(:thickness), Positive(:width), IntegerField(:num_wires),
        Positive(:num_wires), Nonnegative(:lay_ratio),
        OneOf(:lay_direction, (-1, 1)),
        PhysicalFillLimit(:num_wires, (:r_in, :width))
    )
end

function validate(layer::RectStrands)
    Validation.check(RectStrands,
        (
            r_in = layer.r_in,
            r_ex = layer.r_ex,
            thickness = layer.shape.thickness,
            width = layer.shape.width,
            num_wires = layer.shape.num_wires,
            lay_ratio = layer.shape.lay_ratio,
            lay_direction = layer.shape.lay_direction
        ))
    expected = sqrt(
        layer.r_in^2 +
        layer.shape.num_wires * layer.shape.width *
        layer.shape.thickness / (one(layer.r_in) * π),
    )
    isapprox(
        layer.r_ex, expected;
        rtol = sqrt(eps(typeof(float(nominal(layer.r_ex))))),
        atol = sqrt(eps(Float64))
    ) || throw(DomainError(
        layer.r_ex,
        "RectStrands outer radius is inconsistent with its area-preserving geometry"
    ))
    return layer
end

function RectStrands(
        r_in::Real,
        r_ex::Real,
        thickness::Real,
        width::Real,
        num_wires::Integer,
        lay_ratio::Real,
        material::Material;
        lay_direction::Integer = 1
)
    T = promote_type(
        Float32,
        typeof(float(r_in)), typeof(float(r_ex)), typeof(float(thickness)),
        typeof(float(width)), typeof(float(lay_ratio)), eltype(material)
    )
    rin, rex, thick, w, ratio = convert.(
        T, (r_in, r_ex, thickness, width, lay_ratio)
    )
    props = convert(Material{T}, material)
    candidate = (; r_in = rin, r_ex = rex, thickness = thick, width = w, num_wires,
        lay_ratio = ratio, lay_direction)
    Validation.check(RectStrands, candidate)
    expected = sqrt(rin^2 + num_wires * w * thick / (one(rin) * π))
    isapprox(rex, expected; rtol = sqrt(eps(T)), atol = sqrt(eps(Float64))) ||
        throw(DomainError(
            rex,
            "RectStrands outer radius must equal the area-preserving radius $expected"
        ))
    mean_diameter, pitch_length, overlength = helix(rin, rex, ratio)
    cross_section = num_wires * w * thick
    shape = RectStrandsShape(
        thick, w, num_wires, ratio, lay_direction, mean_diameter,
        pitch_length, cross_section
    )
    resistance = strip_resistance(
        thick, w, props.rho, props.alpha, props.T0, props.T0
    ) * overlength / num_wires
    return validate(RectStrands(
        rin, rex, props, resistance, tubular_gmr(rex, rin, props.mu_r), shape
    ))
end

function maxfill(::Type{RectStrands}, lay_radius::Real, width::Real)
    lay_radius >= zero(lay_radius) ||
        throw(DomainError(lay_radius, "lay radius must be nonnegative"))
    width > zero(width) || throw(DomainError(width, "width must be positive"))
    count = 2π * lay_radius / width
    value = float(nominal(count))
    return floor(Int, value + 8eps(value))
end

function Base.convert(
        ::Type{AbstractConductorPart{T}},
        layer::RectStrands
) where {T <: Real}
    rin = convert(T, layer.r_in)
    thickness = convert(T, layer.shape.thickness)
    width = convert(T, layer.shape.width)
    return RectStrands(
        rin, convert(T, layer.r_ex), thickness, width,
        layer.shape.num_wires, convert(T, layer.shape.lay_ratio),
        convert(Material{T}, layer.material_props);
        lay_direction = layer.shape.lay_direction
    )
end
