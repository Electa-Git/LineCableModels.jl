"""
$(TYPEDEF)

Represent a helical conductive strip at material reference temperature.

$(TYPEDFIELDS)
"""
struct Strip{T <: Real} <: AbstractConductorPart{T}
    r_in::T
    r_ex::T
    thickness::T
    width::T
    lay_ratio::T
    mean_diameter::T
    pitch_length::T
    lay_direction::Int
    material_props::Material{T}
    cross_section::T
    resistance::T
    gmr::T
end

function Validation.rules(::Type{Strip})
    (
        Finite(:r_in), Finite(:r_ex), Finite(:width), Finite(:lay_ratio),
        Nonnegative(:r_in), Positive(:r_ex), Less(:r_in, :r_ex), Positive(:width),
        Nonnegative(:lay_ratio), OneOf(:lay_direction, (-1, 1))
    )
end
validate(layer::Strip) = Validation.check(Strip, layer)

function Strip(
        r_in::Real,
        r_ex::Real,
        width::Real,
        lay_ratio::Real,
        material::Material;
        lay_direction::Integer = 1
)
    T = promote_type(
        Float32,
        typeof(float(r_in)), typeof(float(r_ex)), typeof(float(width)),
        typeof(float(lay_ratio)), eltype(material)
    )
    rin, rex, w, ratio = convert.(T, (r_in, r_ex, width, lay_ratio))
    props = convert(Material{T}, material)
    Validation.check(Strip, (; r_in = rin, r_ex = rex, width = w, lay_ratio = ratio,
        lay_direction))
    thickness = rex - rin
    mean_diameter, pitch_length, overlength = helix(rin, rex, ratio)
    resistance = strip_resistance(
        thickness, w, props.rho, props.alpha, props.T0, props.T0
    ) * overlength
    return validate(Strip(
        rin, rex, thickness, w, ratio, mean_diameter, pitch_length,
        Int(lay_direction), props, thickness * w, resistance,
        tubular_gmr(rex, rin, props.mu_r)
    ))
end

function Base.convert(
        ::Type{AbstractConductorPart{T}},
        layer::Strip
) where {T <: Real}
    return Strip(
        convert(T, layer.r_in), convert(T, layer.r_ex), convert(T, layer.width),
        convert(T, layer.lay_ratio), convert(Material{T}, layer.material_props);
        lay_direction = layer.lay_direction
    )
end
