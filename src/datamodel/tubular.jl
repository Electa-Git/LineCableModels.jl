"""
$(TYPEDEF)

Represent a tubular or solid conductor at material reference temperature.

$(TYPEDFIELDS)
"""
struct Tubular{T <: Real} <: AbstractConductorPart{T}
    r_in::T
    r_ex::T
    material_props::Material{T}
    cross_section::T
    resistance::T
    gmr::T
end

_preview_layer_name(::Tubular) = "tubular conductor"
preview_materials(layer::Tubular) = (layer.material_props,)
preview_shapes(layer::Tubular, context) = _annular_preview_shapes(layer, context)

function Validation.rules(::Type{<:Tubular})
    (Finite(:r_in), Finite(:r_ex), Nonnegative(:r_in), Positive(:r_ex),
        Less(:r_in, :r_ex))
end

function Tubular(r_in::Real, r_ex::Real, material::Material)
    T = promote_type(
        Float32, typeof(float(r_in)), typeof(float(r_ex)), eltype(material)
    )
    rin, rex = convert.(T, (r_in, r_ex))
    props = convert(Material{T}, material)
    Validation.check(Tubular, (; r_in = rin, r_ex = rex))
    area = (one(rin) * π) * (rex^2 - rin^2)
    resistance = tubular_resistance(
        rin, rex, props.rho, props.alpha, props.T0, props.T0
    )
    return validate(Tubular(
        rin, rex, props, area, resistance, tubular_gmr(rex, rin, props.mu_r)
    ))
end

function Base.convert(
        ::Type{AbstractConductorPart{T}},
        layer::Tubular
) where {T <: Real}
    return Tubular(
        convert(T, layer.r_in), convert(T, layer.r_ex),
        convert(Material{T}, layer.material_props)
    )
end
