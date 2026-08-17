"""
$(TYPEDEF)

Represent one semiconducting layer at material reference temperature.

$(TYPEDFIELDS)
"""
struct Semicon{T <: Real} <: AbstractInsulatorPart{T}
    r_in::T
    r_ex::T
    material_props::Material{T}
    cross_section::T
    resistance::T
    gmr::T
    shunt_capacitance::T
    "Shunt conductance per unit length \\[S/m\\]."
    shunt_conductance::T
end

function Validation.rules(::Type{Semicon})
    (Finite(:r_in), Finite(:r_ex), Nonnegative(:r_in), Positive(:r_ex),
        Less(:r_in, :r_ex))
end
validate(layer::Semicon) = Validation.check(Semicon, layer)

function Semicon(r_in::Real, r_ex::Real, material::Material)
    T = promote_type(
        Float32, typeof(float(r_in)), typeof(float(r_ex)), eltype(material)
    )
    rin, rex = convert.(T, (r_in, r_ex))
    props = convert(Material{T}, material)
    Validation.check(Semicon, (; r_in = rin, r_ex = rex))
    return validate(Semicon(
        rin,
        rex,
        props,
        (one(rin) * π) * (rex^2 - rin^2),
        tubular_resistance(rin, rex, props.rho, props.alpha, props.T0, props.T0),
        tubular_gmr(rex, rin, props.mu_r),
        shunt_capacitance(rin, rex, props.eps_r),
        shunt_conductance(rin, rex, props.rho)
    ))
end

function Base.convert(
        ::Type{AbstractInsulatorPart{T}},
        layer::Semicon
) where {T <: Real}
    return Semicon(
        convert(T, layer.r_in), convert(T, layer.r_ex),
        convert(Material{T}, layer.material_props)
    )
end
