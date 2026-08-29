"""
$(TYPEDEF)

Store the static physical properties of one earth or air layer.

$(TYPEDFIELDS)
"""
struct EarthLayer{T <: Real}
    "Electrical resistivity \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T
    "Layer thickness \\[m\\]."
    thickness::T

    function EarthLayer{T}(rho::T, eps_r::T, mu_r::T, thickness::T) where {T <: Real}
        return validate(new{T}(rho, eps_r, mu_r, thickness))
    end
end

Base.eltype(::EarthLayer{T}) where {T} = T
Base.eltype(::Type{EarthLayer{T}}) where {T} = T

function _check_earth_layer(layer::EarthLayer)
    isnan(layer.rho) && throw(DomainError(layer.rho, "resistivity cannot be NaN"))
    layer.rho > zero(layer.rho) ||
        throw(DomainError(layer.rho, "resistivity must be positive"))
    isfinite(layer.eps_r) && layer.eps_r > zero(layer.eps_r) ||
        throw(DomainError(layer.eps_r, "relative permittivity must be positive and finite"))
    isfinite(layer.mu_r) && layer.mu_r > zero(layer.mu_r) ||
        throw(DomainError(layer.mu_r, "relative permeability must be positive and finite"))
    layer.thickness > zero(layer.thickness) ||
        throw(DomainError(layer.thickness, "thickness must be positive"))
    return nothing
end

function Validation.rules(::Type{<:EarthLayer})
    (Validation.OwnerRule(:earth_layer_properties, _check_earth_layer),)
end

"""
$(TYPEDSIGNATURES)

Construct one static earth layer after floating and promoting its properties.
"""
function EarthLayer(rho::Real, eps_r::Real, mu_r::Real, thickness::Real)
    values = promote(float(rho), float(eps_r), float(mu_r), float(thickness))
    T = typeof(first(values))
    return EarthLayer{T}(values...)
end

function EarthLayer(rho::Real, eps_r::Real, mu_r::Real)
    values = promote(float(rho), float(eps_r), float(mu_r))
    T = typeof(first(values))
    return EarthLayer{T}(values..., T(Inf))
end

function Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer) where {T <: Real}
    return EarthLayer{T}(
        convert(T, layer.rho), convert(T, layer.eps_r),
        convert(T, layer.mu_r), convert(T, layer.thickness)
    )
end

Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer{T}) where {T <: Real} = layer

"Construct the ephemeral electromagnetic material represented by an earth layer."
EarthMaterial(layer::EarthLayer{T}) where {T <: Real} =
    EarthMaterial{T}(layer.rho, layer.eps_r, layer.mu_r)
