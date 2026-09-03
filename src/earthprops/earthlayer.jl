"""
$(TYPEDEF)

Store the static physical properties of one earth or air layer.

$(TYPEDFIELDS)
"""
struct EarthLayer{T <: Real} <: AbstractEarthModel
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

function validate(layer::EarthLayer)
    isnan(layer.rho) && throw(DomainError(
        layer.rho,
        "EarthLayer.rho must not be NaN"
    ))
    layer.rho > zero(layer.rho) ||
        throw(DomainError(layer.rho, "EarthLayer.rho must be positive"))
    isfinite(layer.eps_r) && layer.eps_r > zero(layer.eps_r) ||
        throw(DomainError(
            layer.eps_r,
            "EarthLayer.eps_r must be positive and finite"
        ))
    isfinite(layer.mu_r) && layer.mu_r > zero(layer.mu_r) ||
        throw(DomainError(
            layer.mu_r,
            "EarthLayer.mu_r must be positive and finite"
        ))
    layer.thickness > zero(layer.thickness) ||
        throw(DomainError(layer.thickness, "EarthLayer.thickness must be positive"))
    return layer
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

function _earth_layer(rho, eps_r, mu_r, thickness)
    permittivity = eps_r === nothing ? one(float(rho)) : eps_r
    permeability = mu_r === nothing ? one(float(rho)) : mu_r
    return thickness === nothing ?
           EarthLayer(rho, permittivity, permeability) :
           EarthLayer(rho, permittivity, permeability, thickness)
end

"""
$(TYPEDSIGNATURES)

Declare one static earth layer directly or as an explicit finite space.

# Keywords

- `rho`: Electrical resistivity \\[Ω·m\\].
- `eps_r=nothing`: Relative permittivity \\[dimensionless\\]; `nothing`
  selects unity in the resistivity scalar type.
- `mu_r=nothing`: Relative permeability \\[dimensionless\\]; `nothing`
  selects unity in the resistivity scalar type.
- `thickness=nothing`: Layer thickness \\[m\\]; `nothing` selects a
  semi-infinite layer.
- `combine=:product`: Gridspace composition rule.

# Returns

- An `EarthLayer`, or a `Gridspace{EarthLayer}` when an explicit finite source
  is supplied.
"""
function layer(;
        rho,
        eps_r = nothing,
        mu_r = nothing,
        thickness = nothing,
        combine::Symbol = :product
)
    values = (rho, eps_r, mu_r, thickness)
    return parameterize(EarthLayer, _earth_layer, values; combine)
end

function Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer) where {T <: Real}
    return EarthLayer{T}(
        convert(T, layer.rho), convert(T, layer.eps_r),
        convert(T, layer.mu_r), convert(T, layer.thickness)
    )
end

Base.convert(::Type{EarthLayer{T}}, layer::EarthLayer{T}) where {T <: Real} = layer

"Construct the ephemeral electromagnetic material represented by an earth layer."
function EarthMaterial(layer::EarthLayer{T}) where {T <: Real}
    EarthMaterial{T}(layer.rho, layer.eps_r, layer.mu_r)
end
