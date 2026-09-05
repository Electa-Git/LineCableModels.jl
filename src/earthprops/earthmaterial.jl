"""
$(TYPEDEF)

Store the electromagnetic properties of one earth material at one frequency.

`EarthMaterial` is an ephemeral constitutive value created while the Coaxial
backend evaluates an [`EarthLayer`](@ref). It is not stored in an
[`EarthModel`](@ref) or in a materials library.

$(TYPEDFIELDS)
"""
struct EarthMaterial{T <: Real} <: AbstractMaterial
    "Electrical resistivity \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]; may be negative for an artificial EHEM."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T

    @inline function EarthMaterial{T}(rho::T, eps_r::T, mu_r::T) where {T <: Real}
        return validate(new{T}(rho, eps_r, mu_r))
    end
end

Base.eltype(::EarthMaterial{T}) where {T} = T
Base.eltype(::Type{EarthMaterial{T}}) where {T} = T

function validate(material::EarthMaterial)
    isnan(material.rho) && throw(DomainError(
        material.rho,
        "EarthMaterial.rho must not be NaN"
    ))
    material.rho > zero(material.rho) || throw(DomainError(
        material.rho,
        "EarthMaterial.rho must be positive"
    ))
    isfinite(material.eps_r) && !iszero(material.eps_r) || throw(DomainError(
        material.eps_r,
        "EarthMaterial.eps_r must be nonzero and finite"
    ))
    isfinite(material.mu_r) && material.mu_r > zero(material.mu_r) ||
        throw(DomainError(
            material.mu_r,
            "EarthMaterial.mu_r must be positive and finite"
        ))
    return material
end

"""
$(TYPEDSIGNATURES)

Construct one ephemeral earth material after floating and promoting its
electromagnetic properties.

# Arguments

- `rho`: Electrical resistivity \\[Ω·m\\].
- `eps_r`: Relative permittivity \\[dimensionless\\]. Artificial equivalent
  earth models may produce a negative value.
- `mu_r`: Relative permeability \\[dimensionless\\].
"""
@inline function EarthMaterial(rho::Real, eps_r::Real, mu_r::Real)
    values = promote(float(rho), float(eps_r), float(mu_r))
    T = typeof(first(values))
    return EarthMaterial{T}(values...)
end

function Base.convert(::Type{EarthMaterial{T}}, material::EarthMaterial) where {T <: Real}
    return EarthMaterial{T}(
        convert(T, material.rho),
        convert(T, material.eps_r),
        convert(T, material.mu_r)
    )
end

function Base.convert(
        ::Type{EarthMaterial{T}}, material::EarthMaterial{T}
) where {T <: Real}
    material
end
