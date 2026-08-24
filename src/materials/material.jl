"""
$(TYPEDEF)

Store electromagnetic and thermal material properties at a reference
temperature.

$(TYPEDFIELDS)
"""
struct Material{T <: Real}
    "Electrical resistivity of the material \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T
    "Reference temperature for property evaluations \\[°C\\]."
    T0::T
    "Temperature coefficient of resistivity \\[1/°C\\]."
    alpha::T

    @inline function Material{T}(
            rho::T,
            eps_r::T,
            mu_r::T,
            T0::T,
            alpha::T
    ) where {T <: Real}
        return validate(new{T}(rho, eps_r, mu_r, T0, alpha))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a material after floating and promoting its real-valued properties.

# Arguments
- `rho`: Resistivity \\[Ω·m\\].
- `eps_r`: Relative permittivity \\[dimensionless\\].
- `mu_r`: Relative permeability \\[dimensionless\\].
- `T0`: Reference temperature \\[°C\\].
- `alpha`: Temperature coefficient of resistivity \\[1/°C\\].

# Returns
- A validated `Material` whose scalar type is the promoted floating type of
  its inputs.

# Errors

- Throws `DomainError` when a property is outside its accepted physical
  domain or is not finite, except that infinite resistivity is accepted.
"""
@inline function Material(rho, eps_r, mu_r, T0, alpha)
    values = promote(float(rho), float(eps_r), float(mu_r), float(T0), float(alpha))
    return Material{typeof(first(values))}(values...)
end

function _check_material(material::Material)
    isnan(material.rho) && throw(DomainError(material.rho, "resistivity cannot be NaN"))
    material.rho > zero(material.rho) ||
        throw(DomainError(material.rho, "resistivity must be positive"))
    isfinite(material.eps_r) && material.eps_r >= zero(material.eps_r) ||
        throw(DomainError(material.eps_r, "relative permittivity must be nonnegative and finite"))
    isfinite(material.mu_r) && material.mu_r > zero(material.mu_r) ||
        throw(DomainError(material.mu_r, "relative permeability must be positive and finite"))
    isfinite(material.T0) ||
        throw(DomainError(material.T0, "reference temperature must be finite"))
    isfinite(material.alpha) ||
        throw(DomainError(material.alpha, "temperature coefficient must be finite"))
    return nothing
end

function Validation.rules(::Type{<:Material})
    (Validation.OwnerRule(:material_properties, _check_material),)
end
