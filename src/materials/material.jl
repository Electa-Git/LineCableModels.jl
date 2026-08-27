"""
$(TYPEDEF)

Store one broad material class and its electromagnetic and thermal properties.

$(TYPEDFIELDS)
"""
struct Material{T <: Real}
    "Broad physical class used by formulation dispatch."
    kind::Symbol
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
    "Thermal resistivity \\[K·m/W\\]."
    rho_thermal::T
    "Maximum continuous operating temperature \\[°C\\]."
    theta_max::T
    "Dielectric loss tangent \\[dimensionless\\]."
    tan_delta::T
    "Solar-absorption coefficient \\[dimensionless\\]."
    sigma_solar::T

    @inline function Material{T}(
            kind::Symbol,
            rho::T,
            eps_r::T,
            mu_r::T,
            T0::T,
            alpha::T,
            rho_thermal::T,
            theta_max::T,
            tan_delta::T,
            sigma_solar::T
    ) where {T <: Real}
        return validate(new{T}(
            kind, rho, eps_r, mu_r, T0, alpha, rho_thermal,
            theta_max, tan_delta, sigma_solar
        ))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a material after floating and promoting its real-valued properties.

# Arguments
- `kind`: Broad physical class, such as `:conductor`, `:insulator`, or
  `:semicon`.
- `rho`: Resistivity \\[Ω·m\\].
- `eps_r`: Relative permittivity \\[dimensionless\\].
- `mu_r`: Relative permeability \\[dimensionless\\].
- `T0`: Reference temperature \\[°C\\].
- `alpha`: Temperature coefficient of resistivity \\[1/°C\\].

# Keywords

- `rho_thermal=0`: Thermal resistivity \\[K·m/W\\].
- `theta_max=90`: Maximum continuous operating temperature \\[°C\\].
- `tan_delta=0`: Dielectric loss tangent.
- `sigma_solar=0`: Solar-absorption coefficient.

# Returns
- A validated `Material` whose scalar type is the promoted floating type of
  its inputs.

# Errors

- Throws `DomainError` when a property is outside its accepted physical
  domain or is not finite, except that infinite resistivity is accepted.
"""
@inline function Material(
        kind::Symbol,
        rho,
        eps_r,
        mu_r,
        T0,
        alpha,
        rho_thermal,
        theta_max,
        tan_delta,
        sigma_solar
)
    values = map(float, promote(
        rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max,
        tan_delta, sigma_solar
    ))
    return Material{typeof(first(values))}(kind, values...)
end

@inline function Material(
        kind::Symbol,
        rho,
        eps_r = 1,
        mu_r = 1,
        T0 = 20,
        alpha = 0;
        rho_thermal = 0,
        theta_max = 90,
        tan_delta = 0,
        sigma_solar = 0
)
    return Material(
        kind, rho, eps_r, mu_r, T0, alpha, rho_thermal,
        theta_max, tan_delta, sigma_solar
    )
end

function _check_material(material::Material)
    isempty(string(material.kind)) && throw(ArgumentError("material kind cannot be empty"))
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
    isfinite(material.rho_thermal) && material.rho_thermal >= zero(material.rho_thermal) ||
        throw(DomainError(
            material.rho_thermal,
            "thermal resistivity must be nonnegative and finite"
        ))
    isfinite(material.theta_max) ||
        throw(DomainError(material.theta_max, "maximum temperature must be finite"))
    isfinite(material.tan_delta) && material.tan_delta >= zero(material.tan_delta) ||
        throw(DomainError(
            material.tan_delta,
            "dielectric loss tangent must be nonnegative and finite"
        ))
    isfinite(material.sigma_solar) && material.sigma_solar >= zero(material.sigma_solar) ||
        throw(DomainError(
            material.sigma_solar,
            "solar-absorption coefficient must be nonnegative and finite"
        ))
    return nothing
end

function Validation.rules(::Type{<:Material})
    (Validation.OwnerRule(:material_properties, _check_material),)
end
