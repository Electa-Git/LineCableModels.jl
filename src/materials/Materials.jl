"""
    LineCableModels.Materials

Define the electromagnetic and thermal properties of the materials of a cable.
"""
module Materials

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..Commons: validate

export Material, material

"""
$(TYPEDEF)

Store the electromagnetic and thermal properties of one material. Its class, `:conductor`,
`:insulator` or `:semicon`, is the type parameter `K`, so laws and formulations dispatch on
it.

$(TYPEDFIELDS)
"""
struct Material{T <: Real, K}
    "Electrical resistivity at the reference temperature `T0` \\[Ω·m\\]."
    rho::T
    "Relative permittivity \\[dimensionless\\]."
    eps_r::T
    "Relative permeability \\[dimensionless\\]."
    mu_r::T
    "Reference temperature \\[°C\\]."
    T0::T
    "Temperature coefficient of resistivity \\[1/°C\\]."
    alpha::T
    "Thermal resistivity \\[K·m/W\\]."
    rho_thermal::T
    "Maximum continuous operating temperature \\[°C\\]."
    theta_max::T
    "Polarization loss tangent, apart from the conduction that `rho` gives \\[dimensionless\\]."
    tan_delta::T
    "Solar-absorption coefficient \\[dimensionless\\]."
    sigma_solar::T
    function Material{T, K}(rho::T, eps_r::T, mu_r::T, T0::T, alpha::T, rho_thermal::T,
            theta_max::T, tan_delta::T, sigma_solar::T) where {T <: Real, K}
        return validate(new{T, K}(
            rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max, tan_delta, sigma_solar))
    end
end

"""
$(TYPEDSIGNATURES)

Construct a material of class `K`, given as `Val(K)` so that the compiler knows the class. The arguments are the fields of [`Material`](@ref), in their order.

# Returns

- A `Material{T, K}`, where `T` is the floating type that the properties promote to.

# Errors

- Throws `DomainError` when a property is outside its physical domain or is not finite. The
  resistivity alone can be infinite.
- Throws `DomainError` when the class is not `:conductor`, `:insulator` or `:semicon`.
"""
function Material(::Val{K}, rho::Real, eps_r::Real, mu_r::Real, T0::Real, alpha::Real,
        rho_thermal::Real, theta_max::Real, tan_delta::Real, sigma_solar::Real) where {K}
    values = map(float,
        promote(
            rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max, tan_delta, sigma_solar))
    return Material{eltype(values), K}(values...)
end

"""
$(TYPEDSIGNATURES)

Construct a material of class `kind`. The arguments and keywords are the fields of
[`Material`](@ref).

# Returns

- A `Material{T, kind}`, where `T` is the floating type that the properties promote to.

# Errors

- Throws `DomainError` when a property is outside its physical domain or is not finite. The
  resistivity alone can be infinite.
- Throws `DomainError` when the class is not `:conductor`, `:insulator` or `:semicon`.
"""
function Material(kind::Symbol, rho::Real, eps_r::Real = 1, mu_r::Real = 1, T0::Real = 20,
        alpha::Real = 0; rho_thermal::Real = 0, theta_max::Real = 90, tan_delta::Real = 0,
        sigma_solar::Real = 0)
    return Material(Val(kind), rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max,
        tan_delta, sigma_solar)
end

"""
$(TYPEDSIGNATURES)

Build the [`Material`](@ref) of class `kind` from its properties, given as keywords with the
names and defaults of the positional constructor.

# Returns

- A `Material{T, kind}`, where `T` is the floating type that the properties promote to.
"""
function material(;
        kind::Symbol, rho::Real, eps_r::Real = 1, mu_r::Real = 1, T0::Real = 20,
        alpha::Real = 0, rho_thermal::Real = 0, theta_max::Real = 90, tan_delta::Real = 0,
        sigma_solar::Real = 0)
    return Material(Val(kind), rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max,
        tan_delta, sigma_solar)
end

function Base.map(f, value::Material{T, K}) where {T, K}
    return material(; kind = K, rho = f(value.rho), eps_r = f(value.eps_r),
        mu_r = f(value.mu_r), T0 = f(value.T0), alpha = f(value.alpha),
        rho_thermal = f(value.rho_thermal), theta_max = f(value.theta_max),
        tan_delta = f(value.tan_delta), sigma_solar = f(value.sigma_solar))
end

Base.eltype(::Material{T, K}) where {T, K} = T

function Base.convert(::Type{Material{T}}, material::Material{
        <:Real, K}) where {T <: Real, K}
    return Material{T, K}(convert(T, material.rho), convert(T, material.eps_r),
        convert(T, material.mu_r), convert(T, material.T0), convert(T, material.alpha),
        convert(T, material.rho_thermal), convert(T, material.theta_max),
        convert(T, material.tan_delta), convert(T, material.sigma_solar))
end

function validate(material::Material{<:Real, K}) where {K}
    throw(DomainError(material,
        "Material class :$K must be :conductor, :insulator or :semicon"))
end

function validate(material::Union{Material{T, :conductor}, Material{T, :insulator},
        Material{T, :semicon}}) where {T}
    isnan(material.rho) && throw(DomainError(material.rho, "Material.rho must not be NaN"))
    material.rho > zero(material.rho) ||
        throw(DomainError(material.rho, "Material.rho must be positive"))
    isfinite(material.eps_r) && material.eps_r >= zero(material.eps_r) ||
        throw(DomainError(material.eps_r, "Material.eps_r must be nonnegative and finite"))
    isfinite(material.mu_r) && material.mu_r > zero(material.mu_r) ||
        throw(DomainError(material.mu_r, "Material.mu_r must be positive and finite"))
    isfinite(material.T0) || throw(DomainError(material.T0, "Material.T0 must be finite"))
    isfinite(material.alpha) ||
        throw(DomainError(material.alpha, "Material.alpha must be finite"))
    isfinite(material.rho_thermal) && material.rho_thermal >= zero(material.rho_thermal) ||
        throw(DomainError(material.rho_thermal,
            "Material.rho_thermal must be nonnegative and finite"))
    isfinite(material.theta_max) ||
        throw(DomainError(material.theta_max, "Material.theta_max must be finite"))
    isfinite(material.tan_delta) && material.tan_delta >= zero(material.tan_delta) ||
        throw(DomainError(material.tan_delta,
            "Material.tan_delta must be nonnegative and finite"))
    isfinite(material.sigma_solar) && material.sigma_solar >= zero(material.sigma_solar) ||
        throw(DomainError(material.sigma_solar,
            "Material.sigma_solar must be nonnegative and finite"))
    return material
end

end
