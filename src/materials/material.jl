"""
Supertype for package-owned electromagnetic material values.
"""
abstract type AbstractMaterial end

"""
$(TYPEDEF)

Store one broad material class and its electromagnetic and thermal properties. The class
`K`, such as `:conductor`, `:insulator` or `:semicon`, is a type parameter, so laws and
formulations dispatch on it. `material.kind` returns it.

$(TYPEDFIELDS)
"""
struct Material{T <: Real, K} <: AbstractMaterial
    "Electrical resistivity at reference temperature T0 \\[Ω·m\\]."
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
    "Polarization loss tangent, excluding conduction represented by rho \\[dimensionless\\]."
    tan_delta::T
    "Solar-absorption coefficient \\[dimensionless\\]."
    sigma_solar::T

    @inline function Material{T, K}(
            rho::T,
            eps_r::T,
            mu_r::T,
            T0::T,
            alpha::T,
            rho_thermal::T,
            theta_max::T,
            tan_delta::T,
            sigma_solar::T
    ) where {T <: Real, K}
        return validate(new{T, K}(
            rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max, tan_delta, sigma_solar
        ))
    end
end

@inline function Material{T}(kind::Symbol, rho::T, eps_r::T, mu_r::T, T0::T, alpha::T,
        rho_thermal::T, theta_max::T, tan_delta::T, sigma_solar::T) where {T <: Real}
    return Material{T, kind}(
        rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max, tan_delta, sigma_solar)
end

# The class is the type parameter: `material.kind` reads it.
@inline function Base.getproperty(material::Material{T, K}, name::Symbol) where {T, K}
    return name === :kind ? K : getfield(material, name)
end
Base.propertynames(::Material) = (:kind, fieldnames(Material)...)

"""
$(TYPEDSIGNATURES)

Construct a material after floating and promoting its real-valued properties.

# Arguments
- `kind`: broad physical class, such as `:conductor`, `:insulator`, or
  `:semicon`.
- `rho`: resistivity \\[Ω·m\\].
- `eps_r`: relative permittivity \\[dimensionless\\].
- `mu_r`: relative permeability \\[dimensionless\\].
- `T0`: reference temperature \\[°C\\].
- `alpha`: temperature coefficient of resistivity \\[1/°C\\].

# Keywords

- `rho_thermal=0`: thermal resistivity \\[K·m/W\\].
- `theta_max=90`: maximum continuous operating temperature \\[°C\\].
- `tan_delta=0`: polarization loss tangent \\[dimensionless\\], excluding the
  conduction contribution already represented by `rho`. A measured total loss
  tangent must be separated into those contributions before supplying both.
  Lossless constitutive selections explicitly suppress this input and conductivity.
- `sigma_solar=0`: solar-absorption coefficient.

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
    return Material(Val(kind), rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max,
        tan_delta, sigma_solar)
end

# Given as `Val(kind)`, the class becomes the type parameter at compile time. A builder
# with a fixed class then infers the concrete material type.
@inline function Material(::Val{K}, rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max,
        tan_delta, sigma_solar) where {K}
    values = map(float, promote(
        rho, eps_r, mu_r, T0, alpha, rho_thermal, theta_max,
        tan_delta, sigma_solar
    ))
    return Material{typeof(first(values)), K}(values...)
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

"""
$(TYPEDSIGNATURES)

Construct electromagnetic and thermal material properties. Scalar property
inputs return a [`Material`](@ref) directly. An explicit `Grid` or nested
`Gridspace` lifts the construction to a finite space.

# Keywords

- `rho`: electrical resistivity \\[Ω·m\\].
- `kind`: broad physical class. Deterministic symbol grids are accepted.
- `eps_r=1`: relative permittivity \\[dimensionless\\].
- `mu_r=1`: relative permeability \\[dimensionless\\].
- `T0=20`: reference temperature \\[°C\\].
- `alpha=0`: temperature coefficient of resistivity \\[1/°C\\].
- `rho_thermal=0`: thermal resistivity \\[K·m/W\\].
- `theta_max=90`: maximum continuous operating temperature \\[°C\\].
- `tan_delta=0`: dielectric loss tangent.
- `sigma_solar=0`: solar-absorption coefficient.
- `combine=:product`: local composition rule when an input varies.

# Returns

- A [`Material`](@ref) for scalar inputs, or a `Gridspace{Material}` when at
  least one direct input is a `Grid` or nested `Gridspace`.
"""
function Material(;
        kind,
        rho,
        eps_r = 1,
        mu_r = 1,
        T0 = 20,
        alpha = 0,
        rho_thermal = 0,
        theta_max = 90,
        tan_delta = 0,
        sigma_solar = 0,
        combine::Symbol = :product
)
    # A fixed class enters as `Val(kind)`, so a Gridspace over the other properties infers
    # its concrete material type.
    values = (
        kind isa Symbol ? Val(kind) : kind, rho, eps_r, mu_r, T0, alpha, rho_thermal,
        theta_max, tan_delta, sigma_solar
    )
    return parameterize(Material, Material, values; combine)
end

function Material(
        material::Material;
        kind = material.kind,
        rho = material.rho,
        eps_r = material.eps_r,
        mu_r = material.mu_r,
        T0 = material.T0,
        alpha = material.alpha,
        rho_thermal = material.rho_thermal,
        theta_max = material.theta_max,
        tan_delta = material.tan_delta,
        sigma_solar = material.sigma_solar,
        combine::Symbol = :product
)
    return Material(;
        kind, rho, eps_r, mu_r, T0, alpha, rho_thermal,
        theta_max, tan_delta, sigma_solar, combine
    )
end

function Material(unexpected; kwargs...)
    throw(ArgumentError(
        "keyword Material construction accepts no positional arguments; " *
        "got $(repr(unexpected)). To declare uncertainty, keep the error " *
        "inside Grid: `property = Grid(values, relative_error)` or " *
        "`property = Grid(values, AbsoluteError(error))`.",
    ))
end

function validate(material::Material)
    validate(material.kind, Material)
    isnan(material.rho) && throw(DomainError(
        material.rho,
        "Material.rho must not be NaN"
    ))
    material.rho > zero(material.rho) ||
        throw(DomainError(material.rho, "Material.rho must be positive"))
    isfinite(material.eps_r) && material.eps_r >= zero(material.eps_r) ||
        throw(DomainError(
            material.eps_r,
            "Material.eps_r must be nonnegative and finite"
        ))
    isfinite(material.mu_r) && material.mu_r > zero(material.mu_r) ||
        throw(DomainError(
            material.mu_r,
            "Material.mu_r must be positive and finite"
        ))
    isfinite(material.T0) ||
        throw(DomainError(material.T0, "Material.T0 must be finite"))
    isfinite(material.alpha) ||
        throw(DomainError(material.alpha, "Material.alpha must be finite"))
    isfinite(material.rho_thermal) && material.rho_thermal >= zero(material.rho_thermal) ||
        throw(DomainError(
            material.rho_thermal,
            "Material.rho_thermal must be nonnegative and finite"
        ))
    isfinite(material.theta_max) ||
        throw(DomainError(material.theta_max, "Material.theta_max must be finite"))
    isfinite(material.tan_delta) && material.tan_delta >= zero(material.tan_delta) ||
        throw(DomainError(
            material.tan_delta,
            "Material.tan_delta must be nonnegative and finite"
        ))
    isfinite(material.sigma_solar) && material.sigma_solar >= zero(material.sigma_solar) ||
        throw(DomainError(
            material.sigma_solar,
            "Material.sigma_solar must be nonnegative and finite"
        ))
    return material
end

# The kind of a material is a nonempty Symbol.
function validate(kind::Symbol, ::Type{Material})
    kind === Symbol("") && throw(ArgumentError(
        "Material.kind must be a nonempty Symbol; received $(repr(kind))"
    ))
    return kind
end
function validate(kind, ::Type{Material})
    throw(ArgumentError("Material.kind must be a nonempty Symbol; received $(repr(kind))"))
end

Commons.input_fields(::Type{<:Material}) = (rho=(name="electrical resistivity",unit="Ω·m"),
    eps_r=(name="relative permittivity",unit=""),mu_r=(name="relative permeability",unit=""),
    T0=(name="reference temperature",unit="°C"),alpha=(name="temperature coefficient",unit="1/°C"),
    rho_thermal=(name="thermal resistivity",unit="K·m/W"),theta_max=(name="maximum temperature",unit="°C"),
    tan_delta=(name="loss tangent",unit=""),sigma_solar=(name="solar absorption coefficient",unit=""))
