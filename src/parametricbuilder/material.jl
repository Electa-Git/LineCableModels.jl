"""
$(TYPEDSIGNATURES)

Construct electromagnetic and thermal material properties. Scalar property
inputs return a [`Material`](@ref) directly. An explicit [`Grid`](@ref) or
nested [`Gridspace`](@ref) lifts the construction to a finite space.

# Keywords

- `rho`: Electrical resistivity \\[Ω·m\\].
- `kind`: Broad physical class. Deterministic symbol grids are accepted.
- `eps_r=1`: Relative permittivity \\[dimensionless\\].
- `mu_r=1`: Relative permeability \\[dimensionless\\].
- `T0=20`: Reference temperature \\[°C\\].
- `alpha=0`: Temperature coefficient of resistivity \\[1/°C\\].
- `rho_thermal=0`: Thermal resistivity \\[K·m/W\\].
- `theta_max=90`: Maximum continuous operating temperature \\[°C\\].
- `tan_delta=0`: Dielectric loss tangent.
- `sigma_solar=0`: Solar-absorption coefficient.
- `combine=:product`: Local composition rule when an input varies.

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
    values = (
        kind, rho, eps_r, mu_r, T0, alpha, rho_thermal,
        theta_max, tan_delta, sigma_solar
    )
    return _construction(
        Materials.Material, Materials.Material, values; combine
    )
end

function Material(
        material::Materials.Material;
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

function Material(
        library::Materials.MaterialsLibrary,
        name::Union{AbstractString, Symbol};
        kwargs...
)
    material = get(library, String(name), nothing)
    material === nothing && throw(KeyError(String(name)))
    return Material(material; kwargs...)
end

"""
$(TYPEDSIGNATURES)

Add the single deterministic material represented by `space` to a material
library.

# Errors

- Throws `ArgumentError` when `space` contains uncertainty or describes more
  than one material.
"""
function add!(
        library::Materials.MaterialsLibrary,
        name::Union{AbstractString, Symbol},
        space::Gridspace{Materials.Material}
)
    has_uncertainty(space) && throw(ArgumentError(
        "a reusable material-library entry must be deterministic",
    ))
    length(space) == 1 || throw(ArgumentError(
        "a reusable material-library entry must describe exactly one material; got $(length(space)) points",
    ))
    return add!(library, String(name), only(space))
end
