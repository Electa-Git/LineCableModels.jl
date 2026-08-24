"""
$(TYPEDSIGNATURES)

Construct electromagnetic and thermal material properties. Scalar property
inputs return a [`Material`](@ref) directly. An explicit [`Grid`](@ref) or
nested [`Gridspace`](@ref) lifts the construction to a finite space.

# Keywords

- `rho`: Electrical resistivity \\[Ω·m\\].
- `eps_r=1.0`: Relative permittivity \\[dimensionless\\].
- `mu_r=1.0`: Relative permeability \\[dimensionless\\].
- `T0=20.0`: Reference temperature \\[°C\\].
- `alpha=0.0`: Temperature coefficient of resistivity \\[1/°C\\].
- `combine=:product`: Local composition rule when an input varies.

# Returns

- A [`Material`](@ref) for scalar inputs, or a `Gridspace{Material}` when at
  least one direct input is a `Grid` or nested `Gridspace`.
"""
function Material(;
        rho,
        eps_r = 1.0,
        mu_r = 1.0,
        T0 = 20.0,
        alpha = 0.0,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip; got :$combine"))
    values = (rho, eps_r, mu_r, T0, alpha)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{Materials.Material}(Materials.Material, grids; combine)
    end
    return Materials.Material(values...)
end

function Material(
        material::Materials.Material;
        rho = material.rho,
        eps_r = material.eps_r,
        mu_r = material.mu_r,
        T0 = material.T0,
        alpha = material.alpha,
        combine::Symbol = :product
)
    return Material(; rho, eps_r, mu_r, T0, alpha, combine)
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
