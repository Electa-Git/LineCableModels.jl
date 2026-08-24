struct MaterialDefinition{R, E, M, T, A, C} <: _AbstractDefinition{Materials.Material}
    rho::R
    eps_r::E
    mu_r::M
    T0::T
    alpha::A
    combine::C
end

function Gridspace(spec::MaterialDefinition)
    return Gridspace{Materials.Material}(
        Materials.Material,
        map(_gridspace_axis, (
            spec.rho, spec.eps_r, spec.mu_r, spec.T0, spec.alpha
        )),
        (:rho, :eps_r, :mu_r, :T0, :alpha);
        combine = _valof(spec.combine)
    )
end

"""
$(TYPEDSIGNATURES)

Construct electromagnetic and thermal material properties. Scalar property
inputs return a [`Material`](@ref) directly. An explicit [`Grid`](@ref) or
nested declarative input lifts the same declaration to a
[`Gridspace{Material}`](@ref).

# Keywords

- `rho`: Electrical resistivity \\[Ω·m\\].
- `eps_r=1.0`: Relative permittivity \\[dimensionless\\].
- `mu_r=1.0`: Relative permeability \\[dimensionless\\].
- `T0=20.0`: Reference temperature \\[°C\\].
- `alpha=0.0`: Temperature coefficient of resistivity \\[1/°C\\].
- `combine=:product`: Local composition rule when an input varies.

# Returns

- A [`Material`](@ref) for scalar inputs, or a [`Gridspace{Material}`](@ref)
  when at least one input is a `Grid` or `_AbstractDefinition`.
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
    any(value -> value isa Union{AbstractGrid, Gridspace, _AbstractDefinition}, values) &&
        return Gridspace(MaterialDefinition(rho, eps_r, mu_r, T0, alpha, Val(combine)))
    return Materials.Material(rho, eps_r, mu_r, T0, alpha)
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

Add the single deterministic material described by `spec` to a
[`LineCableModels.MaterialsLibrary`](@ref). The definition is materialized
through the same [`Gridspace`](@ref) grammar used by the cable builder.

# Arguments

- `library`: Material library to modify.
- `name`: Material key.
- `spec`: Deterministic parameterized material definition.

# Returns

- The modified material library.

# Errors

- Throws `ArgumentError` when `spec` contains uncertainty or describes more
  than one material configuration.
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
        "a reusable material-library entry must describe exactly one material; got $(length(space)) configurations",
    ))
    return add!(library, String(name), only(space))
end
