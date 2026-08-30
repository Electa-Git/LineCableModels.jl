"Return the reconstruction assumptions of the Xue et al. EHEM."
assumptions(::Val{:Xue2021}) = (layer = -1,)

"""
$(TYPEDSIGNATURES)

**Identification.** Equivalent transverse propagation constant, from which
both conductivity and relative permittivity are reconstructed.

**Expression.** With
``\\gamma_{e,k}^2=\\gamma_k^2-\\gamma_0^2``, the upward recursion is

```math
\\gamma_{e,eq,k}=\\gamma_{e,k}
\\frac{\\gamma_{e,k}+\\gamma_{e,eq,k+1}-
(\\gamma_{e,k}-\\gamma_{e,eq,k+1})e^{-2h_k\\gamma_{e,k}}}
{\\gamma_{e,k}+\\gamma_{e,eq,k+1}+
(\\gamma_{e,k}-\\gamma_{e,eq,k+1})e^{-2h_k\\gamma_{e,k}}}.
```

The effective material is recovered from

```math
\\varepsilon_{r,eq}=1-
\\frac{\\Re\\{\\gamma_{e,eq}^2\\}}{\\omega^2\\mu_0\\varepsilon_0},
\\qquad
\\sigma_{eq}=\\frac{\\Im\\{\\gamma_{e,eq}^2\\}}{\\omega\\mu_0}.
```

**Reference.** H. Xue et al., “Generalized Formulation and Surge Analysis on
Overhead Lines: Impedance/Admittance of a Multi-Layer Earth,” *IEEE
Transactions on Power Delivery*, 36(6), 3834–3845, 2021.
DOI: 10.1109/TPWRD.2021.3049595.
"""
description(::Formula{:Xue2021}) =
    "Xue et al. 2021 equivalent propagation constant"

@inline function xue_gamma(
        resistivity,
        relative_permittivity,
        angular_frequency,
        vacuum_permittivity,
        vacuum_permeability,
        air_gamma_squared
)
    jω = complex(zero(angular_frequency), angular_frequency)
    conductivity = inv(resistivity)
    material_gamma_squared = jω * vacuum_permeability *
                             (conductivity + jω * vacuum_permittivity *
                              relative_permittivity)
    return sqrt(material_gamma_squared - air_gamma_squared)
end

"""
$(TYPEDSIGNATURES)

Reduce a horizontally layered, nonmagnetic earth to the equivalent propagation
constant proposed by Xue et al. The recursion is evaluated from the bottommost
soil layer to the surface. Because the transverse layer constant is defined by
``\\gamma_{e,k}^2 = \\gamma_k^2 - \\gamma_0^2``, relative permittivity is
reconstructed as

```math
\\epsilon_{r,\\mathrm{eq}} = 1 -
\\frac{\\Re\\{\\gamma_{e,\\mathrm{eq}}^2\\}}
     {\\omega^2\\mu_0\\epsilon_0}, \\qquad
\\sigma_{\\mathrm{eq}} =
\\frac{\\Im\\{\\gamma_{e,\\mathrm{eq}}^2\\}}{\\omega\\mu_0}.
```

The formulation defines equivalent conductivity and permittivity. Relative
permeability is inherited from the selected reconstruction layer and must be
one for every soil layer.

# Notes

This registered route supports overhead conductor pairs, matching the scope of
H. Xue et al., “Generalized Formulation and Surge Analysis on Overhead Lines:
Impedance/Admittance of a Multi-Layer Earth,” IEEE Transactions on Power
Delivery, 36(6), 3834–3845, 2021. DOI: 10.1109/TPWRD.2021.3049595.
"""
function xue(
        ::Val{:overhead},
        rho::AbstractVector{T},
        eps_r::AbstractVector{T},
        mu_r::AbstractVector{T},
        model::EarthModel{T},
        pair,
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    _horizontal(model, :Xue2021)
    _nonmagnetic(mu_r, :Xue2021)
    base = _material(rho, eps_r, mu_r, model, values.layer)
    unit = one(frequency)
    epsilon0 = unit * 88541878128 * (unit * 10)^(-22)
    mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
    omega = 2 * (unit * π) * frequency
    air_gamma_squared = -mu0 * epsilon0 * omega^2

    bottom = lastindex(rho)
    gamma_equivalent = xue_gamma(
        rho[bottom], eps_r[bottom], omega, epsilon0, mu0, air_gamma_squared
    )
    @inbounds for layer in (bottom - 1):-1:2
        gamma_top = xue_gamma(
            rho[layer], eps_r[layer], omega, epsilon0, mu0, air_gamma_squared
        )
        decay = exp(-2 * model.layers[layer].thickness * gamma_top)
        difference = gamma_top - gamma_equivalent
        sum = gamma_top + gamma_equivalent
        gamma_equivalent = gamma_top *
                           (sum - difference * decay) /
                           (sum + difference * decay)
    end

    squared = gamma_equivalent^2
    conductivity = imag(squared) / (omega * mu0)
    relative_permittivity = one(frequency) -
                            real(squared) / (omega^2 * mu0 * epsilon0)
    return EarthMaterial(inv(conductivity), relative_permittivity, base.mu_r)
end

@inline function (functor::Functor{:Xue2021})(
        layout::Val{:overhead},
        rho::AbstractVector{T},
        eps_r::AbstractVector{T},
        mu_r::AbstractVector{T},
        model::EarthModel{T},
        pair,
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    return xue(layout, rho, eps_r, mu_r, model, pair, frequency, values)
end

# Return the stable discovery identifier.
:Xue2021
