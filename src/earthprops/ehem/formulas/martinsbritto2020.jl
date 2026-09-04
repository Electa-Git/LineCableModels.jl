"Return the reconstruction assumptions of the Martins–Britto et al. EHEM."
assumptions(::Val{:MartinsBritto2020}) = (layer = -1,)

"""
$(TYPEDSIGNATURES)

**Identification.** Equivalent real conductivity; relative permittivity and
permeability are inherited from the selected reconstruction layer.

**Expression.** For layer ``k`` of thickness ``h_k``,

```math
\\sigma_{eq,k}=\\sigma_k\\left[
\\frac{\\sqrt{\\sigma_k}+\\sqrt{\\sigma_{eq,k+1}}-
(\\sqrt{\\sigma_k}-\\sqrt{\\sigma_{eq,k+1}})
e^{-2h_k\\sqrt{\\pi f\\mu_0\\sigma_k}}}
{\\sqrt{\\sigma_k}+\\sqrt{\\sigma_{eq,k+1}}+
(\\sqrt{\\sigma_k}-\\sqrt{\\sigma_{eq,k+1}})
e^{-2h_k\\sqrt{\\pi f\\mu_0\\sigma_k}}}
\\right]^2.
```

**Reference.** A. G. Martins-Britto, F. V. Lopes, and S. R. M. J.
Rondineau, “Multilayer Earth Structure Approximation by a Homogeneous
Conductivity Soil for Ground Return Impedance Calculations,” *IEEE
Transactions on Power Delivery*, 35(2), 881–891, 2020.
DOI: 10.1109/TPWRD.2019.2930406.
"""
description(::Formula{:MartinsBritto2020}) =
    "Martins–Britto et al. equivalent-conductivity earth reduction (2020)"

@inline function equivalent_conductivity_step(
        ::Val{:MartinsBritto2020},
        conductivity_top,
        conductivity_bottom,
        thickness,
        frequency,
        vacuum_permeability
)
    root_top = sqrt(conductivity_top)
    root_bottom = sqrt(conductivity_bottom)
    decay = exp(
        -2 * thickness * sqrt(
            (one(frequency) * π) * frequency * vacuum_permeability * conductivity_top
        )
    )
    difference = root_top - root_bottom
    sum = root_top + root_bottom
    ratio = (sum - difference * decay) / (sum + difference * decay)
    return conductivity_top * ratio^2
end

"""
$(TYPEDSIGNATURES)

Reduce a horizontally layered, nonmagnetic earth to the real equivalent
conductivity proposed by Martins–Britto et al. The recursion is evaluated from
the bottommost soil layer to the surface:

```math
\\sigma_{\\mathrm{eq},k} = \\sigma_k
\\left[
\\frac{\\sqrt{\\sigma_k}+\\sqrt{\\sigma_{\\mathrm{eq},k+1}}-
      (\\sqrt{\\sigma_k}-\\sqrt{\\sigma_{\\mathrm{eq},k+1}})
      e^{-2h_k\\sqrt{\\pi f\\mu_0\\sigma_k}}}
     {\\sqrt{\\sigma_k}+\\sqrt{\\sigma_{\\mathrm{eq},k+1}}+
      (\\sqrt{\\sigma_k}-\\sqrt{\\sigma_{\\mathrm{eq},k+1}})
      e^{-2h_k\\sqrt{\\pi f\\mu_0\\sigma_k}}}
\\right]^2.
```

The formulation defines conductivity only. Relative permittivity and
permeability are inherited unchanged from the selected reconstruction layer,
which is the bottommost layer by default.

# Notes

This registered route supports overhead conductor pairs, matching the scope of
A. G. Martins–Britto, F. V. Lopes, and S. R. M. J. Rondineau, “Multilayer Earth
Structure Approximation by a Homogeneous Conductivity Soil for Ground Return
Impedance Calculations,” IEEE Transactions on Power Delivery, 35(2),
881–891, 2020. DOI: 10.1109/TPWRD.2019.2930406.
"""
function equivalent_material(
        ::Val{:MartinsBritto2020},
        ::Val{:overhead},
        rho::AbstractVector{T},
        eps_r::AbstractVector{T},
        mu_r::AbstractVector{T},
        model::EarthModel{T},
        pair,
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    _horizontal(model, :MartinsBritto2020)
    _nonmagnetic(mu_r, :MartinsBritto2020)
    base = _material(rho, eps_r, mu_r, model, values.layer)
    unit = one(frequency)
    mu0 = unit * 4 * (unit * π) * (unit * 10)^(-7)
    bottom = lastindex(rho)
    conductivity_equivalent = inv(rho[bottom])
    @inbounds for layer in (bottom - 1):-1:2
        conductivity_equivalent = equivalent_conductivity_step(
            Val(:MartinsBritto2020),
            inv(rho[layer]),
            conductivity_equivalent,
            model.layers[layer].thickness,
            frequency,
            mu0
        )
    end
    return EarthMaterial(
        inv(conductivity_equivalent), base.eps_r, base.mu_r
    )
end

# Return the stable discovery identifier.
:MartinsBritto2020
