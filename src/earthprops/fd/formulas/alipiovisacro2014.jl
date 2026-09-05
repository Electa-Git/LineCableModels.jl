"Return the fitted assumptions of the Alipio–Visacro causal soil model."
assumptions(::Val{:AlipioVisacro2014}) = (
    exponent = 0.54,
    epsilon_infinity = 12.0,
    scale = 1.26,
    conductivity_exponent = -0.73
)

"""
$(TYPEDSIGNATURES)

**Identification.** Causal power-law soil model fitted from measured soil
dispersion.

**Expression.** Let ``\\widehat\\sigma_0=1000/\\rho_0`` in mS/m,
``\\gamma=0.54``, and ``D=1.26\\widehat\\sigma_0^{-0.73}``. The implemented
relation is

```math
\\widehat\\sigma(f)=\\widehat\\sigma_0
\\left[1+D\\left(\\frac{f}{10^6}\\right)^\\gamma\\right],
```

```math
\\varepsilon_r(f)=12+
\\tan\\left(\\frac{\\pi\\gamma}{2}\\right)
\\frac{10^{-3}\\widehat\\sigma_0D f^{\\gamma-1}}
{2\\pi\\varepsilon_0\\cdot10^{6\\gamma}},\\qquad
\\sigma(f)=10^{-3}\\widehat\\sigma(f).
```

**Reference.** R. Alipio and S. Visacro, “Modeling the Frequency Dependence
of Electrical Parameters of Soil,” *IEEE Transactions on Electromagnetic
Compatibility*, 56(5), 2014. DOI: 10.1109/TEMC.2014.2313977.
"""
description(::Formula{:AlipioVisacro2014}) =
    "Alipio–Visacro causal soil dispersion (2014)"

#=
Evaluate the Alipio–Visacro causal soil relation.

The reference conductivity is converted from S/m to mS/m for the published
fit. The returned resistivity is converted back to Ω·m and relative
permittivity is dimensionless.

# Reference

R. Alipio and S. Visacro, *Modeling the Frequency Dependence of Electrical
Parameters of Soil*, IEEE Transactions on Electromagnetic Compatibility,
56(5), 2014. DOI: 10.1109/TEMC.2014.2313977.
=#
function earth_material(
        ::Val{:AlipioVisacro2014},
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    gamma = convert(T, values.exponent)
    epsilon_infinity = convert(T, values.epsilon_infinity)
    scale = convert(T, values.scale)
    conductivity_exponent = convert(T, values.conductivity_exponent)
    thousand = convert(T, 1000)
    million = convert(T, 1e6)
    conductivity_reference = thousand / material.rho
    dispersion = scale * conductivity_reference^conductivity_exponent
    epsilon0 = vacuum_permittivity(frequency)
    relative_permittivity = epsilon_infinity +
                            tan((one(frequency) * π) * gamma / 2) *
                            convert(T, 1e-3) * conductivity_reference * dispersion *
                            frequency^(gamma - one(gamma)) /
                            (
                                2 * (one(frequency) * π) * epsilon0 *
                                convert(T, 10)^(6 * gamma)
                            )
    conductivity = conductivity_reference *
                   (one(frequency) + dispersion * (frequency / million)^gamma)
    return EarthMaterial{T}(
        thousand / conductivity,
        relative_permittivity,
        material.mu_r
    )
end

# Return the stable discovery identifier.
:AlipioVisacro2014
