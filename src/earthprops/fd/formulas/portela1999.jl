"Return the fitted dispersion assumptions of the Portela relation."
assumptions(::Val{:Portela1999}) = (beta = 0.1, exponent = 0.72)

"""
$(TYPEDSIGNATURES)

**Identification.** Causal power-law soil dispersion.

**Expression.** For ``\\omega=2\\pi f``, ``\\beta=0.1``, and ``\\alpha=0.72``,

```math
\\sigma(f)=\\sigma_0+\\beta10^{-6}\\omega^\\alpha,
\\qquad
\\varepsilon_r(f)=\\frac{\\beta10^{-6}\\tan(\\pi\\alpha/2)
\\omega^{\\alpha-1}}{\\varepsilon_0}.
```

**Reference.** C. M. Portela, “Measurement and Modeling of Soil
Electromagnetic Behavior,” *IEEE International Symposium on Electromagnetic
Compatibility*, 1004–1009, 1999.
"""
description(::Formula{:Portela1999}) = "Portela 1999 (power-law soil dispersion)"

#=
Evaluate the Portela power-law soil-dispersion relation.

The legacy μS/cm calculation is expressed directly in SI: the fitted `beta`
coefficient contributes `beta × 10⁻⁶ × ω^exponent` S/m. Returned resistivity
uses Ω·m and relative permittivity is dimensionless.

# Reference

C. M. Portela, *Measurement and Modeling of Soil Electromagnetic Behavior*,
IEEE International Symposium on Electromagnetic Compatibility, 1999,
pp. 1004–1009.
=#
function (::Functor{:Portela1999})(
        material::EarthMaterial{T},
        frequency::T,
        values::NamedTuple
) where {T <: Real}
    beta = convert(T, values.beta)
    exponent = convert(T, values.exponent)
    angular_frequency = 2 * (one(frequency) * π) * frequency
    fitted_scale = beta * convert(T, 1e-6)
    conductivity = inv(material.rho) + fitted_scale * angular_frequency^exponent
    relative_permittivity = fitted_scale *
                            tan((one(frequency) * π) * exponent / 2) *
                            angular_frequency^(exponent - one(exponent)) /
                            vacuum_permittivity(frequency)
    return EarthMaterial{T}(inv(conductivity), relative_permittivity, material.mu_r)
end

# Return the stable discovery identifier.
:Portela1999
