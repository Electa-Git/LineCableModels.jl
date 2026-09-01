function routes(identifier::Val{:Carson1926})
    return (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Carson1926})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Carson1926}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Classical homogeneous, conductive-earth overhead
impedance. Displacement currents and longitudinal propagation are neglected.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{\\lambda+\\sqrt{\\lambda^2+\\gamma_g^2}}d\\lambda\\right],
\\qquad \\gamma_g^2=j\\omega\\mu_0\\sigma_g.
```

**Reference.** J. R. Carson, “Wave Propagation in Overhead Wires with Ground
Return,” *Bell System Technical Journal*, 5, 539–554, 1926.
"""
description(::Formula{:Carson1926}) = "Carson homogeneous-earth overhead impedance (1926)"

function propagation_constant(::Val{:Carson1926}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Carson1926})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Carson1926), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Carson's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty
\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{\lambda+\sqrt{\lambda^2+\gamma_g^2}}\,d\lambda\right],
```

where ``\gamma_g^2=j\omega\mu_0\sigma_g``. The original Carson
assumptions neglect earth displacement current, air displacement current in
the correction, and longitudinal propagation.

This same overhead leaf is reused by pair-complete recipes such as
`:Pollaczek1926`, `:Ametani2009`, and `:Lucca1994`; the numerical kernel is
defined only here.

# Reference

J. R. Carson, "Wave propagation in overhead wires with ground return,"
*Bell System Technical Journal*, vol. 5, pp. 539-554, 1926.
"""
function earth_impedance(
        ::Val{:Carson1926}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    gamma_squared = state.gamma_medium_squared[2]
    integral = _quadrature(state) do lambda
        attenuation = sqrt(lambda^2 + gamma_squared)
        exp(-geometry.H * lambda) * cos(geometry.y_ij * lambda) /
        (lambda + attenuation)
    end
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

:Carson1926
