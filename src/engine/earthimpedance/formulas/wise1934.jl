routes(::Val{:Wise1934}) = (
    self = wise1934,
    mutual = wise1934,
    Γ = wise1934_gamma
)

assumptions(::Val{:Wise1934}) = (
    air = _full,
    earth = _full,
    permeability = _material
)

propagation(::Val{:Wise1934}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Homogeneous-earth wideband overhead integral retaining
earth displacement current and magnetic permeability.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty
\\frac{\\mu_1e^{-\\lambda H}}
{\\lambda\\mu_1+a_1\\mu_0}\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\quad a_1=\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}.
```

**Reference.** W. H. Wise, “Propagation of High-Frequency Currents in Ground
Return Circuits,” *Proceedings of the Institute of Radio Engineers*, 22,
522–527, 1934.
"""
description(::Formula{:Wise1934}) = "Wise homogeneous-earth overhead impedance (1934)"

wise1934_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Wise1934})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:Wise1934), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Wise's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty F_{ij}^{W}(\lambda)\cos(y_{ij}\lambda)\,d\lambda\right],
```

```math
F_{ij}^{W}=\frac{\mu_1e^{-\lambda(h_i+h_j)}}
{\lambda\mu_1+a_1\mu_0},\qquad
a_1=\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}.
```
"""
function wise1934(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    contrast = state.gamma_medium_squared[2] - state.gamma_medium_squared[1]
    integral = _quadrature(state) do lambda
        a_1 = sqrt(lambda^2 + contrast)
        state.mu[2] * exp(-lambda * geometry.H) * cos(lambda * geometry.y_ij) /
        (lambda * state.mu[2] + a_1 * state.mu[1])
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

:Wise1934
