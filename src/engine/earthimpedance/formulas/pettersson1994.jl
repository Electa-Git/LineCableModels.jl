routes(::Val{:Pettersson1994}) = (
    self = pettersson1994,
    mutual = pettersson1994,
    Γ = pettersson1994_gamma
)

function assumptions(::Val{:Pettersson1994})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Pettersson1994}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Wideband complex-image approximation for conductors above,
on, or below homogeneous earth.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+\\ln\\frac{
\\sqrt{(H+2/\\beta_\\gamma)^2+y_{ij}^2}}{D_{ij}}\\right],
\\qquad \\beta_\\gamma=\\sqrt{\\gamma_g^2-\\gamma_0^2}.
```

**Reference.** P. Pettersson, “Image Representation of Wave Propagation on
Wires Above, On and Under Ground,” *IEEE Transactions on Power Delivery*, 9,
1049–1055, 1994. DOI: 10.1109/61.296290.
"""
description(::Formula{:Pettersson1994}) = "Pettersson wideband image approximation (1994)"

pettersson1994_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Pettersson1994})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Pettersson1994), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Pettersson's wideband image approximation for overhead earth-return
impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{D_{ij}}{d_{ij}}+M_{e,ij}\right],
\qquad
M_{e,ij}=\ln\frac{
\sqrt{(H+2/\beta_\gamma)^2+y_{ij}^2}}{D_{ij}},
```

```math
\beta_\gamma=\sqrt{\gamma_g^2-\gamma_0^2},\qquad
\gamma_m^2=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m).
```

# Reference

P. Pettersson, "Image representation of wave propagation on wires above,
on and under ground," *IEEE Transactions on Power Delivery*, vol. 9,
pp. 1049-1055, 1994. DOI: 10.1109/61.296290.
"""
function pettersson1994(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    βγ = sqrt(
        state.gamma_medium_squared[2] - state.gamma_medium_squared[1]
    )
    image = sqrt((geometry.H + 2 / βγ)^2 + geometry.y_ij^2)
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) * log(image / geometry.d_ij)
end

:Pettersson1994
