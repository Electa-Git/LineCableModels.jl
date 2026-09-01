function routes(identifier::Val{:Gary1976})
    return (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Gary1976})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Gary1976}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Complex-depth logarithmic approximation for overhead
conductors.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\ln\\frac{S_{ij}}{d_{ij}},\\qquad
S_{ij}=\\sqrt{(H+2h_e)^2+y_{ij}^2},\\qquad
h_e=(j\\omega\\mu_0\\sigma_1)^{-1/2}.
```

**Reference.** C. Gary, “Approche complète de la propagation multifilaire en
haute fréquence par utilisation des matrices complexes,” *EDF Bulletin de la
Direction des Études et Recherches*, série B, 1976; formula as reproduced in
Ametani et al., IET, 2021.
"""
description(::Formula{:Gary1976}) = "Gary complex-depth approximation (1976)"

function propagation_constant(::Val{:Gary1976}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Gary1976})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:Gary1976), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Gary's complex-depth overhead earth-return impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\ln\frac{S_{ij}}{d_{ij}},\qquad
S_{ij}=\sqrt{(h_i+h_j+2h_e)^2+y_{ij}^2},\qquad
h_e=\frac{1}{\sqrt{j\omega\mu_0\sigma_1}}.
```

The direct distance is
``d_{ij}=\sqrt{(h_i-h_j)^2+y_{ij}^2}``. For self terms, `pair`
supplies the outer radius as ``y_{ii}=r_i``.
"""
function earth_impedance(
        ::Val{:Gary1976}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    h_e = inv(state.gamma[2])
    S_ij = sqrt((geometry.H + 2h_e)^2 + geometry.y_ij^2)
    return state.jω * state.mu[1] / (2π) * log(S_ij / geometry.d_ij)
end

:Gary1976
