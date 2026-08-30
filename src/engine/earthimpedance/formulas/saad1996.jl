function routes(::Val{:Saad1996})
    (
        self = saad1996,
        mutual = saad1996,
        Γ = saad1996_gamma
    )
end

function assumptions(::Val{:Saad1996})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Saad1996}) = Val(:zero)
description(::Formula{:Saad1996}) = "Saad underground closed form (1996)"

saad1996_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Saad1996})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Saad1996), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate the Saad et al. underground approximation:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\gamma_1R_{ab})+
\frac{2e^{-(h_i+h_j)\gamma_1}}{4+\gamma_1^2R_{ab}^2}\right].
```

For a self term ``R_{ab}=r_i``; for a mutual term it is the horizontal
center separation supplied by `pair`.

# Reference

O. Saad, G. Gaba, and M. Giroux, "A closed-form approximation for ground
return impedance of underground cables," *IEEE Transactions on Power
Delivery*, vol. 11, no. 3, pp. 1536-1545, 1996.
"""
function saad1996(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    radius = geometry.y_ij
    correction = 2exp(-geometry.H * gamma) / (4 + gamma^2 * radius^2)
    direct = _complex_result(
        state.jω, special_besselk(0, gamma * radius)
    )
    πT = one(radius) * π
    return state.jω * state.mu[1] / (2πT) * (direct + correction)
end

:Saad1996
