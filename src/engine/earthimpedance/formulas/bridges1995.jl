function routes(::Val{:Bridges1995})
    (
        self = bridges1995_self,
        mutual = bridges1995_self,
        Γ = bridges1995_gamma
    )
end

function assumptions(::Val{:Bridges1995})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Bridges1995}) = Val(:zero)
description(::Formula{:Bridges1995}) = "Bridges buried single-cable approximation (1995)"

bridges1995_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Bridges1995})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Bridges1995), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Bridges' buried single-cable impedance:

```math
Z_{e,ii}=-\frac{j\omega\mu_0}{2\pi}
\ln\left(\frac{e_c\gamma_1r_i}{2}\right),\qquad
e_c=1.7811,\quad \gamma_1^2=j\omega\mu_0\sigma_1.
```

The same radial kernel supplies mutual terms by replacing ``r_i`` with the
horizontal pair separation; self evaluation is the standard
``y_{ii}=r_i`` specialization.
"""
function bridges1995_self(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    e_c = oftype(pair.separation, 1.7811)
    return -state.jω * state.mu[1] / (2π) *
           log(e_c * state.gamma[2] * pair.separation / 2)
end

:Bridges1995
