function routes(::Val{:Vance1978})
    (
        self = vance1978_self,
        mutual = vance1978_self,
        Γ = vance1978_gamma
    )
end

function assumptions(::Val{:Vance1978})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Vance1978}) = Val(:zero)
description(::Formula{:Vance1978}) = "Vance lossy-cylinder underground self term (1978)"

vance1978_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Vance1978})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:Vance1978), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Vance's lossy cylindrical-dielectric self impedance:

```math
Z_{e,ii}=\frac{\omega\mu_0}{2\pi\gamma_1r_i}
\frac{H_0^{(1)}(j\gamma_1r_i)}{H_1^{(1)}(j\gamma_1r_i)},\qquad
\gamma_1^2=j\omega\mu_0\sigma_1.
```

The same radial kernel supplies mutual terms by replacing ``r_i`` with the
horizontal pair separation; self evaluation is the standard
``y_{ii}=r_i`` specialization.
"""
function vance1978_self(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    radius = pair.separation
    argument = im * state.gamma[2] * radius
    omega = imag(state.jω)
    return omega * state.mu[1] / (2π * state.gamma[2] * radius) *
           hankelh1(0, argument) / hankelh1(1, argument)
end

:Vance1978
