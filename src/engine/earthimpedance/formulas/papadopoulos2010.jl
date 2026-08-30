function routes(::Val{:Papadopoulos2010})
    (
        self = papadopoulosetal2010,
        mutual = papadopoulosetal2010,
        Γ = papadopoulosetal2010_gamma
    )
end

function assumptions(::Val{:Papadopoulos2010})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2010}) = Val(:explicit)
function description(::Formula{:Papadopoulos2010})
    "Papadopoulos et al. homogeneous-earth underground impedance (2010)"
end

function papadopoulosetal2010_gamma(jω, permeability, permittivity)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2010})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Papadopoulos2010), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Papadopoulos et al. homogeneous-earth underground impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}(\Delta_1+2S),
```

```math
\Delta_1=\int_0^\infty\frac{e^{-|h_i-h_j|\alpha_1}-e^{-(h_i+h_j)\alpha_1}}
{\alpha_1}\cos(y_{ij}\lambda)d\lambda,\qquad
S=\int_0^\infty\frac{e^{-(h_i+h_j)\alpha_1}}
{\alpha_0+\alpha_1}\cos(y_{ij}\lambda)d\lambda,
```

where ``\alpha_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}`` and
``k_x=\omega\sqrt{\mu_0\varepsilon_1}`` by default. An explicit `Γ`
replaces ``k_x`` without changing the leaf routes.
"""
function papadopoulosetal2010(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    alpha_squared = (
        state.gamma_medium_squared[1] + state.gamma_squared,
        state.gamma_medium_squared[2] + state.gamma_squared
    )
    radial = sqrt(alpha_squared[2])
    delta_1 = special_besselk(0, radial * geometry.d_ij) -
              special_besselk(0, radial * geometry.D_ij)
    S = _quadrature(state) do lambda
        alpha_0 = sqrt(lambda^2 + alpha_squared[1])
        alpha_1 = sqrt(lambda^2 + alpha_squared[2])
        exp(-geometry.H * alpha_1) /
        (alpha_0 + alpha_1) * cos(geometry.y_ij * lambda)
    end
    return state.jω * state.mu[1] / (2π) * (delta_1 + 2 * S)
end

:Papadopoulos2010
