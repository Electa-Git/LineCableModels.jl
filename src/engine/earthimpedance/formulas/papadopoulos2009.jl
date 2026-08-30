function routes(::Val{:Papadopoulos2009})
    (
        self = papadopoulosetal2009,
        mutual = papadopoulosetal2009,
        Γ = papadopoulosetal2009_gamma
    )
end

function assumptions(::Val{:Papadopoulos2009})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2009}) = Val(:explicit)
media(::Formula{:Papadopoulos2009}) = Val(:stratified)
function description(::Formula{:Papadopoulos2009})
    "Papadopoulos et al. general two-layer overhead impedance (2009)"
end

function papadopoulosetal2009_gamma(jω, permeability, permittivity)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2009})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) == 3 || throw(DimensionMismatch(
        ":Papadopoulos2009 requires air and exactly two earth layers"
    ))
    k_0 = Γ === nothing ? sqrt(oftype(jω, (-jω^2) * mu[1] * epsilon[1])) : Γ
    return _stratified_functor(
        Val(:Papadopoulos2009), formula,
        rho, epsilon, mu, jω, k_0, segments, thickness
    )
end

raw"""
Evaluate the Papadopoulos et al. general two-layer overhead kernel:

```math
F_{ij}^{P}=\mu_1\frac{s_{12}+d_{12}e^{-2a_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2a_1d}}e^{-\lambda(h_i+h_j)},
```

```math
s_{mn}=a_m\mu_n+a_n\mu_m,\qquad
d_{mn}=a_m\mu_n-a_n\mu_m,\qquad
a_m=\sqrt{\lambda^2+\gamma_m^2+k_0^2}.
```

The contribution is
``Z_{e,ij}=(j\omega\mu_0/(2\pi))[\ln(D_{ij}/d_{ij})+2\int F_{ij}^P
\cos(y_{ij}\lambda)d\lambda]``. By default ``k_0`` is the selected `Γ`.
"""
function papadopoulosetal2009(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    integral = _quadrature(state) do lambda
        a0 = sqrt(lambda^2 + state.gamma_medium_squared[1] + state.gamma_squared)
        a1 = sqrt(lambda^2 + state.gamma_medium_squared[2] + state.gamma_squared)
        a2 = sqrt(lambda^2 + state.gamma_medium_squared[3] + state.gamma_squared)
        s01 = a0 * state.mu[2] + a1 * state.mu[1]
        d01 = a0 * state.mu[2] - a1 * state.mu[1]
        s12 = a1 * state.mu[3] + a2 * state.mu[2]
        d12 = a1 * state.mu[3] - a2 * state.mu[2]
        decay = exp(-2a1 * d)
        state.mu[2] * (s12 + d12 * decay) /
        (s01 * s12 + d01 * d12 * decay) *
        exp(-lambda * geometry.H) * cos(lambda * geometry.y_ij)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

:Papadopoulos2009
