function routes(::Val{:Papadopoulos2009})
    (
        self = papadopoulosetaly2009,
        mutual = papadopoulosetaly2009,
        Γ = papadopoulosetaly2009_gamma
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
    "Papadopoulos et al. two-layer overhead potential coefficient (2009)"
end

function papadopoulosetaly2009_gamma(jω, permeability, permittivity)
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
Evaluate the Papadopoulos et al. two-layer overhead potential coefficient:

```math
P_{e,ij}=\frac1{2\pi\varepsilon_0}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty(F_{ij}^{P}+G_{ij}^{P})\cos(y_{ij}\lambda)d\lambda\right].
```

```math
F_{ij}^{P}=\mu_1\frac{s_{12}+d_{12}e^{-2a_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2a_1d}}e^{-\lambda(h_i+h_j)},
```

```math
G_{ij}^{P}=\lambda\frac{
\mu_0\mu_1(\gamma_0^2-\gamma_1^2)(s_{12}+d_{12}e^{-2a_1d})
(S_{12}+D_{12}e^{-2a_1d})-
4\mu_0\mu_1^2\mu_2a_1^2\gamma_0^2(\gamma_2^2-\gamma_1^2)e^{-2a_1d}}
{(S_{01}S_{12}+D_{01}D_{12}e^{-2a_1d})
(s_{01}s_{12}+d_{01}d_{12}e^{-2a_1d})}e^{-\lambda(h_i+h_j)}.
```

Adjacent parenthesized factors are products, exactly as in the corpus.
"""
function papadopoulosetaly2009(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    integral = _quadrature(state) do lambda
        a0 = sqrt(lambda^2 + state.gamma_medium_squared[1] + state.gamma_squared)
        a1 = sqrt(lambda^2 + state.gamma_medium_squared[2] + state.gamma_squared)
        a2 = sqrt(lambda^2 + state.gamma_medium_squared[3] + state.gamma_squared)
        mu0, mu1, mu2 = state.mu
        gamma0, gamma1, gamma2 = state.gamma_medium_squared
        s01 = a0 * mu1 + a1 * mu0
        d01 = a0 * mu1 - a1 * mu0
        s12 = a1 * mu2 + a2 * mu1
        d12 = a1 * mu2 - a2 * mu1
        S01 = mu0 * gamma1 * a0 + mu1 * gamma0 * a1
        D01 = mu0 * gamma1 * a0 - mu1 * gamma0 * a1
        S12 = mu1 * gamma2 * a1 + mu2 * gamma1 * a2
        D12 = mu1 * gamma2 * a1 - mu2 * gamma1 * a2
        decay = exp(-2a1 * d)
        F = mu1 * (s12 + d12 * decay) /
            (s01 * s12 + d01 * d12 * decay)
        numerator = mu0 * mu1 * (gamma0 - gamma1) *
                    (s12 + d12 * decay) * (S12 + D12 * decay) -
                    4mu0 * mu1^2 * mu2 * a1^2 * gamma0 *
                    (gamma2 - gamma1) * decay
        denominator = (S01 * S12 + D01 * D12 * decay) *
                      (s01 * s12 + d01 * d12 * decay)
        G = lambda * numerator / denominator
        (F + G) * exp(-lambda * geometry.H) * cos(lambda * geometry.y_ij)
    end
    return (log(geometry.D_ij / geometry.d_ij) + 2 * integral) /
           (2π * state.epsilon[1])
end

:Papadopoulos2009
