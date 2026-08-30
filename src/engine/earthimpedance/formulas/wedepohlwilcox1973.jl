function routes(::Val{:WedepohlWilcox1973})
    (
        self = wedepohlwilcox1973_self,
        mutual = wedepohlwilcox1973_mutual,
        Γ = wedepohlwilcox1973_gamma
    )
end

function assumptions(::Val{:WedepohlWilcox1973})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:WedepohlWilcox1973}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Low-frequency underground expansion for conductive,
nonmagnetic earth.

**Expression.**

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1r_i}{2}\\right)+\\frac12-
\\frac43\\gamma_1h_i\\right],
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[-\\ln
\\left(\\frac{e_c\\gamma_1y_{ij}}{2}\\right)+\\frac12-
\\frac23\\gamma_1H\\right],\\qquad e_c=1.7811.
```

**Reference.** L. M. Wedepohl and D. J. Wilcox, “Transient Analysis of
Underground Power-Transmission Systems: System-Model and Wave-Propagation
Characteristics,” *Proceedings of the IEE*, 120, 253–260, 1973.
"""
function description(::Formula{:WedepohlWilcox1973})
    "Wedepohl-Wilcox low-frequency underground approximation (1973)"
end

function wedepohlwilcox1973_gamma(jω, permeability, permittivity)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:WedepohlWilcox1973})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:WedepohlWilcox1973), formula,
        rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Wedepohl-Wilcox low-frequency underground terms:

```math
Z_{e,ii}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1r_i}{2}\right)+\frac12
-\frac43\gamma_1h_i\right],
```

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[-\ln\left(\frac{e_c\gamma_1y_{ij}}{2}\right)+\frac12
-\frac23\gamma_1(h_i+h_j)\right],\qquad e_c=1.7811.
```
"""
function wedepohlwilcox1973_self(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    e_c = oftype(geometry.h_i, 1.7811)
    bracket = -log(e_c * state.gamma[2] * geometry.y_ij / 2) +
              one(geometry.h_i) / 2 -
              (4 * one(geometry.h_i) / 3) * state.gamma[2] * geometry.h_i
    return state.jω * state.mu[1] / (2π) * bracket
end

function wedepohlwilcox1973_mutual(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    e_c = oftype(geometry.H, 1.7811)
    bracket = -log(e_c * state.gamma[2] * geometry.y_ij / 2) +
              one(geometry.H) / 2 -
              (2 * one(geometry.H) / 3) * state.gamma[2] * geometry.H
    return state.jω * state.mu[1] / (2π) * bracket
end

:WedepohlWilcox1973
