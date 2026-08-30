routes(::Val{:Noda2006}) = (
    self = noda2006,
    mutual = noda2006,
    Γ = noda2006_gamma
)

function assumptions(::Val{:Noda2006})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Noda2006}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Double-logarithmic approximation to Carson's overhead
integral.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+A\\ln\\frac{S_a}{D_{ij}}+
(1-A)\\ln\\frac{S_\\beta}{D_{ij}}\\right],
```

```math
S_a=\\sqrt{(H+2ap_g)^2+y_{ij}^2},\\quad
S_\\beta=\\sqrt{(H+2\\beta p_g)^2+y_{ij}^2},\\quad
\\beta=\\frac{1-Aa}{1-A}.
```

The piecewise ``A`` and ``a`` coefficients use
``\\theta=\\tan^{-1}(y_{ij}/H)`` as specified in the paper.

**Reference.** T. Noda, “A Double Logarithmic Approximation of Carson's
Ground-Return Impedance,” *IEEE Transactions on Power Delivery*, 21,
472–479, 2006. DOI: 10.1109/TPWRD.2005.852307.
"""
description(::Formula{:Noda2006}) = "Noda double-logarithmic approximation (2006)"

noda2006_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Noda2006})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Noda2006), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Noda's double-logarithmic overhead earth-return approximation:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[
\ln\frac{D_{ij}}{d_{ij}}
+A\ln\frac{S_a}{D_{ij}}
+(1-A)\ln\frac{S_\beta}{D_{ij}}\right],
```

```math
S_a=\sqrt{(H+2ap_g)^2+y_{ij}^2},\qquad
S_\beta=\sqrt{(H+2\beta p_g)^2+y_{ij}^2},\qquad
\beta=\frac{1-Aa}{1-A},
```

with ``p_g=[j\omega\mu_0\sigma_g]^{-1/2}`` and
``\theta=\tan^{-1}(y_{ij}/H)`` in degrees. The empirical coefficients are

```math
A=\begin{cases}
0.07360,&\theta\le50.45^\circ,\\
0.002474\theta-0.05127,&\theta>50.45^\circ,
\end{cases}\qquad
a=\begin{cases}
0.1500,&\theta\le50.45^\circ,\\
0.004726\theta-0.08852,&\theta>50.45^\circ.
\end{cases}
```

# Reference

T. Noda, "A double logarithmic approximation of Carson's ground-return
impedance," *IEEE Transactions on Power Delivery*, vol. 21,
pp. 472-479, 2006. DOI: 10.1109/TPWRD.2005.852307.
"""
function noda2006(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    R = typeof(geometry.H)
    θ = atan(geometry.y_ij, geometry.H) * R(180) / R(π)
    threshold = R(50.45)
    if θ <= threshold
        A = R(0.07360)
        a = R(0.1500)
    else
        A = R(0.002474) * θ - R(0.05127)
        a = R(0.004726) * θ - R(0.08852)
    end
    β = (1 - A * a) / (1 - A)
    p_g = inv(state.gamma[2])
    S_a = sqrt((geometry.H + 2a * p_g)^2 + geometry.y_ij^2)
    S_β = sqrt((geometry.H + 2β * p_g)^2 + geometry.y_ij^2)
    correction = A * log(S_a / geometry.D_ij) +
                 (1 - A) * log(S_β / geometry.D_ij)
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) *
           (log(geometry.D_ij / geometry.d_ij) + correction)
end

:Noda2006
