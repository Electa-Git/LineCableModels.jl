function routes(identifier::Val{:Noda2006})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

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

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Straight parallel cylindrical conductors with radii ``r_i``, heights ``h_i``, and horizontal separation ``x_{ij}``; no insulation layer. |
| Calculated quantities | Per-unit-length self and mutual ground-return impedances for parallel overhead cylindrical conductors |
| Earth structure | Homogeneous lossy half-space below a plane interface. |
| Model and approximation | Noda begins from Carson's integral (2), uses a two-logarithm form, enforces ``A+B=1`` and the infinite-frequency condition ``A\\alpha+B\\beta=1``, then deliberately omits the exact zero-frequency matching condition because the minimax fit is better between the limits. Piecewise-linear ``A(\\theta)`` and ``\\alpha(\\theta)`` are least-squares fits to minimax solutions; ``B=1-A`` and ``\\beta=(1-A\\alpha)/(1-A)`` follow analytically. The fitted formula is an approximation, not a full-wave exact result. |
| Main source | Taku Noda (2006) |
| Citation key(s) | `:Noda2006` |
| Evidence status | Original PDF page images checked |

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

function propagation_constant(::Val{:Noda2006}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

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
function earth_impedance(
        ::Val{:Noda2006}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    R = typeof(geometry.H)
    lateral=pair.row==pair.column ? zero(R) : geometry.y_ij
    θ = atan(lateral, geometry.H) * R(180) / R(π)
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
    S_a = sqrt((geometry.H + 2a * p_g)^2 + lateral^2)
    S_β = sqrt((geometry.H + 2β * p_g)^2 + lateral^2)
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) *
           (A*log(S_a/geometry.d_ij)+(1-A)*log(S_β/geometry.d_ij))
end

:Noda2006
