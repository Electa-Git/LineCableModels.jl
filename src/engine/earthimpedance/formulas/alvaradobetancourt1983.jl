function routes(identifier::Val{:Alvarado1983})
    return (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Alvarado1983})
    (
        air = _lossless,
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Alvarado1983}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Thin overhead self or mutual conductor pair with height sum H and horizontal separation x. |
| Calculated quantities | Closed-form overhead earth-return impedance |
| Earth structure | Homogeneous conducting nonmagnetic earth below lossless air. |
| Model and approximation | Complex-image approximation with a cubic residual correction; conduction-only earth. |
| Main source | F. L. Alvarado and R. Betancourt (1983) |
| Citation key(s) | Primary: `:Alvarado1983`; equation witness: `:Papadopoulos2020` |
| Evidence status | Published comparison page image checked, equation (8) and the stated self substitution. |

**Description.** Closed-form approximation to Carson's overhead integral
for conductive, nonmagnetic earth.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\left[\\ln\\frac{D_{ij}}{d_{ij}}+J_{m,ij}\\right],
```

```math
\\begin{aligned}
J_{m,ij}={}&\\frac12\\ln\\frac{[1+p_g/(H/2)]^2+(y_{ij}/H)^2}
{1+(y_{ij}/H)^2}\\\\
&-\\frac1{24}\\sum_{s\\in\\{-1,1\\}}
\\left[1+\\frac{H}{2p_g}\\left(1+sj\\frac{y_{ij}}H\\right)\\right]^{-3},
\\qquad p_g=(j\\omega\\mu_0\\sigma_g)^{-1/2}.
\\end{aligned}
```

**Reference.** F. L. Alvarado and R. Betancourt, “An Accurate Closed-Form
Approximation for Ground Return Impedance Calculations,” *Proceedings of the
IEEE*, 71, 279–280, 1983. DOI: 10.1109/PROC.1983.12573.
"""
function description(::Formula{:Alvarado1983})
    "Alvarado-Betancourt closed-form Carson correction (1983)"
end

function propagation_constant(
        ::Val{:Alvarado1983}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Alvarado1983})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Alvarado1983), formula,
        rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Alvarado and Betancourt's closed-form approximation to Carson's
overhead earth-return correction:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[\ln\frac{D_{ij}}{d_{ij}}+J_{m,ij}\right],
```

```math
\begin{aligned}
J_{m,ij}={}&\frac12\ln\left[
\frac{\left(1+p_g/(H/2)\right)^2+(y_{ij}/H)^2}
{1+(y_{ij}/H)^2}\right]\\
&-\frac1{24}\left\{
\left[1+\frac{H}{2p_g}\left(1+j\frac{y_{ij}}H\right)\right]^{-3}
+\left[1+\frac{H}{2p_g}\left(1-j\frac{y_{ij}}H\right)\right]^{-3}
\right\},
\end{aligned}
```

where ``p_g=[j\omega\mu_0\sigma_g]^{-1/2}``, ``H=h_i+h_j``, and
``d_{ij}`` and ``D_{ij}`` are the direct and image distances.

# Reference

F. L. Alvarado and R. Betancourt, "An accurate closed-form approximation
for ground return impedance calculations," *Proceedings of the IEEE*, vol.
71, pp. 279-280, 1983. DOI: 10.1109/PROC.1983.12573.
"""
function earth_impedance(
        ::Val{:Alvarado1983}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    p_g = inv(state.gamma[2])
    horizontal = pair.row == pair.column ? zero(geometry.H) : geometry.y_ij
    horizontal_ratio = horizontal / geometry.H
    depth_ratio = geometry.H / (2p_g)
    logarithmic = log(
        ((1 + inv(depth_ratio))^2 + horizontal_ratio^2) /
        (1 + horizontal_ratio^2)
    ) / 2
    positive = 1 + depth_ratio * (1 + im * horizontal_ratio)
    negative = 1 + depth_ratio * (1 - im * horizontal_ratio)
    correction = logarithmic - (inv(positive^3) + inv(negative^3)) / 24
    πT = one(geometry.H) * π
    return state.jω * state.mu[1] / (2πT) *
           (log(hypot(horizontal, geometry.H) / geometry.d_ij) + correction)
end

:Alvarado1983
