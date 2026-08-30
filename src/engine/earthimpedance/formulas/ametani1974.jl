function routes(::Val{:Ametani1974})
    (
        self = ametanitwolayer1974,
        mutual = ametanitwolayer1974,
        Γ = ametanitwolayer1974_gamma
    )
end

function assumptions(::Val{:Ametani1974})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Ametani1974}) = Val(:zero)
media(::Formula{:Ametani1974}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

**Identification.** Overhead-conductor impedance above two magnetic earth
layers.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty F_{ij}^{A}(\\lambda)
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
F_{ij}^{A}=\\frac{b_1+b_2+(b_1-b_2)e^{-2a_1d}}
{(\\lambda+\\mu_0b_1)(b_1+b_2)+(\\lambda-\\mu_0b_1)
(b_1-b_2)e^{-2a_1d}}e^{-\\lambda H},
\\quad b_m=\\frac{a_m}{\\mu_m},\\quad
a_m=\\sqrt{\\lambda^2+\\gamma_m^2-\\gamma_0^2}.
```

**Reference.** A. Ametani, “Wave Propagation Characteristics of Cables,”
*IEEE Transactions on Power Apparatus and Systems*, PAS-93, 499–505, 1974;
equation transcription follows Ametani et al., IET, 2021.
"""
function description(::Formula{:Ametani1974})
    "Ametani two-layer magnetic-earth overhead impedance (1974)"
end

function ametanitwolayer1974_gamma(jω, permeability, permittivity)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Ametani1974})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) == 3 || throw(DimensionMismatch(
        ":Ametani1974 requires air and exactly two earth layers"
    ))
    return _stratified_functor(
        Val(:Ametani1974), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

raw"""
Evaluate Ametani's two-layer magnetic-earth overhead impedance:

```math
F_{ij}^{A}=\frac{b_1+b_2+(b_1-b_2)e^{-2a_1d}}
{(\lambda+\mu_0b_1)(b_1+b_2)+
(\lambda-\mu_0b_1)(b_1-b_2)e^{-2a_1d}}e^{-\lambda(h_i+h_j)},
```

where ``b_m=a_m/\mu_m`` and
``a_m=\sqrt{\lambda^2+\gamma_m^2-\gamma_0^2}``. The complete contribution is
``Z_{e,ij}=(j\omega\mu_0/(2\pi))[\ln(D_{ij}/d_{ij})+2\int F_{ij}^A
\cos(y_{ij}\lambda)d\lambda]``.
"""
function ametanitwolayer1974(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    gamma_0_squared = state.gamma_medium_squared[1]
    integral = _quadrature(state) do lambda
        a_1 = sqrt(lambda^2 + state.gamma_medium_squared[2] - gamma_0_squared)
        a_2 = sqrt(lambda^2 + state.gamma_medium_squared[3] - gamma_0_squared)
        b_1 = a_1 / state.mu[2]
        b_2 = a_2 / state.mu[3]
        decay = exp(-2a_1 * d)
        numerator = b_1 + b_2 + (b_1 - b_2) * decay
        denominator = (lambda + state.mu[1] * b_1) * (b_1 + b_2) +
                      (lambda - state.mu[1] * b_1) * (b_1 - b_2) * decay
        numerator / denominator * exp(-lambda * geometry.H) *
        cos(lambda * geometry.y_ij)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

:Ametani1974
