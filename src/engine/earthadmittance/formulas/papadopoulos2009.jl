function routes(identifier::Val{:Papadopoulos2009})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
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
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Two infinite, electrically thin overhead conductors; the self term uses the source-prescribed radius substitution. Conductor internal and insulation terms are excluded. |
| Calculated quantities | Mutual earth-return potential correction, self correction, homogeneous reduction, and admittance-matrix assembly. |
| Earth structure | Air above a finite earth layer of thickness ``d`` and a lower earth half-space. |
| Model and approximation | Integral representation under the thin-conductor quasi-TEM model. The 2010 derivation replaces the unknown longitudinal constant by the air value and applies a Bessel transform; no additional matrix inversion is applied to (2b). |
| Main source | T. A. Papadopoulos, G. K. Papagiannis, and D. A. Labridis (2009); the 2010 paper supplies the auxiliary derivation. |
| Citation key(s) | Primary: `:Papadopoulos2009`; auxiliary derivation: `:Papadopoulos2010a` |
| Evidence status | Both publications checked against page images. Prime/index conflicts, the root branch, and scalar-inverse interpretation remain unresolved. |

**Expression.**

```math
P_{e,ij}=\\frac1{2\\pi\\varepsilon_0}\\left[\\ln\\frac{D_{ij}}{d_{ij}}+
2\\int_0^\\infty(F_{ij}^{P}+G_{ij}^{P})
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
F_{ij}^{P}=\\mu_1\\frac{s_{12}+d_{12}e^{-2a_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2a_1d}}e^{-\\lambda H},
```

with the electric interface kernel

```math
G_{ij}^{P}=\\lambda
\\frac{\\mu_0\\mu_1(\\gamma_0^2-\\gamma_1^2)
(s_{12}+d_{12}e^{-2a_1d})(S_{12}+D_{12}e^{-2a_1d})-
4\\mu_0\\mu_1^2\\mu_2a_1^2\\gamma_0^2(\\gamma_2^2-\\gamma_1^2)e^{-2a_1d}}
{(S_{01}S_{12}+D_{01}D_{12}e^{-2a_1d})
(s_{01}s_{12}+d_{01}d_{12}e^{-2a_1d})}e^{-\\lambda H}.
```

**Reference.** T. A. Papadopoulos, G. K. Papagiannis, and D. P. Labridis,
“Wave Propagation Characteristics of Overhead Conductors Above Imperfect
Stratified Earth for a Wide Frequency Range,” *IEEE Transactions on
Magnetics*, 45(3), 1064–1067, 2009.
"""
function description(::Formula{:Papadopoulos2009})
    "Papadopoulos et al. two-layer overhead potential coefficient (2009)"
end

function propagation_constant(
        ::Val{:Papadopoulos2009}, jω, permeability, permittivity
)
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
function earth_potential_coefficient(
        ::Val{:Papadopoulos2009}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    # Common permeability and bulk-square scales cancel from F and G.
    # Removing them avoids products below Float32's normal range.
    mu0,mu1,mu2=state.mu ./ state.mu[1]
    gamma0,gamma1,gamma2=state.gamma_medium_squared ./
        maximum(abs,state.gamma_medium_squared)
    integral = _quadrature(state) do lambda
        a0 = spectral_root(lambda^2 + state.gamma_medium_squared[1] + state.gamma_squared,state.jω)
        a1 = spectral_root(lambda^2 + state.gamma_medium_squared[2] + state.gamma_squared,state.jω)
        a2 = spectral_root(lambda^2 + state.gamma_medium_squared[3] + state.gamma_squared,state.jω)
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
           (2*(one(geometry.H)*π) * state.epsilon[1])
end

:Papadopoulos2009
