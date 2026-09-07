function routes(identifier::Val{:Papadopoulos2011})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Papadopoulos2011})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2011}) = Val(:explicit)
media(::Formula{:Papadopoulos2011}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel single-core cable axes; outer radius enters self substitution. The earth coefficient is separate from dielectric-insulation and conductor contributions; the source references their assembly to [9,15]. |
| Calculated quantities | Mutual earth-return potential coefficient and source-printed admittance relation; self earth term by radius/depth substitutions |
| Earth structure | Air half-space, finite upper earth layer of thickness ``d``, infinite-depth lower earth half-space. General different-layer/mixed/arbitrary-layer arrangements are not supplied by this final kernel. |
| Model and approximation | The source replaces the unknown longitudinal propagation constant by the dielectric value of the upper earth layer and applies transform (5). ``G`` is retained; no spectral truncation or tail approximation is specified. |
| Main source | T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis (2011), DOI `10.1049/iet-gtd.2010.0228` |
| Citation key(s) | `:Papadopoulos2011` |
| Evidence status | Original PDF equations checked; the transformed-factor index, root branch, and scalar/matrix interpretation remain unresolved. |

**Expression.**

```math
P_{e,ij}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\int_0^\\infty[F_{ij}^{strat}+G_{ij}^{strat}]
\\cos(y_{ij}\\lambda)d\\lambda,
```

```math
G_{ij}^{strat}=2\\alpha_1(G_{1,ij}+G_{2,ij}+G_{3,ij}+G_{4,ij}),
```

where ``F_{ij}^{strat}`` is the four-exponential impedance kernel and

```math
A_{mn}=\\alpha_n\\gamma_m^2\\mu_n+\\alpha_m\\gamma_n^2\\mu_m,\\qquad
\\Delta_{mn}=\\alpha_n\\gamma_m^2\\mu_n-\\alpha_m\\gamma_n^2\\mu_m.
```

**Reference.** T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis,
“Earth Return Admittances and Impedances of Underground Cables in
Non-Homogeneous Earth,” *IET Generation, Transmission & Distribution*, 5(2),
161–171, 2011.
"""
function description(::Formula{:Papadopoulos2011})
    "Papadopoulos et al. two-layer underground potential coefficient (2011)"
end

function propagation_constant(
        ::Val{:Papadopoulos2011}, jω, permeability, permittivity
)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2011})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) == 3 || throw(DimensionMismatch(
        ":Papadopoulos2011 requires air and exactly two earth layers"
    ))
    return _stratified_functor(
        Val(:Papadopoulos2011), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

raw"""
Evaluate the Papadopoulos et al. two-layer underground potential coefficient:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\int_0^\infty[F_{ij}^{strat}+G_{ij}^{strat}]
\cos(y_{ij}\lambda)d\lambda,
```

where ``F_{ij}^{strat}`` is the four-exponential impedance kernel and
``G_{ij}^{strat}=2\alpha_1(G_{1,ij}+G_{2,ij}+G_{3,ij}+G_{4,ij})``. The four
interface terms, and the definitions

```math
A_{mn}=\alpha_n\gamma_m^2\mu_n+\alpha_m\gamma_n^2\mu_m,\qquad
\Delta_{mn}=\alpha_n\gamma_m^2\mu_n-\alpha_m\gamma_n^2\mu_m,
```

are implemented without algebraic contraction so the corpus symbols remain
traceable in code. Conductors must lie in the finite top earth layer.
"""
function earth_potential_coefficient(
        ::Val{:Papadopoulos2011}, ::Val{:mutual}, functor, pair
)
    pair.layers == (2, 2) || throw(ArgumentError(
        ":Papadopoulos2011 requires both conductors in the top earth layer"
    ))
    state = functor.state
    geometry = _geometry(pair)
    hi, hj = geometry.h_i, geometry.h_j
    difference = abs(hi - hj)
    d = state.thickness[2]
    isfinite(d) && d>0 && max(geometry.h_i,geometry.h_j)<=d ||
        throw(DomainError(d,":Papadopoulos2011 requires both depths within a positive finite top-layer thickness"))
    radial = spectral_root(
        state.gamma_medium_squared[2] + state.gamma_squared,state.jω
    )
    direct = special_besselk(0, radial * geometry.d_ij)
    # Common permeability and bulk-square scales cancel from F and G.
    # Removing them avoids products below Float32's normal range.
    mu0,mu1,mu2=state.mu ./ state.mu[1]
    gamma0,gamma1,gamma2=state.gamma_medium_squared ./
        maximum(abs,state.gamma_medium_squared)
    integral = _quadrature(state) do lambda
        alpha0 = spectral_root(lambda^2 + state.gamma_medium_squared[1] + state.gamma_squared,state.jω)
        alpha1 = spectral_root(lambda^2 + state.gamma_medium_squared[2] + state.gamma_squared,state.jω)
        alpha2 = spectral_root(lambda^2 + state.gamma_medium_squared[3] + state.gamma_squared,state.jω)
        s10 = mu0 * alpha1 + mu1 * alpha0
        d10 = mu1 * alpha0 - mu0 * alpha1
        s21 = mu1 * alpha2 + mu2 * alpha1
        d21 = mu2 * alpha1 - mu1 * alpha2
        A10 = alpha0 * gamma1 * mu0 + alpha1 * gamma0 * mu1
        Delta10 = alpha0 * gamma1 * mu0 - alpha1 * gamma0 * mu1
        A12 = alpha2 * gamma1 * mu2 + alpha1 * gamma2 * mu1
        Delta12 = alpha2 * gamma1 * mu2 - alpha1 * gamma2 * mu1
        decay = exp(-2alpha1 * d)
        magnetic = s10 * s21 + d10 * d21 * decay
        electric = A10 * A12 - Delta10 * Delta12 * decay
        common = electric * magnetic
        F = (
            s10 * s21 * exp(-alpha1 * difference) +
            s10 * d21 * exp(-alpha1 * (2d - hi - hj)) -
            d10 * s21 * exp(-alpha1 * (hi + hj)) -
            d10 * d21 * exp(-alpha1 * (2d - difference))
        ) / (alpha1 * magnetic)
        G1 = mu1 * mu2 * (gamma1 - gamma2) *
             (
                 s10 * A10 * exp(-alpha1 * (2d - hi - hj)) -
                 d10 * A10 * exp(-alpha1 * (2d + hi - hj))
             ) / common
        G2 = mu1 * mu2 * (gamma1 - gamma2) *
             (
                 s10 * Delta10 * exp(-alpha1 * (2d + hj - hi)) -
                 d10 * Delta10 * exp(-alpha1 * (2d + hi + hj))
             ) / common
        G3 = mu1 * mu0 * (gamma1 - gamma0) *
             (
                 s21 * Delta12 * exp(-alpha1 * (2d + hi - hj)) +
                 d21 * Delta12 * exp(-alpha1 * (4d - hi - hj))
             ) / common
        G4 = mu1 * mu0 * (gamma1 - gamma0) *
             (
                 s21 * A12 * exp(-alpha1 * (hi + hj)) +
                 d21 * A12 * exp(-alpha1 * (2d + hj - hi))
             ) / common
        direct_spectrum = exp(-alpha1 * difference) / alpha1
        (F - direct_spectrum + 2alpha1 * (G1 + G2 + G3 + G4)) *
        cos(geometry.y_ij * lambda)
    end
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return _complex_result(state.jω,state.jω /
        (2*(one(geometry.H)*π) * kappa) * (direct + integral))
end

:Papadopoulos2011
