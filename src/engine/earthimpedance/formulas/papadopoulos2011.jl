function routes(identifier::Val{:Papadopoulos2011})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
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
| Family | External impedance |
| Geometry | Infinite parallel single-core cable axes at positive burial depths ``h_i,h_j`` with horizontal separation ``y_{ij}``. Finite radius enters the self substitution. Internal conductor and insulation contributions are separate. |
| Calculated quantities | Per-unit-length mutual earth-return impedance; self earth term by the published radius/depth substitutions |
| Earth structure | Air above a horizontal finite first earth layer of thickness ``d``; second earth layer is the terminal infinite-depth half-space. This is not an arbitrary-layer or different-layer-source kernel. |
| Model and approximation | The source replaces the unknown longitudinal propagation constant by the dielectric value of the upper earth layer and applies identity (5). Equation (6) remains a spectral integral; no analytical truncation is stated. |
| Main source | T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis, 2011, DOI `10.1049/iet-gtd.2010.0228` |
| Citation key(s) | `:Papadopoulos2011` |
| Evidence status | Original PDF equations checked; the transformed-factor index and root branch remain unresolved. |

**Expression.**

```math
Z'_{e,ij}=\\frac{j\\omega\\mu_1}{2\\pi}\\int_0^\\infty
F_{ij}^{strat}(\\lambda)\\cos(y_{ij}\\lambda)d\\lambda,
```

```math
F_{ij}^{strat}=\\frac{
s_{10}s_{21}e^{-\\alpha_1|h_i-h_j|}+
s_{10}d_{21}e^{-\\alpha_1(2d-h_i-h_j)}-
d_{10}s_{21}e^{-\\alpha_1H}-
d_{10}d_{21}e^{-\\alpha_1(2d-|h_i-h_j|)}}
{\\alpha_1(s_{10}s_{21}+d_{10}d_{21}e^{-2\\alpha_1d})}.
```

**Reference.** T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis,
“Earth Return Admittances and Impedances of Underground Cables in
Non-Homogeneous Earth,” *IET Generation, Transmission & Distribution*, 5(2),
161–171, 2011.
"""
function description(::Formula{:Papadopoulos2011})
    "Papadopoulos et al. two-layer underground mutual impedance (2011)"
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
Evaluate the Papadopoulos et al. two-layer underground impedance:

```math
Z'_{e,ij}=\frac{j\omega\mu_1}{2\pi}\int_0^\infty
F_{ij}^{strat}(\lambda)\cos(y_{ij}\lambda)d\lambda,
```

```math
F_{ij}^{strat}=\frac{
s_{10}s_{21}e^{-\alpha_1|h_i-h_j|}+
s_{10}d_{21}e^{-\alpha_1(2d-h_i-h_j)}-
d_{10}s_{21}e^{-\alpha_1(h_i+h_j)}-
d_{10}d_{21}e^{-\alpha_1(2d-|h_i-h_j|)}}
{\alpha_1(s_{10}s_{21}+d_{10}d_{21}e^{-2\alpha_1d})}.
```

Here ``s_{mn}=\mu_n\alpha_m+\mu_m\alpha_n``, 
``d_{mn}=\mu_m\alpha_n-\mu_n\alpha_m``, and
``\alpha_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}``. Conductors must lie in the
finite top earth layer.
"""
function earth_impedance(
        ::Val{:Papadopoulos2011}, ::Val{:mutual}, functor, pair
)
    pair.layers == (2, 2) || throw(ArgumentError(
        ":Papadopoulos2011 requires both conductors in the top earth layer"
    ))
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    isfinite(d) && d>0 && max(geometry.h_i,geometry.h_j)<=d ||
        throw(DomainError(d,":Papadopoulos2011 requires both depths within a positive finite top-layer thickness"))
    difference = abs(geometry.h_i - geometry.h_j)
    radial = spectral_root(
        state.gamma_medium_squared[2] + state.gamma_squared,state.jω
    )
    direct = special_besselk(0, radial * geometry.d_ij)
    integral = _quadrature(state) do lambda
        alpha0 = spectral_root(lambda^2 + state.gamma_medium_squared[1] + state.gamma_squared,state.jω)
        alpha1 = spectral_root(lambda^2 + state.gamma_medium_squared[2] + state.gamma_squared,state.jω)
        alpha2 = spectral_root(lambda^2 + state.gamma_medium_squared[3] + state.gamma_squared,state.jω)
        s10 = state.mu[1] * alpha1 + state.mu[2] * alpha0
        d10 = state.mu[2] * alpha0 - state.mu[1] * alpha1
        s21 = state.mu[2] * alpha2 + state.mu[3] * alpha1
        d21 = state.mu[3] * alpha1 - state.mu[2] * alpha2
        numerator = s10 * s21 * exp(-alpha1 * difference) +
                    s10 * d21 * exp(-alpha1 * (2d - geometry.H)) -
                    d10 * s21 * exp(-alpha1 * geometry.H) -
                    d10 * d21 * exp(-alpha1 * (2d - difference))
        denominator = alpha1 * (
            s10 * s21 + d10 * d21 * exp(-2alpha1 * d)
        )
        F = numerator / denominator
        direct_spectrum = exp(-alpha1 * difference) / alpha1
        (F - direct_spectrum) * cos(geometry.y_ij * lambda)
    end
    return _complex_result(state.jω,state.jω * state.mu[2] /
        (2*(one(geometry.H)*π)) * (direct + integral))
end

:Papadopoulos2011
