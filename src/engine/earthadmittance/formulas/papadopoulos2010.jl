function routes(identifier::Val{:Papadopoulos2010b})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Papadopoulos2010b})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Papadopoulos2010b}) = Val(:explicit)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel single-core cables represented by their axes for mutual interaction; the outermost radius enters the self substitution. Insulation admittance is separate and is combined through the potential-coefficient matrix in Appendix C. |
| Calculated quantities | Per-unit-length mutual earth-return potential coefficient and admittance correction; self term by the source-prescribed substitution |
| Earth structure | Homogeneous earth half-space ``z\\geq0`` with air in ``z<0`` and planar interface ``z=0``. |
| Model and approximation | The final spectral expression follows the source's lossless longitudinal-propagation approximation ``\\gamma_x\\simeq j\\omega\\sqrt{\\mu_1\\varepsilon_1}`` and the transform identity (4). No analytical truncation of the semi-infinite integral is stated. ``G(\\lambda)`` is retained; dropping it is a different approximation explicitly associated by the source with ``\\gamma_x=\\gamma_1``. |
| Main source | Theofilos A. Papadopoulos, Dimitrios A. Tsiamitros, and Grigoris K. Papagiannis (2010) |
| Citation key(s) | `:Papadopoulos2010b` |
| Evidence status | Verified visually against the original PDF page image |

**Expression.**

```math
P_{e,ij}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
(\\Delta_1+2S_P),
```

```math
\\Delta_1=\\int_0^\\infty
\\frac{e^{-|h_i-h_j|\\alpha_1}-e^{-H\\alpha_1}}{\\alpha_1}
\\cos(y\\lambda)d\\lambda,\\quad
S_P=\\int_0^\\infty
\\frac{e^{-H\\alpha_1}(\\alpha_0+r\\alpha_1)}
{(\\alpha_1+r\\alpha_0)
 (\\alpha_0+r(\\gamma_0^2/\\gamma_1^2)\\alpha_1)}
\\cos(y\\lambda)d\\lambda,\\qquad r=\\mu_1/\\mu_0.
```

**Reference.** T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis,
“Impedances and Admittances of Underground Cables for the Homogeneous Earth
Case,” *IEEE Transactions on Power Delivery*, 25(2), 961–969, 2010.
DOI: 10.1109/TPWRD.2009.2034797.
"""
function description(::Formula{:Papadopoulos2010b})
    "Papadopoulos et al. homogeneous-earth underground potential coefficient (2010)"
end

function propagation_constant(
        ::Val{:Papadopoulos2010b}, jω, permeability, permittivity
)
    squared = oftype(jω, (-jω^2) * permeability * permittivity)
    return (Γ = sqrt(squared), squared)
end

function (formula::Formula{:Papadopoulos2010b})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Papadopoulos2010b), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Papadopoulos et al. underground potential coefficient:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
(\Delta_1+2S_P),
```

```math
\Delta_1=\int_0^\infty\frac{e^{-|h_i-h_j|\alpha_1}-e^{-H\alpha_1}}
{\alpha_1}\cos(y\lambda)d\lambda,\qquad
S_P=\int_0^\infty
\frac{e^{-H\alpha_1}(\alpha_0+r\alpha_1)}
{(\alpha_1+r\alpha_0)
 (\alpha_0+r(\gamma_0^2/\gamma_1^2)\alpha_1)}
\cos(y\lambda)d\lambda,\qquad r=\mu_1/\mu_0,
```

where ``\alpha_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}``.

The permeability-dependent weight is obtained by combining the complete
interface addends in source (5b) and (6c). Its simpler single-denominator
form is valid only when the two permeabilities are equal.
"""
function earth_potential_coefficient(
        ::Val{:Papadopoulos2010b}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    radial_squared = (
        state.gamma_medium_squared[1] + state.gamma_squared,
        state.gamma_medium_squared[2] + state.gamma_squared
    )
    radial = spectral_root(radial_squared[2],state.jω)
    delta_1 = special_besselk(0, radial * geometry.d_ij) -
              special_besselk(0, radial * geometry.D_ij)
    S_P = _quadrature(state) do lambda
        alpha_0 = spectral_root(lambda^2 + radial_squared[1],state.jω)
        alpha_1 = spectral_root(lambda^2 + radial_squared[2],state.jω)
        ratio=state.mu[2]/state.mu[1]
        wave_ratio=state.gamma_medium_squared[1]/state.gamma_medium_squared[2]
        magnetic=alpha_1+ratio*alpha_0
        electric=alpha_0+ratio*wave_ratio*alpha_1
        # Half of the complete F+G interface contribution in (5b),(6c).
        weight=(alpha_0+ratio*alpha_1)/(magnetic*electric)
        exp(-geometry.H * alpha_1) * weight * cos(geometry.y_ij * lambda)
    end
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return _complex_result(state.jω,state.jω /
        (2*(one(geometry.H)*π) * kappa) * (delta_1 + 2 * S_P))
end

:Papadopoulos2010b
