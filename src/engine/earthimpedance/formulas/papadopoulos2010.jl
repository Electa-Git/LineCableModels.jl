function routes(identifier::Val{:Papadopoulos2010b})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
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
| Family | External impedance |
| Geometry | Infinite parallel single-core cables represented by their axes for mutual interaction; depths ``h_i,h_j`` and horizontal separation ``y_{ij}``; the outermost cable radius is used in the stated self substitution. Internal conductor and insulation terms are assembled separately. |
| Calculated quantities | Per-unit-length mutual earth-return impedance correction; self term by the source-prescribed substitution |
| Earth structure | Homogeneous earth half-space ``z\\geq0`` with air in ``z<0`` and a planar interface at ``z=0``. |
| Model and approximation | Equation (5) applies the lossless longitudinal-propagation approximation ``\\gamma_x\\simeq j\\omega\\sqrt{\\mu_1\\varepsilon_1}`` to the dipole-field starting solution and uses ``u^2-k_x^2=\\lambda^2`` with identity (4). No series truncation is applied; the semi-infinite integral is evaluated numerically. |
| Main source | Theofilos A. Papadopoulos, Dimitrios A. Tsiamitros, and Grigoris K. Papagiannis (2010) |
| Citation key(s) | `:Papadopoulos2010b` |
| Evidence status | Verified visually against the original PDF page image |

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_1}{2\\pi}(\\Delta_1+2S),
```

```math
\\Delta_1=\\int_0^\\infty
\\frac{e^{-|h_i-h_j|\\alpha_1}-e^{-H\\alpha_1}}{\\alpha_1}
\\cos(y_{ij}\\lambda)d\\lambda,\\quad
S=\\int_0^\\infty\\frac{e^{-H\\alpha_1}}{\\alpha_1+(\\mu_1/\\mu_0)\\alpha_0}
\\cos(y_{ij}\\lambda)d\\lambda,
```

where ``\\alpha_m=\\sqrt{\\lambda^2+\\gamma_m^2+k_x^2}``.

**Reference.** T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis,
“Impedances and Admittances of Underground Cables for the Homogeneous Earth
Case,” *IEEE Transactions on Power Delivery*, 25(2), 961–969, 2010.
DOI: 10.1109/TPWRD.2009.2034797.
"""
function description(::Formula{:Papadopoulos2010b})
    "Papadopoulos et al. homogeneous-earth underground impedance (2010)"
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
Evaluate the Papadopoulos et al. homogeneous-earth underground impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_1}{2\pi}(\Delta_1+2S),
```

```math
\Delta_1=\int_0^\infty\frac{e^{-|h_i-h_j|\alpha_1}-e^{-(h_i+h_j)\alpha_1}}
{\alpha_1}\cos(y_{ij}\lambda)d\lambda,\qquad
S=\int_0^\infty\frac{e^{-(h_i+h_j)\alpha_1}}
{\alpha_1+(\mu_1/\mu_0)\alpha_0}\cos(y_{ij}\lambda)d\lambda,
```

where ``\alpha_m=\sqrt{\lambda^2+\gamma_m^2+k_x^2}`` and
``k_x=\omega\sqrt{\mu_1\varepsilon_1}`` by default. An explicit `Γ`
replaces ``k_x`` without changing the leaf routes.
"""
function earth_impedance(
        ::Val{:Papadopoulos2010b}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    alpha_squared = (
        state.gamma_medium_squared[1] + state.gamma_squared,
        state.gamma_medium_squared[2] + state.gamma_squared
    )
    radial = spectral_root(alpha_squared[2],state.jω)
    delta_1 = special_besselk(0, radial * geometry.d_ij) -
              special_besselk(0, radial * geometry.D_ij)
    S = _quadrature(state) do lambda
        alpha_0 = spectral_root(lambda^2 + alpha_squared[1],state.jω)
        alpha_1 = spectral_root(lambda^2 + alpha_squared[2],state.jω)
        exp(-geometry.H * alpha_1) /
        (alpha_1 + (state.mu[2]/state.mu[1])*alpha_0) * cos(geometry.y_ij * lambda)
    end
    return _complex_result(state.jω,state.jω * state.mu[2] /
        (2*(one(geometry.H)*π)) * (delta_1 + 2 * S))
end

:Papadopoulos2010b
