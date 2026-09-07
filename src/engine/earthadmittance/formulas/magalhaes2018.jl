function routes(identifier::Val{:Magalhaes2018})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Magalhaes2018})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Magalhaes2018}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Insulated circular cables represented by their axes; the outer insulation radius enters the self term. Each core/sheath pair receives identical earth-return entries. |
| Calculated quantities | Single-cable and multiconductor per-unit-length earth-return admittance matrices; self and mutual entries of the intermediate matrix ``\\Lambda-T``; source-provided series combination with insulation admittance |
| Earth structure | Homogeneous medium 1 beneath homogeneous medium 2, separated by a plane interface. |
| Model and approximation | Integral and matrix representation under the source's quasi-TEM reduction of the cited full-wave model. The paper does not reproduce the parent longitudinal-propagation prescription. |
| Main source | A. P. C. Magalhães et al. (2018). |
| Citation key(s) | `:Magalhaes2018` |
| Evidence status | Original publication and accepted-manuscript page images checked. |

**Expression.**

```math
P_{e,ij}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\left[\\Lambda_{ij}-2\\int_0^\\infty I_{ij}^{M}(\\lambda)
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
I_{ij}^{M}=\\frac{u_0}{u_1}
\\frac{e^{-u_1H/2}-e^{-u_1H}}
{u_0+(\\gamma_0^2/\\gamma_1^2)u_1},\\qquad
\\Lambda_{ij}=K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij}).
```

**Reference.** A. P. C. Magalhães, M. T. C. de Barros, and A. C. S. Lima,
“Earth Return Admittance Effect on Underground Cable System Modeling,” *IEEE
Transactions on Power Delivery*, 33(2), 662–670, 2018; as reproduced in A.
Ametani, H. Xue, T. Ohno, and H. Khalilnezhad, *Electromagnetic Transients in
Large HV Cable Networks*, Section 2.6.2.2, equations (2.70)--(2.71).
"""
function description(::Formula{:Magalhaes2018})
    "Magalhaes earth-surface-reference potential coefficient (2018)"
end

function propagation_constant(
        ::Val{:Magalhaes2018}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Magalhaes2018})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Magalhaes2018), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Magalhaes et al. underground potential coefficient:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[\Lambda_{ij}-2\int_0^\infty I_{ij}^{M}(\lambda)
\cos(y_{ij}\lambda)d\lambda\right],
```

```math
I_{ij}^{M}=\frac{u_0}{u_1}
\frac{e^{-u_1H/2}-e^{-u_1H}}
{u_0+(\gamma_0^2/\gamma_1^2)u_1},\qquad
\Lambda_{ij}=K_0(\gamma_1d_{ij})-K_0(\gamma_1D_{ij}).
```

For the same material assumptions this expression is algebraically equivalent
to the Xue et al. earth-surface-reference formulation.
"""
function earth_potential_coefficient(
        ::Val{:Magalhaes2018}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    integral = _quadrature(state) do lambda
        u_0 = spectral_root(lambda^2 + gamma_0_squared,state.jω)
        u_1 = spectral_root(lambda^2 + gamma_1_squared,state.jω)
        attenuation = exp(-u_1 * geometry.H / 2) * (-expm1(-u_1 * geometry.H / 2))
        u_0 / u_1 * attenuation / (u_0 + ratio * u_1) *
        cos(geometry.y_ij * lambda)
    end
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return _complex_result(state.jω,state.jω /
        (2*(one(geometry.H)*π) * kappa) * (direct - 2 * integral))
end

:Magalhaes2018
