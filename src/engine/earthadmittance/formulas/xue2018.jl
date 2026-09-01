function routes(identifier::Val{:Xue2018})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        infinite = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:infinite_depth)
        ),
        surface = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:surface_reference)
        ),
        penetration = FormulaMethod(
            identifier, earth_potential_coefficient, Val(:penetration_depth)
        ),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Xue2018})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Xue2018}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Generalized underground potential coefficient. The
registered expression is referenced to infinite earth depth; surface and
penetration-depth corrections are retained as alternative physical forms.

**Expression.**

```math
P_{e,ij}^{\\infty}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\left[K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2S_{12}^c+
2\\gamma_1^2S_{13}^c\\right],
```

```math
S_{12}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\lambda^2\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)[u_0+(\\gamma_0^2/\\gamma_1^2)u_1]}d\\lambda,
\\quad
S_{13}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda.
```

The surface reference is
``P_{e,ij}^{0}=P_{e,ij}^{\\infty}-P_{ec,ij}``; the penetration-depth reference
uses the same structure at the exact lossy-medium penetration depth.

**Reference.** H. Xue, *Electromagnetic Transients in Large HV Cable
Networks*, doctoral thesis, Delft University of Technology, 2018; equations
as consolidated in Ametani et al., IET, 2021.
"""
function description(::Formula{:Xue2018})
    "Xue et al. generalized underground potential coefficient (2018)"
end

function propagation_constant(::Val{:Xue2018}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Xue2018})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Xue2018), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Xue et al. underground potential coefficient referred to
infinite earth depth:

```math
P_{e,ij}^{\infty}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[K_0(\gamma_1d_{ij})-K_0(\gamma_1D_{ij})+2S_{12}^c+
2\gamma_1^2S_{13}^c\right],
```

```math
S_{12}^c=\int_0^\infty\frac{e^{-Hu_1}\lambda^2\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)[u_0+(\gamma_0^2/\gamma_1^2)u_1]}d\lambda,
\quad
S_{13}^c=\int_0^\infty\frac{e^{-Hu_1}\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)(u_0+u_1)}d\lambda.
```

This infinite-depth reference is the registered default. The surface and
penetration-depth references are support leaves of the same `:Xue2018`
entry, available through [`routes`](@ref) without creating more formula IDs.
"""
function integral_terms(
        ::Val{:Xue2018}, state, geometry, height = geometry.H
)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    S12 = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-height * u_1) * lambda^2 * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + ratio * u_1))
    end
    S13 = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-height * u_1) * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
    end
    return S12, S13
end

function earth_potential_coefficient(
        identifier::Val{:Xue2018}, ::Val{:mutual}, functor, pair
)
    return earth_potential_coefficient(
        identifier, Val(:infinite_depth), functor, pair
    )
end

function earth_potential_coefficient(
        identifier::Val{:Xue2018}, ::Val{:infinite_depth}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    S12, S13 = integral_terms(identifier, state, geometry)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) *
           (direct + 2 * S12 + 2 * gamma_1^2 * S13)
end

raw"""
Evaluate the surface-referenced support coefficient from Xue et al.:

```math
P_{e,ij}^{0}=P_{e,ij}^{\infty}-P_{ec,ij},\qquad
P_{ec,ij}=\frac{j\omega}{\pi(\sigma_1+j\omega\varepsilon_1)}
(S_{14}^c+\gamma_1^2S_{15}^c),
```

where ``S_{14}^c`` and ``S_{15}^c`` use the infinite-depth denominators and
the attenuation ``e^{-\frac12(h_i+h_j)u_1}``.
"""
function earth_potential_coefficient(
        identifier::Val{:Xue2018}, ::Val{:surface_reference}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    S12, S13 = integral_terms(identifier, state, geometry)
    S14, S15 = integral_terms(identifier, state, geometry, geometry.H / 2)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    infinite = state.jω / (2π * kappa) *
               (direct + 2 * S12 + 2 * gamma_1^2 * S13)
    correction = state.jω / (π * kappa) * (S14 + gamma_1^2 * S15)
    return infinite - correction
end

raw"""
Evaluate the penetration-depth-referenced support coefficient from Xue et al.:

```math
P_{e,ij}^{\delta}=P_{e,ij}^{\infty}-P_{\delta_e,ij},
```

```math
h_r=-\left\{\frac{\omega^2\varepsilon_1\mu_0}{2}
\left[\sqrt{1+(\sigma_1/(\omega\varepsilon_1))^2}-1\right]
\right\}^{-1/2}.
```

The correction retains the corpus distances ``d_{\delta,ij}``,
``D_{\delta,ij}`` and integrals ``S_{16}^c``, ``S_{17}^c`` verbatim.
"""
function earth_potential_coefficient(
        identifier::Val{:Xue2018}, ::Val{:penetration_depth}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    S12, S13 = integral_terms(identifier, state, geometry)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    omega = imag(state.jω)
    bracket = sqrt(1 + (state.sigma[2] / (omega * state.epsilon[2]))^2) - 1
    h_r = -inv(sqrt(omega^2 * state.epsilon[2] * state.mu[1] / 2 * bracket))
    exponent_height = geometry.H / 2 - h_r
    S16 = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-exponent_height * u_1) * lambda^2 * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + ratio * u_1))
    end
    S17 = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-exponent_height * u_1) * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
    end
    d_delta = hypot(geometry.y_ij, h_r - geometry.H / 2)
    D_delta = hypot(geometry.y_ij, h_r + geometry.H / 2)
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    reference = special_besselk(0, gamma_1 * d_delta) -
                special_besselk(0, gamma_1 * D_delta) +
                2 * S16 + 2 * gamma_1_squared * S17
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    infinite = direct + 2 * S12 + 2 * gamma_1_squared * S13
    return state.jω / (2π * kappa) * (infinite - reference)
end

:Xue2018
