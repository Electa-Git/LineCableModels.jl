function routes(identifier::Val{:Xue2018})
    (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
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

**Identification.** Generalized homogeneous-earth underground wideband
impedance.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2S_{11}^c+
2\\gamma_1^2S_{13}^c\\right],
```

```math
S_{11}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\lambda^2\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda,\\qquad
S_{13}^c=\\int_0^\\infty\\frac{e^{-Hu_1}\\cos(y\\lambda)}
{(\\lambda^2+\\gamma_1^2)(u_0+u_1)}d\\lambda.
```

**Reference.** H. Xue, *Electromagnetic Transients in Large HV Cable
Networks*, doctoral thesis, Delft University of Technology, 2018; equations
as consolidated in Ametani et al., IET, 2021.
"""
function description(::Formula{:Xue2018})
    "Xue et al. generalized homogeneous-earth underground impedance (2018)"
end

function propagation_constant(::Val{:Xue2018}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Xue2018})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Xue2018), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate the Xue et al. generalized homogeneous-earth underground impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\gamma_1d_{ij})-K_0(\gamma_1D_{ij})+2S_{11}^c+
2\gamma_1^2S_{13}^c\right],
```

```math
S_{11}^c=\int_0^\infty\frac{e^{-Hu_1}\lambda^2\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)(u_0+u_1)}d\lambda,\qquad
S_{13}^c=\int_0^\infty\frac{e^{-Hu_1}\cos(y\lambda)}
{(\lambda^2+\gamma_1^2)(u_0+u_1)}d\lambda,
```

where ``u_m=\sqrt{\lambda^2+\gamma_m^2}``.
"""
function earth_impedance(
        ::Val{:Xue2018}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    S11 = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-geometry.H * u_1) * lambda^2 * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
    end
    S13 = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-geometry.H * u_1) * cos(geometry.y_ij * lambda) /
        ((lambda^2 + gamma_1_squared) * (u_0 + u_1))
    end
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    return state.jω * state.mu[1] / (2π) *
           (direct + 2 * S11 + 2 * gamma_1_squared * S13)
end

:Xue2018
