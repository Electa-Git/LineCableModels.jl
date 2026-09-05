function routes(identifier::Val{:Magalhaes2018})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Magalhaes2018})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Magalhaes2018}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Companion homogeneous-earth underground impedance used
with the Magalhães earth-admittance formulation. Under the registered
homogeneous assumptions this kernel is algebraically equivalent to the
Xue2018 generalized impedance route.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_1}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2\\int_0^\\infty
\\frac{e^{-u_1H}}{u_0+u_1}\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\qquad u_m=\\sqrt{\\lambda^2+\\gamma_m^2}.
```

**Reference.** A. P. C. Magalhães, M. T. C. de Barros, and A. C. S. Lima,
“Earth Return Admittance Effect on Underground Cable System Modeling,” *IEEE
Transactions on Power Delivery*, 33(2), 662–670, 2018. The companion
impedance identity is documented with the generalized underground formulas in
A. Ametani et al., *Electromagnetic Transients in Large HV Cable Networks*,
IET, 2021.
"""
function description(::Formula{:Magalhaes2018})
    "Magalhaes companion Xue-equivalent homogeneous-earth impedance (2018)"
end

function propagation_constant(
        ::Val{:Magalhaes2018}, jω, permeability, permittivity
)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Magalhaes2018})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Magalhaes2018), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Magalhaes et al. homogeneous-earth underground kernel:

```math
Z_{e,ij}=\frac{j\omega\mu_1}{2\pi}\left[K_0(\gamma_1d_{ij})-
K_0(\gamma_1D_{ij})+2\int_0^\infty
\frac{e^{-u_1(h_i+h_j)}}{u_0+u_1}\cos(y_{ij}\lambda)d\lambda\right],
```

where ``u_m=\sqrt{\lambda^2+\gamma_m^2}``.
"""
function earth_impedance(
        ::Val{:Magalhaes2018}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    integral = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        exp(-u_1 * geometry.H) / (u_0 + u_1) *
        cos(geometry.y_ij * lambda)
    end
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    return state.jω * state.mu[2] / (2π) * (direct + 2 * integral)
end

:Magalhaes2018
