function routes(::Val{:Magalhaes2018})
    (
        self = magalhaesetalz2018,
        mutual = magalhaesetalz2018,
        Γ = magalhaesetalz2018_gamma
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

**Identification.** Wideband homogeneous-earth underground integral.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_1}{2\\pi}\\left[
K_0(\\gamma_1d_{ij})-K_0(\\gamma_1D_{ij})+2\\int_0^\\infty
\\frac{e^{-u_1H}}{u_0+u_1}\\cos(y_{ij}\\lambda)d\\lambda\\right],
\\qquad u_m=\\sqrt{\\lambda^2+\\gamma_m^2}.
```

**Reference.** F. C. R. Magalhães et al., “Closed-Form Expressions for the
Calculation of the Ground-Return Impedance and Admittance of Underground
Cables,” 2018; implemented equation follows the formula corpus and Ametani
et al., IET, 2021.
"""
function description(::Formula{:Magalhaes2018})
    "Magalhaes et al. homogeneous-earth underground impedance (2018)"
end

function magalhaesetalz2018_gamma(jω, permeability, permittivity)
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
function magalhaesetalz2018(functor, pair)
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
