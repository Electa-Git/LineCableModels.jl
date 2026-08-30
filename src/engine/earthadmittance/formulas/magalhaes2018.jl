function routes(::Val{:Magalhaes2018})
    (
        self = magalhaesetal2018,
        mutual = magalhaesetal2018,
        Γ = magalhaesetal2018_gamma
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
function description(::Formula{:Magalhaes2018})
    "Magalhaes et al. underground potential-coefficient kernel (2018)"
end

magalhaesetal2018_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Magalhaes2018})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Magalhaes2018), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate the Magalhaes et al. underground potential coefficient exactly as
printed in the corpus:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[\Lambda_{ij}+2\int_0^\infty I_{ij}^{M}(\lambda)
\cos(y_{ij}\lambda)d\lambda\right],
```

```math
I_{ij}^{M}=\frac{u_0e^{-u_1H}-e^{-u_1H/2}}
{u_0+(\gamma_0^2/\gamma_1^2)u_1},\qquad
\Lambda_{ij}=K_0(\gamma_1d_{ij})-K_0(\gamma_1D_{ij}).
```

The source record rates this extracted kernel as low confidence; no
dimensional repair is applied here.
"""
function magalhaesetal2018(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_0_squared / gamma_1_squared
    integral = _quadrature(state) do lambda
        u_0 = sqrt(lambda^2 + gamma_0_squared)
        u_1 = sqrt(lambda^2 + gamma_1_squared)
        numerator = u_0 * exp(-u_1 * geometry.H) - exp(-u_1 * geometry.H / 2)
        numerator / (u_0 + ratio * u_1) * cos(geometry.y_ij * lambda)
    end
    gamma_1 = state.gamma[2]
    direct = special_besselk(0, gamma_1 * geometry.d_ij) -
             special_besselk(0, gamma_1 * geometry.D_ij)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) * (direct + 2 * integral)
end

:Magalhaes2018
