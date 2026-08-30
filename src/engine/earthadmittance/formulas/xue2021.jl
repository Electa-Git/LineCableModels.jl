function routes(::Val{:Xue2021})
    (
        self = xueclosedform2021,
        mutual = xueclosedform2021,
        Γ = xueclosedform2021_gamma
    )
end

function assumptions(::Val{:Xue2021})
    (
        air = _full,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Xue2021}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Closed-form underground potential coefficient combining
a logarithmic radial term with an interface correction.

**Expression.**

```math
P_{e,ij}=\\frac{j\\omega}{2\\pi(\\sigma_1+j\\omega\\varepsilon_1)}
\\left[\\ln\\left(\\frac{1+\\gamma_1R_{ab}}{\\gamma_1R_{ab}}\\right)+
\\frac{2e^{-H\\gamma_1}}{4+\\gamma_1^2R_{ab}^2}\\right].
```

**Reference.** H. Xue et al., “Generalized Formulation and Surge Analysis on
Overhead Lines: Impedance/Admittance of a Multi-Layer Earth,” *IEEE
Transactions on Power Delivery*, 36(6), 3834–3845, 2021.
DOI: 10.1109/TPWRD.2021.3049595.
"""
function description(::Formula{:Xue2021})
    "Xue closed-form underground earth potential coefficient (2021)"
end

xueclosedform2021_gamma(jω, permeability, permittivity) = (Γ = zero(jω), squared = zero(jω))

function (formula::Formula{:Xue2021})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    return _homogeneous_functor(
        Val(:Xue2021), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

raw"""
Evaluate Xue's closed-form underground potential coefficient:

```math
P_{e,ij}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\left[\ln\left(\frac{1+\gamma_1R_{ab}}{\gamma_1R_{ab}}\right)+
\frac{2e^{-(h_i+h_j)\gamma_1}}{4+\gamma_1^2R_{ab}^2}\right].
```
"""
function xueclosedform2021(functor, pair)
    _require(pair, Val(:underground))
    state = functor.state
    geometry = _geometry(pair)
    gamma = state.gamma[2]
    argument = gamma * geometry.y_ij
    bracket = log((1 + argument) / argument) +
              2exp(-geometry.H * gamma) / (4 + argument^2)
    kappa = state.sigma[2] + state.jω * state.epsilon[2]
    return state.jω / (2π * kappa) * bracket
end

:Xue2021
