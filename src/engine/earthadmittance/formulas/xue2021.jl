function routes(identifier::Val{:Xue2021})
    (
        self = FormulaMethod(identifier, earth_potential_coefficient, Val(:self)),
        mutual = FormulaMethod(identifier, earth_potential_coefficient, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
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

**Reference.** A. Ametani, H. Xue, T. Ohno, and H. Khalilnezhad,
*Electromagnetic Transients in Large HV Cable Networks: Modeling and
Calculations*, IET, 2021, Section 2.6.3, equation (2.73). The book presents
this closed form as Xue's approximation based on the Saad and Petrache
expressions.
"""
function description(::Formula{:Xue2021})
    "Xue closed-form underground earth potential coefficient (2021)"
end

function propagation_constant(::Val{:Xue2021}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

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
function earth_potential_coefficient(
        ::Val{:Xue2021}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    pair.row == pair.column || _require_horizontal_separation(pair)
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
