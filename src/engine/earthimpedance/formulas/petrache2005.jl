function routes(identifier::Val{:Petrache2005})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Petrache2005})
    (
        air = _lossless,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Petrache2005}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Burial-depth-independent logarithmic underground
approximation.

**Expression.**

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\ln\\left(\\frac{1+\\gamma_1R_{ab}}{\\gamma_1R_{ab}}\\right).
```

**Reference.** E. Petrache, F. Rachidi, M. Paolone, C. A. Nucci, V. A.
Rakov, and M. A. Uman, “Lightning Induced Disturbances in Buried Cables—Part
I: Theory,” *IEEE Transactions on Electromagnetic Compatibility*, 47(3),
498–508, 2005. DOI: 10.1109/TEMC.2005.853161.
"""
function description(::Formula{:Petrache2005})
    "Petrache logarithmic underground approximation (2005)"
end

function propagation_constant(::Val{:Petrache2005}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Petrache2005})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(
        Val(:Petrache2005), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate the burial-depth-independent Petrache et al. approximation:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}
\ln\left(\frac{1+\gamma_1R_{ab}}{\gamma_1R_{ab}}\right),\qquad
\gamma_1^2=j\omega\mu_0(\sigma_1+j\omega\varepsilon_1).
```
"""
function earth_impedance(
        ::Val{:Petrache2005}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    pair.row == pair.column || _require_horizontal_separation(pair)
    state = functor.state
    argument = state.gamma[2] * pair.separation
    return state.jω * state.mu[1] / (2π) * log((1 + argument) / argument)
end

:Petrache2005
