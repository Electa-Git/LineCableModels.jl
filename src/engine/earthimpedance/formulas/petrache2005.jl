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

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried cable in an unbounded conducting medium; self uses cable radius and mutual terms use conductor distance. |
| Calculated quantities | Logarithmic approximation for self and mutual earth-return impedance |
| Earth structure | Infinite homogeneous earth. |
| Model and approximation | A closed logarithm in the complex earth propagation-distance product approximates the infinite-earth return term. |
| Main source | Marian Petrache, Florin Rachidi, Mario Paolone, Carlo Alberto Nucci, Vladimir A. Rakov, and M. A. Uman (2005) |
| Citation key(s) | Primary attribution: `:Petrache2005`; equation witness: `:Guneri2018` |
| Evidence status | Equation checked in the accessible comparative publication; publication identity and DOI verified |

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
    state = functor.state
    argument = state.gamma[2] * _geometry(pair).d_ij
    return state.jω * state.mu[1] / (2π) * log((1 + argument) / argument)
end

:Petrache2005
