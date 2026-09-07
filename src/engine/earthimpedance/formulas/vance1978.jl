function routes(identifier::Val{:Vance1978})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Vance1978})
    (
        air = _lossless,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Vance1978}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Buried cable in an unbounded conducting medium; self uses cable radius and mutual terms use conductor distance. |
| Calculated quantities | Hankel-ratio approximation for self and mutual earth-return impedance |
| Earth structure | Infinite homogeneous earth. |
| Model and approximation | The result uses the ratio of first-kind Hankel functions at the complex earth propagation-distance product and does not represent the air–earth boundary separately. |
| Main source | E. F. Vance (1978) |
| Citation key(s) | Attributed source: `:Vance1978`; equation witness: `:Guneri2018` |
| Evidence status | Equation checked in the accessible comparative publication; original book attribution retained |

**Expression.**

```math
Z_{e,ii}=\\frac{\\omega\\mu_0}{2\\pi\\gamma_1r_i}
\\frac{H_0^{(1)}(j\\gamma_1r_i)}{H_1^{(1)}(j\\gamma_1r_i)},
\\qquad \\gamma_1^2=j\\omega\\mu_0(\\sigma_1+j\\omega\\varepsilon_1).
```

The comparative source supplies the mutual distance substitution recorded
in the survey. The implementation uses an equivalent scaled modified-Bessel
ratio to avoid underflow at large arguments.

**Reference.** E. F. Vance, *Coupling to Shielded Cables*, Wiley, 1978.
"""
description(::Formula{:Vance1978}) = "Vance infinite-earth radial impedance (1978)"

function propagation_constant(::Val{:Vance1978}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Vance1978})(rho, epsilon, mu, jω, Γ, segments = nothing)
    return _homogeneous_functor(Val(:Vance1978), formula, rho, epsilon, mu, jω, Γ, segments)
end

raw"""
Evaluate Vance's lossy cylindrical-dielectric self impedance:

```math
Z_{e,ii}=\frac{\omega\mu_0}{2\pi\gamma_1r_i}
\frac{H_0^{(1)}(j\gamma_1r_i)}{H_1^{(1)}(j\gamma_1r_i)},\qquad
\gamma_1^2=j\omega\mu_0(\sigma_1+j\omega\varepsilon_1).
```

The cable outer radius is supplied by the self-pair separation field.
"""
function earth_impedance(
        ::Val{:Vance1978}, ::Val{:self}, functor, pair
)
    return earth_impedance(Val(:Vance1978), Val(:mutual), functor, pair)
end

function earth_impedance(
        ::Val{:Vance1978}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    argument = state.gamma[2] * _geometry(pair).d_ij
    πT = one(real(state.jω)) * π
    ratio = special_besselkx(0, argument) / special_besselkx(1, argument)
    return _complex_result(state.jω, state.jω * state.mu[1] / (2πT * argument) * ratio)
end

:Vance1978
