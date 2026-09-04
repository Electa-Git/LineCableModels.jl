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
        earth = _conductive,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Vance1978}) = Val(:zero)
"""
$(TYPEDSIGNATURES)

**Identification.** Lossy cylindrical-dielectric radial term for a buried
conductor.

**Expression.**

```math
Z_{e,ii}=\\frac{\\omega\\mu_0}{2\\pi\\gamma_1r_i}
\\frac{H_0^{(1)}(j\\gamma_1r_i)}{H_1^{(1)}(j\\gamma_1r_i)},
\\qquad \\gamma_1^2=j\\omega\\mu_0\\sigma_1.
```

This is a single-cable self term. The cited cylindrical model does not
provide a mutual-impedance specialization for a multiconductor system.

**Reference.** E. F. Vance, *Coupling to Shielded Cables*, Wiley, 1978.
"""
description(::Formula{:Vance1978}) = "Vance lossy-cylinder underground self term (1978)"

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
\gamma_1^2=j\omega\mu_0\sigma_1.
```

The cable outer radius is supplied by the self-pair separation field.
"""
function earth_impedance(
        ::Val{:Vance1978}, ::Val{:self}, functor, pair
)
    _require(pair, Val(:underground))
    state = functor.state
    radius = pair.separation
    argument = im * state.gamma[2] * radius
    omega = imag(state.jω)
    return omega * state.mu[1] / (2π * state.gamma[2] * radius) *
           hankelh1(0, argument) / hankelh1(1, argument)
end

function earth_impedance(
        ::Val{:Vance1978}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:underground))
    throw(ArgumentError(
        "Vance1978 supplies only a single-cable self impedance; " *
        "the cited formula has no mutual route"
    ))
end

:Vance1978
